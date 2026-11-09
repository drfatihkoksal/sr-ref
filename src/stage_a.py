"""Stage A cohort (analysis plan 5.1): HEEDB records with PVCs, eligibility, and matched pairs.

Candidates: records whose 12SL v24 statements report PVCs and none of the excluded rhythms (the
rules of `src.audit`), from patients with at least two such records. Pass 1 reads the headers and
keeps 500-Hz records of patients with at least two of them. Pass 2 runs `records.process` (locked
detector v2, lead-quality rule) on each and keeps, for eligible records with an OT-like PVC, the
8-lead PVC and SR templates (float32, -100 to +600 ms around the QRS onset).

Pairs: per patient, eligible OT-like recordings in time order; two recordings match when the
multi-lead correlation of their 8-lead PVC windows (-10 to +180 ms from the QRS onset, per-lead
centring, lag up to ±10 ms) is >= 0.85. The pair is the first matched pair in time order within
each interval stratum: same day (< 24 h), primary (<= 365 days), any interval. Patient and
acquisition time (de-identified, intervals preserved within a patient) come from the file names on
MGH and from the metadata on EUH, whose file names carry neither.

Row-level outputs go to ARTIFACTS/wp2 (never committed); the report holds counts only:

    python -m src.stage_a cohort --site mgh|euh   -> docs/wp2_cohort_<site>.md
"""
import argparse
import os
import re
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import wfdb

from src.audit import EXC_12SL, INC_12SL
from src.beats import MS
from src.config import ARTIFACTS, DOCS, HEEDB
from src.records import load, process
from src.signal import IDX8

SITES = {"mgh": HEEDB / "ECG" / "I0001", "euh": HEEDB / "ECG" / "I0006"}
OUT = ARTIFACTS / "wp2"
MATCH_R = 0.85
WIN = (int(10 * MS), int(180 * MS))       # PVC window around the QRS onset
LAG = int(10 * MS)
STRATA = {"same day": 1.0, "primary (<= 365 days)": 365.0, "any interval": np.inf}


def candidates(site):
    root = SITES[site]
    d = pd.read_csv(root / "12SL_diagnoses" / "diagnoses_dictionary.csv")
    inc = set(d.codes[d.diagnoses.str.contains(INC_12SL, regex=True)])
    exc = set(d.codes[d.diagnoses.str.contains(EXC_12SL, regex=True)])
    rows = []
    for ch in pd.read_csv(root / "12SL_diagnoses" / "diagnoses_v24.csv", chunksize=500_000, dtype=str):
        ch = ch[ch.codes.str.contains(r"\b23[12]\b", regex=True, na=False)]
        for f, c in zip(ch.FileName, ch.codes):
            codes = {int(v) for v in re.findall(r"\d+", c)}
            if codes & inc and not codes & exc:
                rows.append(f.strip().lstrip("."))
    df = pd.DataFrame({"path": rows})
    if site == "mgh":       # MGH file names carry the patient and the acquisition time
        parts = df.path.str.extract(r"de_(\d+)_(\d{14})_")
        df["patient"], df["time"] = parts[0], pd.to_datetime(parts[1], format="%Y%m%d%H%M%S", errors="coerce")
        df["file"] = [str(root / "WFDB" / p.lstrip("/").removesuffix(".hea")) for p in df.path]
    else:                   # EUH file names do not: patient and acquisition time from the metadata
        m = pd.read_csv(root / "metadata" / "metadata.csv", usecols=["BDSPPatientID", "FileID", "ECGAcquisitionTime"],
                        dtype=str, engine="pyarrow").set_index("FileID")
        fid = df.path.str.rsplit("/", n=1).str[-1].str.removesuffix(".hea")
        df["patient"] = m.BDSPPatientID.reindex(fid).to_numpy()
        df["time"] = pd.to_datetime(m.ECGAcquisitionTime.reindex(fid).to_numpy(), errors="coerce")
        df["file"] = [str(root / p.lstrip("/").removesuffix(".hea")) for p in df.path]
        df = df[df.patient.notna()]
    n = df.groupby("patient").path.transform("size")
    return df[(n >= 2) & df.time.notna()].reset_index(drop=True)


def _fs(f):
    try:
        return wfdb.rdheader(f).fs
    except Exception:
        return np.nan


def _process(f):
    try:
        x, fs = load(f)
    except Exception:
        return {"reason": "unreadable"}, None
    if x is None or np.isnan(x).any():
        return {"reason": "missing lead or samples"}, None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = process(x, fs, "heedb")
    row = {"reason": r["reason"] or ("OT-like" if r["ot_like"] else "eligible, not OT-like"),
           "n_pvc": r["pvc"]["n"] if r.get("pvc") else 0, "n_sr": r["sr"]["n"] if r.get("sr") else 0}
    if not (r["eligible"] and r["ot_like"]):
        return row, None
    return row, {"pvc": r["pvc"]["T"][IDX8].astype(np.float32), "sr": r["sr"]["T"][IDX8].astype(np.float32),
                 "pvc_on": r["pvc"]["on"], "pvc_off": r["pvc"]["off"], "sr_on": r["sr"]["on"], "sr_off": r["sr"]["off"]}


def match_r(a, b, on_a, on_b):
    """Multi-lead correlation of two 8-lead PVC windows, best lag within ±LAG."""
    wa = a[:, on_a - WIN[0]:on_a + WIN[1] + 1]
    wa = wa - wa.mean(1, keepdims=True)
    best = -1.0
    for s in range(-LAG, LAG + 1):
        lo, hi = on_b + s - WIN[0], on_b + s + WIN[1] + 1
        if lo < 0 or hi > b.shape[1]:
            continue
        wb = b[:, lo:hi] - b[:, lo:hi].mean(1, keepdims=True)
        best = max(best, float((wa * wb).sum() / np.sqrt((wa * wa).sum() * (wb * wb).sum() + 1e-12)))
    return best


def pairs(rec, tmpl):
    """First matched pair per patient and stratum: list of dicts."""
    out = []
    ok = rec[rec.reason == "OT-like"].sort_values(["patient", "time"])
    for pid, g in ok.groupby("patient"):
        if len(g) < 2:
            continue
        idx, times = g.index.to_list(), g.time.to_list()
        found = {}
        for i in range(len(idx)):
            for j in range(i + 1, len(idx)):
                days = (times[j] - times[i]).total_seconds() / 86400
                todo = [s for s, lim in STRATA.items() if s not in found and days <= lim]
                if not todo:
                    continue
                a, b = tmpl[idx[i]], tmpl[idx[j]]
                r = match_r(a["pvc"], b["pvc"], a["pvc_on"], b["pvc_on"])
                if r >= MATCH_R:
                    for s in todo:
                        found[s] = {"patient": pid, "stratum": s, "i": idx[i], "j": idx[j], "days": days, "r": r}
        out += list(found.values())
    return pd.DataFrame(out)


def cohort(site):
    OUT.mkdir(parents=True, exist_ok=True)
    c = candidates(site)
    with ProcessPoolExecutor(max_workers=max(os.cpu_count() - 4, 1)) as ex:
        c["fs"] = list(ex.map(_fs, c.file, chunksize=500))
    c500 = c[c.fs == 500]
    c500 = c500[c500.groupby("patient").path.transform("size") >= 2].reset_index(drop=True)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    with ProcessPoolExecutor(max_workers=max(os.cpu_count() - 4, 1)) as ex:
        res = list(ex.map(_process, c500.file, chunksize=50))
    rec = pd.concat([c500[["patient", "time", "file"]], pd.DataFrame([r for r, _ in res])], axis=1)
    tmpl = {i: t for i, (_, t) in enumerate(res) if t is not None}
    rec.to_parquet(OUT / f"{site}_records.parquet")
    np.save(OUT / f"{site}_templates.npy", np.array([tmpl.get(i) for i in range(len(rec))], dtype=object),
            allow_pickle=True)
    pr = pairs(rec, tmpl)
    pr.to_parquet(OUT / f"{site}_pairs.parquet")
    ot = rec[rec.reason == "OT-like"].groupby("patient").size()
    lines = [f"# WP2: Stage A cohort, HEEDB-{site.upper()}", "",
             f"Generated by `python -m src.stage_a cohort --site {site}`. Counts only; row-level outputs in",
             "ARTIFACTS/wp2. Definitions: module docstring and analysis plan 5.1.", "",
             f"- Records with PVC statements (no excluded statement) from patients with ≥ 2 such records: "
             f"{len(c)} ({c.patient.nunique()} patients).",
             f"- Sampling rate: 500 Hz {int((c.fs == 500).sum())}, 250 Hz {int((c.fs == 250).sum())}, other or "
             f"unreadable {int((~c.fs.isin([250, 500])).sum())}.",
             f"- 500-Hz records of patients with ≥ 2 of them: {len(c500)} ({c500.patient.nunique()} patients); "
             "all processed.", "", "| Outcome | Records |", "|---|---|"]
    for k, n in rec.reason.value_counts().items():
        lines.append(f"| {k} | {n} |")
    lines += ["", f"- Patients with ≥ 2 eligible OT-like records: {int((ot >= 2).sum())} "
              f"({int(ot[ot >= 2].sum())} records).", "",
              f"Matched pairs (PVC windows r ≥ {MATCH_R}), one per patient and stratum:", "",
              "| Stratum | Patients | Interval, median [IQR] days | r, median [IQR] |", "|---|---|---|---|"]
    for s in STRATA:
        g = pr[pr.stratum == s] if len(pr) else pr
        if len(g):
            q = lambda v, d: f"{np.median(v):.{d}f} [{np.percentile(v, 25):.{d}f}, {np.percentile(v, 75):.{d}f}]"
            lines.append(f"| {s} | {len(g)} | {q(g.days, 1)} | {q(g.r, 2)} |")
        else:
            lines.append(f"| {s} | 0 | — | — |")
    (DOCS / f"wp2_cohort_{site}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["cohort"])
    ap.add_argument("--site", choices=["mgh", "euh"], default="mgh")   # EUH: the confirmatory H4 set (plan 3.3)
    a = ap.parse_args()
    cohort(a.site)


if __name__ == "__main__":
    main()
