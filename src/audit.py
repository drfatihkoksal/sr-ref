"""Template audit, gate K1-T (docs/decisions.md item 11; analysis plan 4.8).

Sample: 75 eligible records with an OT-like PVC from HEEDB-MGH (I0001, 500 Hz) and 75 from
MIMIC-IV-ECG. Candidates are records whose machine statements report PVCs and none of the
excluded rhythms (pacing, AF/flutter, complete bundle branch block, pre-excitation, ventricular
or supraventricular tachycardia, idioventricular or ventricular rhythm); one record per patient.
Candidates are shuffled with `SEED` and processed in that order (`records.process`, locked
rules) until 75 eligible OT-like records per source are found, which equals a random sample of
eligible OT-like records within the statement-filtered pool.

Outputs (ARTIFACTS/wp1/audit, never committed):
- key.csv (case code -> dataset and record path; never shown to a reader), cases.npy;
- reader1.html (Fatih Köksal) and reader2.html (Fatih Levent): self-contained pages, each with its
  own case order, showing per case the 10-s strip (II, V1) with the classified beats and the 12-lead
  PVC and SR templates, without dataset, record, age or sex. Verdicts stay in the browser and are
  exported as <reader>_form.csv.
The pages hold HEEDB waveforms (BDSP DUA): they stay on this machine and are read only by readers
with BDSP HEEDB access. MIMIC-IV-ECG itself is open access (ODbL 1.0).

docs/wp1_audit_sample.md: counts only.

Scoring: both forms (<reader>_form.csv in OUT) -> `disagreements` lists the codes the readers
judged differently (disagreements.csv, to be settled jointly in consensus.csv with the same columns
as the forms); `score` takes agreed verdicts plus consensus.csv, and writes
docs/wp1_audit_result.md: PVC- and SR-template PPV with Wilson 95% CIs (gate: both ≥ 0.95),
percent agreement and Cohen κ per question (reported, not gated).

    python -m src.audit sample         # selects the cases, writes key.csv, cases.npy and the report
    python -m src.audit pages          # builds the two reader pages from cases.npy
    python -m src.audit disagreements  # after both forms are in OUT
    python -m src.audit score          # after consensus.csv
"""
import argparse
import base64
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from src.config import ARTIFACTS, DOCS, HEEDB, LEADS12, MIMIC, SEED
from src.records import load, process

N_PER_SOURCE = 75
OUT = ARTIFACTS / "wp1" / "audit"
PAGE = Path(__file__).with_name("audit.html")
READERS = {"reader1": "Fatih Köksal", "reader2": "Fatih Levent"}
HEEDB_MGH = HEEDB / "ECG" / "I0001"
INC_12SL = r"^premature ventricular (?:complexes|and fusion complexes)$"
EXC_12SL = (r"fibrillation|flutter|pac(?:ed|ing|emaker)|pre-excitation|wolff|^(?:right|left) bundle branch block|"
            r"bundle branch block -or-|idioventricular|ventricular tachycardia|ventricular rhythm|ventricular escape")
INC_MIMIC = r"\bpvc|ventricular premature complex|premature ventricular"
EXC_MIMIC = (r"pac(?:ed|ing|emaker|er)|fibrillation|flutter|(?<!incomplete )(?:right|left) bundle branch block|"
             r"pre-excitation|wpw|wolff|idioventricular|ventricular tachycardia|ventricular rhythm|aberran|multiform")


def heedb_candidates():
    d = pd.read_csv(HEEDB_MGH / "12SL_diagnoses" / "diagnoses_dictionary.csv")
    inc = set(d.codes[d.diagnoses.str.contains(INC_12SL, regex=True)])
    exc = set(d.codes[d.diagnoses.str.contains(EXC_12SL, regex=True)])
    rows = []
    for ch in pd.read_csv(HEEDB_MGH / "12SL_diagnoses" / "diagnoses_v24.csv", chunksize=500_000, dtype=str):
        ch = ch[ch.codes.str.contains(r"\b23[12]\b", regex=True, na=False)]
        for f, c in zip(ch.FileName, ch.codes):
            codes = {int(v) for v in re.findall(r"\d+", c)}
            if codes & inc and not codes & exc:
                rows.append(f.strip().lstrip("."))
    df = pd.DataFrame({"path": rows})
    df["patient"] = df.path.str.extract(r"de_(\d+)_")[0]
    df["file"] = [str(HEEDB_MGH / "WFDB" / p.lstrip("/").removesuffix(".hea")) for p in df.path]
    return df.sample(frac=1, random_state=SEED).drop_duplicates("patient"), inc, exc


def mimic_candidates():
    cols = [f"report_{i}" for i in range(18)]
    mm = pd.read_csv(MIMIC / "machine_measurements.csv", usecols=["subject_id", "study_id"] + cols, dtype=str)
    rep = mm[cols].fillna("").agg(" ".join, axis=1).str.lower()
    ok = rep.str.contains(INC_MIMIC, regex=True) & ~rep.str.contains(EXC_MIMIC, regex=True)
    rl = pd.read_csv(MIMIC / "record_list.csv", dtype=str)[["study_id", "path"]]
    df = mm.loc[ok, ["subject_id", "study_id"]].merge(rl, on="study_id")
    df["file"] = [str(MIMIC / p) for p in df.path]
    df["patient"] = df.subject_id
    return df.sample(frac=1, random_state=SEED).drop_duplicates("patient")


def halve(x):
    n = x.shape[-1] // 2 * 2
    return x[..., :n].reshape(*x.shape[:-1], n // 2, 2).mean(-1)


def select(cands, dataset, n_target):
    """Process candidates in order until n_target eligible OT-like records; reason counts."""
    cases, reasons, n_seen = [], {}, 0
    for f in cands.file:
        if len(cases) >= n_target:
            break
        n_seen += 1
        try:
            x, fs = load(f)
        except Exception:
            reasons["unreadable"] = reasons.get("unreadable", 0) + 1
            continue
        if x is None or np.isnan(x).any():
            reasons["missing lead or samples"] = reasons.get("missing lead or samples", 0) + 1
            continue
        if dataset == "heedb" and fs != 500:
            reasons["not 500 Hz"] = reasons.get("not 500 Hz", 0) + 1
            continue
        r = process(x, fs, dataset)
        key = r["reason"] or ("OT-like" if r["ot_like"] else "eligible, not OT-like")
        reasons[key] = reasons.get(key, 0) + 1
        if r["eligible"] and r["ot_like"]:
            b = r["beats"]
            cases.append({"dataset": dataset, "file": f,
                          "strip": halve(r["x12"][[LEADS12.index("II"), LEADS12.index("V1")]]),
                          "marks": [(int(t) // 2, c) for t, c in zip(b.fid, b.cls) if c in ("PVC", "SR")],
                          "pvc": halve(r["pvc"]["T"]), "sr": halve(r["sr"]["T"]),
                          "pvc_qrs": (r["pvc"]["on"] // 2, r["pvc"]["off"] // 2),
                          "sr_qrs": (r["sr"]["on"] // 2, r["sr"]["off"] // 2)})
    return cases, reasons, n_seen


def sample():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "key.csv").exists():
        raise SystemExit("key.csv exists: the sample is fixed; delete it only with the author's approval")
    hc, inc, exc = heedb_candidates()
    mc = mimic_candidates()
    h_cases, h_reasons, h_seen = select(hc, "heedb", N_PER_SOURCE)
    m_cases, m_reasons, m_seen = select(mc, "mimic", N_PER_SOURCE)
    cases = h_cases + m_cases
    order = np.random.default_rng(SEED).permutation(len(cases))
    codes = [f"K{i + 1:03d}" for i in range(len(cases))]
    cases = [cases[i] for i in order]
    for c, code in zip(cases, codes):
        c["code"] = code
    pd.DataFrame({"code": codes, "dataset": [c["dataset"] for c in cases], "file": [c["file"] for c in cases]}
                 ).to_csv(OUT / "key.csv", index=False)
    np.save(OUT / "cases.npy", np.array([{k: v for k, v in c.items() if k not in ("dataset", "file")}
                                         for c in cases], dtype=object), allow_pickle=True)
    lines = ["# WP1: template audit sample (gate K1-T)", "",
             "Generated by `python -m src.audit sample`. Counts only; the key is outside the repository.", "",
             f"- HEEDB-MGH candidates (12SL v24 codes {sorted(inc)}, none of {len(exc)} excluded codes, one per "
             f"patient): {len(hc)}; processed {h_seen}; selected {len(h_cases)}.",
             f"- MIMIC-IV-ECG candidates (machine statements, one per subject): {len(mc)}; processed {m_seen}; "
             f"selected {len(m_cases)}.", "", "| Outcome | HEEDB-MGH | MIMIC-IV-ECG |", "|---|---|---|"]
    for k in sorted(set(h_reasons) | set(m_reasons)):
        lines.append(f"| {k} | {h_reasons.get(k, 0)} | {m_reasons.get(k, 0)} |")
    (DOCS / "wp1_audit_sample.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def b64(a, dtype):
    return base64.b64encode(np.ascontiguousarray(a, dtype=dtype).tobytes()).decode()


def pages():
    cases = list(np.load(OUT / "cases.npy", allow_pickle=True))
    payload = []
    for c in cases:
        top_t = float(np.abs(np.r_[c["pvc"].ravel(), c["sr"].ravel()]).max()) or 1.0
        top_s = float(np.abs(c["strip"]).max()) or 1.0
        payload.append({"code": c["code"], "fs": 250, "leads": LEADS12,
                        "strip": b64(np.round(c["strip"] / top_s * 127), "<i1"), "strip_scale": top_s / 127,
                        "strip_n": int(c["strip"].shape[1]), "marks": c["marks"],
                        "pvc": b64(np.round(c["pvc"] / top_t * 32000), "<i2"),
                        "sr": b64(np.round(c["sr"] / top_t * 32000), "<i2"),
                        "t_scale": top_t / 32000, "t_n": int(c["pvc"].shape[1]),
                        "pvc_qrs": [int(v) for v in c["pvc_qrs"]], "sr_qrs": [int(v) for v in c["sr_qrs"]]})
    html = PAGE.read_text(encoding="utf-8")
    sample_id = hashlib.sha256((OUT / "key.csv").read_bytes()).hexdigest()[:12]
    by_code = {p["code"]: p for p in payload}
    # reader orders are kept across rebuilds (orders.json): codes no longer in the sample drop out,
    # new codes go to the end in a random order
    orders_file = OUT / "orders.json"
    orders = json.loads(orders_file.read_text()) if orders_file.exists() else {}
    for i, (reader, name) in enumerate(READERS.items(), start=1):
        if reader in orders:
            keep = [c for c in orders[reader] if c in by_code]
            new = [c for c in by_code if c not in orders[reader]]
            orders[reader] = keep + [new[j] for j in np.random.default_rng(SEED + 100 + i).permutation(len(new))]
        else:
            codes = list(by_code)
            orders[reader] = [codes[j] for j in np.random.default_rng(SEED + i).permutation(len(codes))]
        data = {"reader": reader, "name": name, "sample_id": sample_id,
                "cases": [by_code[c] for c in orders[reader]]}
        (OUT / f"{reader}.html").write_text(html.replace("/*__DATA__*/null", json.dumps(data)), encoding="utf-8")
    orders_file.write_text(json.dumps(orders))
    print(f"pages: {', '.join(str(OUT / f'{r}.html') for r in READERS)}")


QUESTIONS = {"pvc_ventricular": "PVC template is a ventricular ectopic beat",
             "sr_sinus": "SR template is a conducted sinus beat"}


def forms():
    f = [pd.read_csv(OUT / f"{r}_form.csv", dtype=str).set_index("code") for r in READERS]
    for r, d in zip(READERS, f):
        missing = d[list(QUESTIONS)].isna().any(axis=1).sum()
        if missing:
            raise SystemExit(f"{r}: {missing} cases without both answers")
    return f


def disagreements():
    a, b = forms()
    diff = a.index[(a[list(QUESTIONS)] != b.loc[a.index, list(QUESTIONS)]).any(axis=1)]
    pd.DataFrame({"code": diff}).to_csv(OUT / "disagreements.csv", index=False)
    print(f"{len(diff)} cases to settle by consensus: {OUT / 'disagreements.csv'}")


def wilson(k, n, z=1.959964):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def kappa(x, y):
    po = np.mean(x == y)
    pe = sum(np.mean(x == v) * np.mean(y == v) for v in ("yes", "no"))
    return po, (po - pe) / (1 - pe) if pe < 1 else np.nan


def score():
    a, b = forms()
    cons = pd.read_csv(OUT / "consensus.csv", dtype=str).set_index("code") if (OUT / "consensus.csv").exists() else None
    lines = ["# WP1: template audit result (gate K1-T)", "",
             "Generated by `python -m src.audit score`. Readers: " + ", ".join(READERS.values()) + ".", ""]
    passed = True
    for q, text in QUESTIONS.items():
        x, y = a[q], b.loc[a.index, q]
        final = x.where(x == y)
        todo = final.index[final.isna()]
        if len(todo):
            if cons is None or not set(todo) <= set(cons.index):
                raise SystemExit(f"{q}: {len(todo)} disagreements without consensus")
            final.loc[todo] = cons.loc[todo, q]
        k, n = int((final == "yes").sum()), len(final)
        lo, hi = wilson(k, n)
        po, kap = kappa(x.to_numpy(), y.to_numpy())
        passed &= k / n >= 0.95
        lines += [f"- **{text}:** {k}/{n} = {k / n:.3f} (Wilson 95% CI {lo:.3f} to {hi:.3f}); "
                  f"agreement {po:.1%}, Cohen κ {kap:.2f}; {len(todo)} settled by consensus."]
    key = pd.read_csv(OUT / "key.csv").set_index("code")
    lines += ["", "By source (consensus verdicts):", "", "| Source | PVC template ventricular | SR template sinus |",
              "|---|---|---|"]
    for src, g in key.groupby("dataset"):
        row = []
        for q in QUESTIONS:
            v = a[q].where(a[q] == b.loc[a.index, q])
            if cons is not None:
                v = v.fillna(cons[q].reindex(v.index))
            v = v.loc[g.index]
            row.append(f"{int((v == 'yes').sum())}/{len(v)}")
        lines.append(f"| {src} | {row[0]} | {row[1]} |")
    lines.append("")
    for r, d in zip(READERS.values(), (a, b)):
        t = pd.to_datetime(d.updated, utc=True).sort_values()
        lines.append(f"- {r}: answers given over {(t.max() - t.min()).total_seconds() / 60:.0f} min, "
                     f"median {t.diff().dt.total_seconds().median():.0f} s per case.")
    lines += ["", f"**K1-T {'passed' if passed else 'not passed'}** (threshold 0.95 for both, point estimates)."]
    (DOCS / "wp1_audit_result.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["sample", "pages", "disagreements", "score"])
    {"sample": sample, "pages": pages, "disagreements": disagreements, "score": score}[ap.parse_args().step]()


if __name__ == "__main__":
    main()
