"""Ningbo IVA pipeline (analysis plan 4.1-4.5, 4.9; 5.6 item 6). Labels are not read.

Per record (raw CSV, 2000 Hz, undeclared integer unit): leads mapped by header name; the first
0.4 s cut (the files rise from zero over up to ~0.33 s; RECON-ORIGIN); pacing stimuli detected at
2000 Hz before any filtering (RECON-ORIGIN `signals.stimuli`, ported unchanged).

Beat classes (docs/decisions.md item 22), `--source`:
- `adjudicated` (primary): the PVC and SR beats of the RECON-ORIGIN adjudication (blinded to label,
  sub-site and sex; two readers, κ 0.97): beat positions from its clusters table, PVC and sinus
  clusters from its assignments. Records the adjudication excluded (G0b, rejected) or marked as
  AF/flutter in SR are excluded. Paced beats (RECON-ORIGIN flag) and their neighbours (±1 beat) are
  dropped. Templates: `records.from_beats` (SR-REF template code and eligibility).
- `auto` (sensitivity): the SR-REF detector (`records.process`) with the dominant PVC cluster and
  the SR fallback (decision 22), and the stimuli above.

Scale (plan 4.1, 5.6 item 6): pass 1 runs on the raw integers; s = median SR-loop RMS of HEEDB-MGH
(objects/sr_loop_rms.json) / median SR-loop RMS of the eligible Ningbo records; pass 2 runs on the
signals multiplied by s (µV-equivalent). The pipeline is scale-free except the flat-lead bound;
records whose outcome differs between the passes are counted.

Agreement with RECON-ORIGIN (plan 4.9, label-free): per record in both, the SR-REF PVC (and SR)
template, low-passed at 40 Hz, against the RECON-ORIGIN adjudicated template (0.5-40 Hz; its -100 to
+200 ms around the peak), best multi-lead correlation over all alignments of the 12 leads; r < 0.90 flags the record for the
visual QC.

Row-level outputs (ARTIFACTS/wp1/ningbo, never committed): records.parquet, templates.npy (12-lead
PVC and SR templates of eligible records), strips.npy (II and V1 at 250 Hz with beat marks, for the
QC pages), scale.json. Report with counts only:

    python -m src.ningbo [--source adjudicated|auto]   -> docs/wp1_ningbo_<source>.md
"""
import argparse
import json
import os
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks, iirnotch, sosfiltfilt

from src.config import ARTIFACTS, DOCS, IVA_RAW, LEADS12, RECON_ORIGIN_ARTIFACTS
from src.h4 import OBJ, sr_loop
from src.records import from_beats, process
from src.signal import IDX8

FS_RAW = 2000
TRIM = int(0.4 * FS_RAW)
OUT = ARTIFACTS / "wp1" / "ningbo"
RO = RECON_ORIGIN_ARTIFACTS / "wp1"
RO_TRIM_S = 0.4                            # RECON-ORIGIN beat positions include its 0.4-s trim


def adjudication():
    """record -> PVC fiducials, SR fiducials, excluded fiducials (500 Hz, in the trimmed signal), or a
    reason the record is excluded."""
    a = pd.read_parquet(RO / "assignments.parquet").set_index("record")
    c = pd.read_parquet(RO / "clusters.parquet")
    out = {}
    for rid, g in c.groupby("record"):
        r = a.loc[rid]
        if r.status != "adjudicated":
            out[str(rid)] = f"adjudication: {r.status}"
            continue
        if bool(r.non_sinus):
            out[str(rid)] = "SR from AF/flutter (adjudication)"
            continue
        g = g.sort_values("pos_s")
        fid = np.round((g.pos_s.to_numpy() - RO_TRIM_S) * 500).astype(int)
        cl, paced = g.cluster.to_numpy(), g.paced.to_numpy()
        near = paced | np.r_[paced[1:], False] | np.r_[False, paced[:-1]]
        parse = lambda v: {int(x) for x in str(v).split(";") if x}
        out[str(rid)] = {"pvc": fid[np.isin(cl, list(parse(r.pvc_clusters)))].tolist(),
                         "sr": fid[np.isin(cl, list(parse(r.sinus_clusters)))].tolist(),
                         "excluded": fid[near].tolist(), "paced": int(paced.sum())}
    return out
AGREE_R = 0.90
LP40 = butter(3, 40.0, "lowpass", fs=500, output="sos")

# ---- pacing stimuli (RECON-ORIGIN signals.py, unchanged) -------------------------------------
SPIKE_HP = butter(4, 80, btype="high", fs=FS_RAW, output="sos")
SPIKE_LF = butter(3, [5, 40], btype="band", fs=FS_RAW, output="sos")
NOTCHES = [iirnotch(f0, 30, FS_RAW) for f0 in (50, 100, 150)]
SPIKE_Z, SPIKE_LEAD_Z, SPIKE_MIN_LEADS = 8, 10, 6
SPIKE_MAX_MS, SPIKE_QUIET, SPIKE_HF_RATIO = 6, 3, 0.05


def stimuli(x):
    """Sample indices (2000 Hz) of pacing-stimulus artefacts; x is samples x leads."""
    y = x - x.mean(0)
    for b, a in NOTCHES:
        y = filtfilt(b, a, y, axis=0)
    h = sosfiltfilt(SPIKE_HP, y, axis=0)
    lf = sosfiltfilt(SPIKE_LF, y, axis=0)
    z = np.abs(h) / (np.median(np.abs(h), 0) + 1e-9)
    e = np.median(z, 1)
    pk, _ = find_peaks(e, height=SPIKE_Z, distance=int(0.1 * FS_RAW))
    out = []
    for p in pk:
        if p < 400 or p > len(e) - 400:
            continue
        lo = hi = p
        while lo > 0 and e[lo] > e[p] / 2:
            lo -= 1
        while hi < len(e) - 1 and e[hi] > e[p] / 2:
            hi += 1
        w = slice(p - 30, p + 30)
        ratio = np.median((h[w] ** 2).sum(0) / ((lf[w] ** 2).sum(0) + 1e-9))
        if ((hi - lo) / FS_RAW * 1000 <= SPIKE_MAX_MS and (z[p] > SPIKE_LEAD_Z).sum() >= SPIKE_MIN_LEADS
                and np.median(e[p - 120:p - 40]) < SPIKE_QUIET and ratio >= SPIKE_HF_RATIO):
            out.append(p)
    return np.array(out, dtype=int)


def load(path):
    """Raw integers (12 x samples, LEADS12 order by header name) after the 0.4-s trim."""
    d = pd.read_csv(path)
    return d[LEADS12].to_numpy(dtype=float)[TRIM:].T


def halve(x):
    n = x.shape[-1] // 2 * 2
    return x[..., :n].reshape(*x.shape[:-1], n // 2, 2).mean(-1)


def run_record(args):
    path, scale, source, adj = args
    x = load(path)
    stim = stimuli(x.T)
    row = {"record": path.stem, "seconds": x.shape[1] / FS_RAW, "n_stim": len(stim)}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if source == "auto":
            r = process(x * scale, FS_RAW, "ningbo", stim=stim // 4, pvc_cluster="dominant", sr_fallback=True)
            b = r.get("beats")
            marks = [(int(t) // 2, c) for t, c in zip(b.fid, b.cls) if c in ("PVC", "SR")] if b is not None else []
            row["n_paced"] = int(b.paced.sum()) if b is not None and "paced" in b else 0
        elif isinstance(adj, str):
            return {**row, "reason": adj, "ot_like": False, "n_pvc": 0, "n_sr": 0, "n_paced": 0}, None, None
        else:
            r = from_beats(x * scale, FS_RAW, "ningbo", adj["pvc"], adj["sr"], adj["excluded"])
            ex = set(adj["excluded"])
            marks = [(f // 2, "PVC") for f in adj["pvc"] if f not in ex] + [(f // 2, "SR") for f in adj["sr"] if f not in ex]
            row["n_paced"] = adj["paced"]
    row.update({"reason": r["reason"] or "eligible", "ot_like": bool(r["ot_like"]),
                "n_pvc": r["pvc"]["n"] if r.get("pvc") else 0, "n_sr": r["sr"]["n"] if r.get("sr") else 0})
    tmpl = strip = None
    if r["eligible"]:
        tmpl = {k: r[k]["T"].astype(np.float32) for k in ("pvc", "sr")}
        tmpl.update({f"{k}_{e}": r[k][e] for k in ("pvc", "sr") for e in ("on", "off")})
        x12 = r["x12"]
        strip = {"strip": halve(x12[[LEADS12.index("II"), LEADS12.index("V1")]]).astype(np.float32),
                 "marks": sorted(marks)}
    return row, tmpl, strip


def run_all(scale, source):
    files = sorted(IVA_RAW.glob("*.csv"))
    adj = adjudication() if source == "adjudicated" else {}
    args = [(f, scale, source, adj.get(f.stem, "not in the adjudication")) for f in files]
    with ProcessPoolExecutor(max_workers=max(os.cpu_count() - 4, 1)) as ex:
        res = list(ex.map(run_record, args, chunksize=4))
    return pd.DataFrame([r for r, _, _ in res]), [t for _, t, _ in res], [s for _, _, s in res]


def loop_rms(tmpl):
    return np.array([np.sqrt((sr_loop({"sr": t["sr"][IDX8], "sr_on": t["sr_on"], "sr_off": t["sr_off"]}) ** 2).mean())
                     for t in tmpl if t is not None])


RO_CROP = slice(50, 200)                   # RECON-ORIGIN window -100 to +200 ms around its peak


def best_r(a, b):
    """Best multi-lead correlation of b slid along a (12 x n). b is the RECON-ORIGIN template cut to
    -100..+200 ms around its peak: its full ±200 ms window starts earlier than the SR-REF template
    (-100 ms from the QRS onset), so the QRS-aligned position would not fit."""
    b = b[:, RO_CROP]
    bb = b - b.mean(1, keepdims=True)
    best = -1.0
    for s in range(a.shape[1] - b.shape[1] + 1):
        w = a[:, s:s + b.shape[1]]
        w = w - w.mean(1, keepdims=True)
        best = max(best, float((w * bb).sum() / np.sqrt((w * w).sum() * (bb * bb).sum() + 1e-12)))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=["adjudicated", "auto"], default="adjudicated")
    source = ap.parse_args().source
    out = OUT / source
    out.mkdir(parents=True, exist_ok=True)
    rec1, tmpl1, _ = run_all(1.0, source)
    mgh = json.loads((OBJ / "sr_loop_rms.json").read_text())["median_uv"]
    scale = mgh / float(np.median(loop_rms(tmpl1)))
    rec, tmpl, strips = run_all(scale, source)
    (out / "scale.json").write_text(json.dumps({"scale_uv_per_unit": scale, "mgh_median_uv": mgh}))
    rec.to_parquet(out / "records.parquet")
    np.save(out / "templates.npy", np.array(tmpl, dtype=object), allow_pickle=True)
    np.save(out / "strips.npy", np.array(strips, dtype=object), allow_pickle=True)
    changed = int((rec.reason != rec1.reason).sum())
    # agreement with RECON-ORIGIN adjudicated templates
    ro = np.load(RECON_ORIGIN_ARTIFACTS / "wp1" / "templates.npz")
    ro_idx = {str(r): i for i, r in enumerate(ro["record"])}
    agree = []
    for (_, row), t in zip(rec.iterrows(), tmpl):
        if t is None or row.record not in ro_idx:
            continue
        k = ro_idx[row.record]
        agree.append({"record": row.record,
                      "r_pvc": best_r(sosfiltfilt(LP40, t["pvc"], axis=-1), ro["pvc_filt"][k]),
                      "r_sr": best_r(sosfiltfilt(LP40, t["sr"], axis=-1), ro["sinus_filt"][k])})
    agree = pd.DataFrame(agree)
    agree.to_parquet(out / "agreement.parquet")
    q = lambda v: f"median {np.median(v):.3f} (IQR {np.percentile(v, 25):.3f} to {np.percentile(v, 75):.3f})"
    lines = [f"# WP1: Ningbo IVA pipeline, beat classes: {source} (label-free)", "",
             f"Generated by `python -m src.ningbo --source {source}`; definitions in the module docstring. Labels",
             f"are not read; row-level outputs in ARTIFACTS/wp1/ningbo/{source}.", "",
             f"- Records: {len(rec)}; duration after the 0.4-s trim median {rec.seconds.median():.1f} s "
             f"(range {rec.seconds.min():.1f} to {rec.seconds.max():.1f}).",
             f"- Pacing stimuli found in {int((rec.n_stim > 0).sum())} records; paced beats {int(rec.n_paced.sum())} "
             f"in {int((rec.n_paced > 0).sum())} records (excluded with their neighbours).",
             f"- Scale factor s = {scale:.4f} µV per raw unit (HEEDB-MGH median SR-loop RMS {mgh:.0f} µV / Ningbo "
             f"median in raw units); outcome changed between pass 1 and pass 2 in {changed} records.", "",
             "| Outcome | Records |", "|---|---|"]
    for k, n in rec.reason.value_counts().items():
        lines.append(f"| {k} | {n} |")
    elig = rec[rec.reason == "eligible"]
    lines += ["", f"- Eligible: {len(elig)} of {len(rec)} ({len(elig) / len(rec):.1%}); OT-like PVC template in "
              f"{int(elig.ot_like.sum())}. PVC beats per template median {elig.n_pvc.median():.0f}, SR beats "
              f"{elig.n_sr.median():.0f}.", "",
              "## Agreement with the RECON-ORIGIN adjudicated templates", "",
              f"- Eligible records also in RECON-ORIGIN: {len(agree)}.",
              f"- PVC template r: {q(agree.r_pvc)}; r ≥ {AGREE_R}: {int((agree.r_pvc >= AGREE_R).sum())}.",
              f"- SR template r: {q(agree.r_sr)}; r ≥ {AGREE_R}: {int((agree.r_sr >= AGREE_R).sum())}.",
              f"- Records flagged for the visual QC (either r < {AGREE_R}): "
              f"{int(((agree.r_pvc < AGREE_R) | (agree.r_sr < AGREE_R)).sum())}."]
    lines += ["", f"- K3 (≥ 95% of 334 with a valid PVC and SR pair): {len(elig)}/334 = {len(elig) / 334:.1%} "
              f"-> {'met' if len(elig) / 334 >= 0.95 else 'not met'} (before the visual check of flagged records)."]
    (DOCS / f"wp1_ningbo_{source}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
