"""INCART detector evaluation (analysis plan 4.8, gate K1).

Records are cut into 14-s windows with a 10-s stride. Each window is pre-processed and
classified on its own, as a short clinical record would be, and only beats in its central 10 s
are scored, so every scored beat has a preceding beat in its window. Detections are matched to
the reference annotations within 150 ms (ANSI/AAMI EC57).

- QRS detection: sensitivity and PPV over all reference beats.
- V class: reference V = {V, E}. A detection classified as PVC that matches a fusion (F) or
  unclassifiable (Q) beat is neither a true nor a false positive (EC57).
- 95% CIs from 2000 record-level bootstrap resamples (gross statistics).

Tuning subset (any number of runs; writes docs/wp1_incart_tuning.md, aggregates only):
    python -m src.incart tuning
Validation subset (once after the rules are locked at K1; one re-validation after a revision):
    python -m src.incart validation --confirm-locked      (detector: config.DETECTOR)
"""
import argparse
import json
import os
import re
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import wfdb

from src.beats import METHODS, P_CORR_MIN, classify, failed_criteria
from src.config import (ARTIFACTS, DETECTOR, DOCS, EC57_MATCH_MS, FS, INCART, INCART_FS, INCART_N_TUNING, LEADS12,
                        N_BOOT, RULES, SEED)
from src.signal import lead8, preprocess

WIN_S, STRIDE_S, GUARD_S = 14, 10, 2
BEATS = set("NLRBAaJSVrFejnE/fQ")
V_REF = {"V", "E"}
NEUTRAL = {"F", "Q"}
OUT = ARTIFACTS / "wp1" / "incart"
VALIDATION_LOG = OUT / "validation_runs.json"
MAX_VALIDATION_RUNS = 2


def records():
    """Record names and their patient numbers (header comment '# patient N')."""
    rows = []
    for hea in sorted(INCART.glob("I*.hea")):
        txt = hea.read_text()
        m = re.search(r"#\s*patient\s+(\d+)", txt)
        rows.append({"record": hea.stem, "patient": int(m.group(1))})
    return pd.DataFrame(rows)


def split():
    """Patient-level split with SEED: INCART_N_TUNING tuning patients, the rest validation."""
    recs = records()
    patients = np.sort(recs.patient.unique())
    tuning = set(np.random.default_rng(SEED).permutation(patients)[:INCART_N_TUNING].tolist())
    recs["subset"] = np.where(recs.patient.isin(tuning), "tuning", "validation")
    return recs


def load(record):
    """12-lead signal (LEADS12 order, µV, native 257 Hz) and beat annotations (sample, symbol)."""
    rec = wfdb.rdrecord(str(INCART / record))
    names = [s.upper() for s in rec.sig_name]
    x = np.stack([rec.p_signal[:, names.index(lead.upper())] for lead in LEADS12]) * 1000.0
    ann = wfdb.rdann(str(INCART / record), "atr")
    keep = np.isin(ann.symbol, list(BEATS))
    return x, np.asarray(ann.sample)[keep], np.asarray(ann.symbol)[keep]


def match(ref_t, det_t, tol):
    """One-to-one nearest matching within tol (same units); returns pairs (i_ref, j_det)."""
    pairs, used = [], set()
    if not len(det_t):
        return pairs
    for i in np.argsort(ref_t):
        d = np.abs(det_t - ref_t[i])
        for j in np.argsort(d):
            if d[j] > tol:
                break
            if j not in used:
                used.add(j)
                pairs.append((i, j))
                break
    return pairs


def evaluate_record(args):
    """Per-record counts and the scored beats for one detector and one rule set."""
    record, method, rules, p_corr_min = args
    x, ann_s, ann_y = load(record)
    n_native = x.shape[1]
    win, stride, guard = WIN_S * INCART_FS, STRIDE_S * INCART_FS, GUARD_S * INCART_FS
    rows = []
    for start in range(0, n_native - win + 1, stride):
        seg = lead8(preprocess(x[:, start:start + win], INCART_FS, "incart"))
        beats = classify(seg, method, rules, p_corr_min)
        lo_t, hi_t = guard / INCART_FS, (guard + stride) / INCART_FS        # central 10 s, seconds
        det_t = beats.fid.to_numpy() / FS if len(beats) else np.array([])
        in_det = (det_t >= lo_t) & (det_t < hi_t)
        r = (ann_s >= start) & (ann_s < start + win)
        ref_t = (ann_s[r] - start) / INCART_FS
        in_ref = (ref_t >= lo_t) & (ref_t < hi_t)
        ref_t, ref_y = ref_t[in_ref], ann_y[r][in_ref]
        dets = beats[in_det].reset_index(drop=True)
        eligible = bool(len(beats)) and bool(beats.eligible.iloc[0])
        pairs = match(ref_t, det_t[in_det], EC57_MATCH_MS / 1000.0)
        ref_of_det = {j: ref_y[i] for i, j in pairs}
        det_of_ref = {i: j for i, j in pairs}
        for i, y in enumerate(ref_y):
            j = det_of_ref.get(i)
            feat = dets.iloc[j] if j is not None else None
            rows.append({"record": record, "eligible": eligible, "ref": y, "det": j is not None,
                         "cls": dets.cls[j] if j is not None else "missed",
                         "qrs_ms": feat.qrs_ms if feat is not None else np.nan,
                         "rr_ratio": feat.rr_ms / feat.med_rr if feat is not None else np.nan,
                         "corr": feat["corr"] if feat is not None else np.nan,
                         "p_score": feat.p_score if feat is not None else np.nan,
                         "fails": ",".join(failed_criteria(dets.iloc[j], rules, p_corr_min))
                         if j is not None and dets.cls[j] not in ("PVC", "unclassified") else ""})
        for j in range(len(dets)):
            if j not in ref_of_det:
                rows.append({"record": record, "eligible": eligible, "ref": "none", "det": True,
                             "cls": dets.cls[j], "fails": "",
                             "qrs_ms": dets.qrs_ms[j], "rr_ratio": dets.rr_ms[j] / dets.med_rr[j],
                             "corr": dets["corr"][j], "p_score": dets.p_score[j]})
    return pd.DataFrame(rows)


def counts(df):
    """Per-record detection counts (all windows) and V-class counts (eligible windows only, i.e.
    windows with a valid SR template; docs/decisions.md item 6)."""
    g = []
    for record, d_all in df.groupby("record"):
        ref_all = d_all[d_all.ref != "none"]
        d = d_all[d_all.eligible]
        ref = d[d.ref != "none"]
        is_v = ref.ref.isin(V_REF)
        pvc = d.cls == "PVC"
        g.append({"record": record,
                  "n_ref": len(ref_all), "n_det": int(d_all.det.sum()), "matched": int(ref_all.det.sum()),
                  "v_all": int(ref_all.ref.isin(V_REF).sum()),
                  "tp": int((is_v & (ref.cls == "PVC")).sum()),
                  "fn": int((is_v & (ref.cls != "PVC")).sum()),
                  "fp": int((pvc & ~d.ref.isin(V_REF | NEUTRAL)).sum())})
    return pd.DataFrame(g)


def gross(c):
    s = c[["n_ref", "n_det", "matched", "v_all", "tp", "fn", "fp"]].sum()
    return {"qrs_se": s.matched / s.n_ref, "qrs_ppv": s.matched / s.n_det,
            "v_se": s.tp / (s.tp + s.fn) if s.tp + s.fn else np.nan,
            "v_ppv": s.tp / (s.tp + s.fp) if s.tp + s.fp else np.nan,
            "n_v": int(s.tp + s.fn), "v_ineligible_share": 1 - (s.tp + s.fn) / s.v_all if s.v_all else np.nan,
            "n_ref": int(s.n_ref)}


def bootstrap(c):
    """Percentile 95% CIs of the gross statistics over record-level resamples."""
    rng = np.random.default_rng(SEED)
    stats = [gross(c.iloc[rng.integers(0, len(c), len(c))]) for _ in range(N_BOOT)]
    return {k: (float(np.nanpercentile([s[k] for s in stats], 2.5)),
                float(np.nanpercentile([s[k] for s in stats], 97.5)))
            for k in ("qrs_se", "qrs_ppv", "v_se", "v_ppv")}


def run(subset, method, rules=RULES, p_corr_min=P_CORR_MIN):
    recs = split()
    names = recs.record[recs.subset == subset].tolist()
    with ProcessPoolExecutor(max_workers=max(os.cpu_count() - 2, 1)) as ex:
        df = pd.concat(ex.map(evaluate_record, [(r, method, rules, p_corr_min) for r in names]))
    return df, counts(df)


def fmt(v, ci):
    return f"{v:.3f} (95% CI {ci[0]:.3f} to {ci[1]:.3f})"


def report_tuning(results, recs):
    tun = recs[recs.subset == "tuning"]
    lines = ["# WP1: detector on the INCART tuning subset", "",
             "Generated by `python -m src.incart tuning`. Tuning subset only; the validation subset",
             "is not touched (gate K1). 14-s windows, 10-s stride, central 10 s scored; EC57 matching",
             "within 150 ms; 95% CIs from 2000 record-level bootstrap resamples.", "",
             f"Tuning subset: {tun.patient.nunique()} patients, {len(tun)} records.", "",
             "Rules: " + ", ".join(f"`{k}` = {v}" for k, v in RULES.items()) + f", `p_corr_min` = {P_CORR_MIN}.", "",
             "## Detection and V class by detector", "",
             "QRS detection over all windows; V class over eligible windows (a valid SR template).", "",
             "| Detector | QRS Se | QRS PPV | V Se | V PPV | V beats (eligible) | V in ineligible windows |",
             "|---|---|---|---|---|---|---|"]
    for method, (df, c) in results.items():
        g, ci = gross(c), bootstrap(c)
        lines.append(f"| {method} | {fmt(g['qrs_se'], ci['qrs_se'])} | {fmt(g['qrs_ppv'], ci['qrs_ppv'])} | "
                     f"{fmt(g['v_se'], ci['v_se'])} | {fmt(g['v_ppv'], ci['v_ppv'])} | {g['n_v']} | "
                     f"{g['v_ineligible_share']:.1%} |")
    lines += ["", "K1 (validation subset, point estimates): V Se ≥ 0.95 and V PPV ≥ 0.95; fallback PPV ≥ 0.95",
              "and Se ≥ 0.85 (docs/decisions.md item 8).", ""]
    for method, (df, c) in results.items():
        ref = df[(df.ref != "none") & df.eligible]
        v = ref[ref.ref.isin(V_REF)]
        fn = v[v.cls != "PVC"]
        lines += [f"## {method}: missed V beats ({len(fn)} of {len(v)})", "",
                  "| Outcome | n |", "|---|---|"]
        outcome = fn.cls.where(fn.cls != "other", "other: " + fn.fails.replace("", "none"))
        for k, n in outcome.value_counts().items():
            lines.append(f"| {k} | {n} |")
        fp = df[df.eligible & (df.cls == "PVC") & ~df.ref.isin(V_REF | NEUTRAL)]
        lines += ["", f"False-positive PVCs by reference class ({len(fp)}):", "",
                  "| Reference | n |", "|---|---|"]
        for k, n in fp.ref.value_counts().items():
            lines.append(f"| {k} | {n} |")
        lines.append("")
    (DOCS / "wp1_incart_tuning.md").write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("subset", choices=["tuning", "validation"])
    ap.add_argument("--method", choices=METHODS)
    ap.add_argument("--confirm-locked", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    recs = split()
    if a.subset == "tuning":
        results = {m: run("tuning", m) for m in ([a.method] if a.method else METHODS)}
        for m, (df, c) in results.items():
            df.to_parquet(OUT / f"tuning_beats_{m}.parquet")
            c.to_csv(OUT / f"tuning_counts_{m}.csv", index=False)
        report_tuning(results, recs)
        return
    a.method = a.method or DETECTOR
    if not a.confirm_locked:
        raise SystemExit("validation needs --confirm-locked (rules locked at K1, docs/decisions.md)")
    log = json.loads(VALIDATION_LOG.read_text()) if VALIDATION_LOG.exists() else []
    if len(log) >= MAX_VALIDATION_RUNS:
        raise SystemExit(f"validation subset already used {len(log)} times (plan item 18)")
    df, c = run("validation", a.method)
    g, ci = gross(c), bootstrap(c)
    log.append({"method": a.method, "rules": RULES, "p_corr_min": P_CORR_MIN, "result": g, "ci": ci})
    VALIDATION_LOG.write_text(json.dumps(log, indent=2))
    df.to_parquet(OUT / f"validation_beats_run{len(log)}.parquet")
    print(json.dumps({"result": g, "ci": ci}, indent=2))


if __name__ == "__main__":
    main()
