"""QRS duration scale (docs/decisions.md item 5).

The slope-envelope QRS method (`beats.qrs_window`, ported from RECON-ORIGIN, threshold 15%)
measures the QRS shorter than clinical measurements. A constant offset puts it on the machine
scale, on which the 120 ms thresholds and the 12SL exclusion statements are defined.

- **Primary reference: MIMIC-IV-ECG machine global QRS** (`qrs_end − qrs_onset`,
  `machine_measurements.csv`). Sample: 500 ECGs drawn with `SEED` among those whose first three
  machine statements say sinus rhythm and mention no block, pacing, AF/flutter, ventricular
  rhythm or pre-excitation. Per ECG, the median beat of the largest cluster is measured. The
  offset is the negative median bias, rounded to 2 ms; `config.QRS_OFFSET_MS` must equal it.
  Label-free; no H4 quantity is involved.
- **Secondary check: LUDB expert global QRS** (earliest onset to latest offset over at least 6
  annotated leads, 200 records). Reported, not used for the offset: the LUDB convention runs
  longer than the machine convention.

    python -m src.qrs_scale   -> docs/wp1_qrs_scale.md
"""
import warnings

import numpy as np
import pandas as pd
import wfdb

from src.beats import CLUSTER_HALF, HALF, cluster, detect, qrs_window
from src.config import DOCS, LEADS12, LUDB, MIMIC, QRS_OFFSET_MS, SEED
from src.signal import lead8, preprocess

N_MIMIC = 500
EXCLUDE = "block|pac|fibril|flutter|ventricular|wpw|pre-excit"
LUDB_LEADS = ["i", "ii", "iii", "avr", "avl", "avf", "v1", "v2", "v3", "v4", "v5", "v6"]
MIN_LEADS = 6
BIN = 25                                  # samples (50 ms): LUDB beats grouped across leads by peak


def _fiducials(x8):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fid = detect(x8, "nk_vcg")
    return fid[(fid - HALF >= 0) & (fid + HALF < x8.shape[1])]


def _template_ms(W):
    lab = cluster(W[:, :, HALF - CLUSTER_HALF:HALF + CLUSTER_HALF + 1])
    T = np.median(W[lab == np.bincount(lab).argmax()], axis=0)
    return np.diff(qrs_window(T))[0] * 2.0


def mimic_sample():
    mm = pd.read_csv(MIMIC / "machine_measurements.csv",
                     usecols=["study_id", "report_0", "report_1", "report_2", "qrs_onset", "qrs_end"])
    mm = mm.merge(pd.read_csv(MIMIC / "record_list.csv")[["study_id", "path"]], on="study_id")
    mm["qrs"] = mm.qrs_end - mm.qrs_onset
    rep = (mm.report_0.fillna("") + " " + mm.report_1.fillna("") + " " + mm.report_2.fillna("")).str.lower()
    ok = (mm.qrs > 40) & (mm.qrs < 250) & rep.str.contains("sinus rhythm") & ~rep.str.contains(EXCLUDE)
    return mm[ok].sample(N_MIMIC, random_state=SEED)


def mimic():
    """(method ms without offset, machine ms) per ECG."""
    out = []
    for _, r in mimic_sample().iterrows():
        rec = wfdb.rdrecord(str(MIMIC / r.path))
        names = [n.upper() for n in rec.sig_name]
        x = np.stack([rec.p_signal[:, names.index(lead.upper())] for lead in LEADS12]) * 1000.0
        if np.isnan(x).any():
            continue
        x8 = lead8(preprocess(x, rec.fs, "mimic"))
        fid = _fiducials(x8)
        if len(fid) >= 4:
            out.append((_template_ms(np.stack([x8[:, f - HALF:f + HALF + 1] for f in fid])), r.qrs))
    return np.array(out)


def ludb_expert(rec):
    on, off = {}, {}
    for lead in LUDB_LEADS:
        a = wfdb.rdann(str(LUDB / str(rec)), lead)
        s, y = np.asarray(a.sample), np.asarray(a.symbol)
        for k in np.flatnonzero(y == "N"):
            if 0 < k < len(y) - 1 and y[k - 1] == "(" and y[k + 1] == ")":
                key = int(round(s[k] / BIN))
                on.setdefault(key, []).append(s[k - 1])
                off.setdefault(key, []).append(s[k + 1])
    return np.array(sorted((min(on[k]), max(off[k])) for k in on if len(on[k]) >= MIN_LEADS))


def ludb():
    """(method ms without offset, expert ms) per record; the template uses the beats inside an
    expert QRS."""
    out = []
    for rec in range(1, 201):
        r = wfdb.rdrecord(str(LUDB / str(rec)))
        names = [n.lower() for n in r.sig_name]
        x8 = lead8(preprocess(np.stack([r.p_signal[:, names.index(lead)] for lead in LUDB_LEADS]) * 1000.0,
                              r.fs, "ptbxl"))
        ex = ludb_expert(rec)
        wins, true_ms = [], []
        for f in _fiducials(x8):
            m = np.flatnonzero((ex[:, 0] <= f) & (ex[:, 1] >= f)) if len(ex) else []
            if len(m):
                wins.append(x8[:, f - HALF:f + HALF + 1])
                true_ms.append((ex[m[0], 1] - ex[m[0], 0]) * 2.0)
        if len(wins) >= 3:
            out.append((np.diff(qrs_window(np.median(wins, axis=0)))[0] * 2.0, np.median(true_ms)))
    return np.array(out)


def q(v):
    return f"median {np.median(v):+.0f} ms (IQR {np.percentile(v, 25):+.0f} to {np.percentile(v, 75):+.0f})"


def main():
    m, lu = mimic(), ludb()
    bias = m[:, 0] - m[:, 1]
    offset = -2.0 * round(np.median(bias) / 2.0)
    resid = bias + offset
    lines = ["# WP1: QRS duration scale", "",
             "Generated by `python -m src.qrs_scale`. Method: slope-envelope QRS on the median beat",
             "(`beats.qrs_window`, RECON-ORIGIN plan 5.2, threshold 15%). See the module docstring.", "",
             "## Primary: MIMIC-IV-ECG machine global QRS", "",
             f"- ECGs measured: {len(m)} of {N_MIMIC} sampled; machine QRS median {np.median(m[:, 1]):.0f} ms "
             f"(IQR {np.percentile(m[:, 1], 25):.0f} to {np.percentile(m[:, 1], 75):.0f}).",
             f"- Bias (method − machine): {q(bias)}.",
             f"- **Offset: +{offset:.0f} ms.** `config.QRS_OFFSET_MS` = {QRS_OFFSET_MS:.0f}"
             + (" (matches)." if offset == QRS_OFFSET_MS else " (**does not match: update config**)."),
             f"- Residual after correction: {q(resid)}; mean absolute error {np.mean(np.abs(resid)):.0f} ms; "
             f"within ±10 ms {np.mean(np.abs(resid) <= 10):.0%}.", "",
             "## Secondary: LUDB expert global QRS", "",
             f"- Records: {len(lu)}; expert QRS median {np.median(lu[:, 1]):.0f} ms "
             f"(IQR {np.percentile(lu[:, 1], 25):.0f} to {np.percentile(lu[:, 1], 75):.0f}).",
             f"- Bias (method − expert), no offset: {q(lu[:, 0] - lu[:, 1])}; with the offset: "
             f"{q(lu[:, 0] + offset - lu[:, 1])}.",
             "- The LUDB convention (union of per-lead bounds) runs longer than the machine convention."]
    (DOCS / "wp1_qrs_scale.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
