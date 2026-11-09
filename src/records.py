"""Record-level pipeline (analysis plan 4.3-4.6): loaders, PVC and SR templates, eligibility.

`process(x12, fs, dataset)` pre-processes a 12-lead record (µV), classifies its beats with the
locked rules (`config.DETECTOR`, `config.RULES`), and builds one median template per class:

- beats are cut from 300 ms before to 700 ms after their fiducial (beats too close to the record
  edges are not used), aligned to the provisional median by cross-correlation of the 8 leads over
  ±120 ms (lags up to ±20 ms), and the sample-wise median is taken (12 leads);
- QRS onset and offset come from `beats.qrs_window` on the 8-lead template (+ `QRS_OFFSET_MS` for
  the duration); the template is cut from 100 ms before to 600 ms after the QRS onset;
- baseline: per lead, the mean of the 20 ms before the QRS onset is subtracted.

Lead quality (plan 4.5, "noisy or missing leads"; docs/decisions.md item 13): a record is excluded
when any of the 12 leads is flat (peak-to-peak < 50 µV) or holds a constant run longer than 0.5 s
(lead off, or leads recorded in sequential 2.5-s groups with the rest padded).

Eligibility: a valid SR template in the record, ≥ `min_pvc` PVC and ≥ `min_sr` SR beats with a
usable window, one PVC cluster (monomorphic), SR template QRS < 120 ms. OT-like: PVC template R
fraction r = R / (R + S) > 0.5 in II, III and aVF (R, S: the largest positive local maximum and the
deepest negative local minimum within the QRS, as in RECON-ORIGIN `features.extremum`).
"""
import warnings

import numpy as np
import wfdb

from src.beats import CLUSTER_HALF, HALF, MS, P_CORR_MIN, TEMPLATE_LAG, classify, qrs_window
from src.config import DETECTOR, LEADS12, QRS_OFFSET_MS, RULES
from src.signal import IDX8, preprocess

CUT = (int(300 * MS), int(700 * MS))      # beat cut around the fiducial
TEMPLATE = (int(100 * MS), int(600 * MS))  # template window around the QRS onset
BASE = int(20 * MS)
INFERIOR = ["II", "III", "aVF"]
FLAT_UV = 50.0
RUN_S = 0.5


def load(path):
    """12-lead record in LEADS12 order (µV) and its sampling rate; None if a lead is missing."""
    rec = wfdb.rdrecord(str(path))
    names = [n.upper() for n in rec.sig_name]
    if not all(lead.upper() in names for lead in LEADS12):
        return None, None
    x = np.stack([rec.p_signal[:, names.index(lead.upper())] for lead in LEADS12]).astype(float)
    units = [u.lower() for u in rec.units]
    scale = 1.0 if units and units[0] in ("uv", "µv") else 1000.0
    return x * scale, rec.fs


def lead_quality(x, fs):
    """Per lead: flat (peak-to-peak < FLAT_UV) or a constant run longer than RUN_S, on the raw
    record (before filtering)."""
    flat = np.ptp(x, axis=1) < FLAT_UV
    run = np.zeros(len(x), dtype=bool)
    for i, y in enumerate(x):
        same = np.r_[False, np.diff(y) == 0, False]
        idx = np.flatnonzero(np.diff(same.astype(int)))
        run[i] = bool(len(idx)) and (idx[1::2] - idx[::2]).max() > RUN_S * fs
    return flat | run


def extremum(y, on, off, sign):
    """Largest positive local maximum (sign +1) or deepest negative local minimum (sign -1) of y
    within [on, off]; 0 if there is none (RECON-ORIGIN features.extremum)."""
    z = sign * y
    t = np.arange(max(on, 1), min(off, len(y) - 2) + 1)
    t = t[(z[t] > z[t - 1]) & (z[t] >= z[t + 1]) & (z[t] > 0)]
    return float(z[t].max()) if len(t) else 0.0


def template(x12, fid):
    """Median template (12 leads, -100 to +600 ms around the QRS onset), QRS duration (ms), and
    the number of beats used; None if no beat has a full window."""
    n = x12.shape[1]
    fid = [f for f in fid if f - CUT[0] - TEMPLATE_LAG >= 0 and f + CUT[1] + TEMPLATE_LAG < n]
    if not fid:
        return None
    W = np.stack([x12[:, f - CUT[0]:f + CUT[1] + 1] for f in fid])
    c, h = CUT[0], CUT[0] - CLUSTER_HALF
    ref = np.median(W, axis=0)[IDX8, h:c + CLUSTER_HALF + 1]
    aligned = []
    for f in fid:
        best, best_s = -np.inf, 0
        for s in range(-TEMPLATE_LAG, TEMPLATE_LAG + 1):
            seg = x12[IDX8, f + s - CLUSTER_HALF:f + s + CLUSTER_HALF + 1]
            a, b = seg - seg.mean(-1, keepdims=True), ref - ref.mean(-1, keepdims=True)
            r = (a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum() + 1e-12)
            if r > best:
                best, best_s = r, s
        aligned.append(x12[:, f + best_s - CUT[0]:f + best_s + CUT[1] + 1])
    T = np.median(np.stack(aligned), axis=0)
    on, off = qrs_window(T[IDX8, c - HALF:c + HALF + 1])
    on, off = on + c - HALF, off + c - HALF
    lo, hi = on - TEMPLATE[0], on + TEMPLATE[1] + 1
    if lo < 0 or hi > T.shape[1]:
        return None
    T = T[:, lo:hi]
    T = T - T[:, TEMPLATE[0] - BASE:TEMPLATE[0]].mean(1, keepdims=True)
    return {"T": T, "on": TEMPLATE[0], "off": TEMPLATE[0] + off - on,
            "qrs_ms": (off - on) / MS + QRS_OFFSET_MS, "n": len(fid)}


def r_fraction(T, on, off, lead):
    y = T[LEADS12.index(lead)]
    r, s = extremum(y, on, off, 1), extremum(y, on, off, -1)
    return r / (r + s) if r + s > 0 else np.nan


PACED_BEFORE = int(150 * MS)               # a stimulus up to 150 ms before a beat's fiducial: paced


def process(x12_raw, fs, dataset, rules=RULES, stim=None, pvc_cluster="single", sr_fallback=False):
    """Beats, templates, eligibility and OT-like flag of one record; see the module docstring.
    stim: pacing-stimulus positions at 500 Hz (Ningbo); paced beats and their neighbours (±1 beat)
    are taken out of the PVC and SR classes (plan 4.5)."""
    if lead_quality(x12_raw, fs).any():
        return {"eligible": False, "reason": "flat or incomplete lead", "ot_like": False}
    x12 = preprocess(x12_raw, fs, dataset)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        beats = classify(x12[IDX8], DETECTOR, rules)
    if stim is not None and len(beats) and len(stim):
        lag = beats.fid.to_numpy()[:, None] - np.asarray(stim)[None, :]
        paced = ((lag >= 0) & (lag <= PACED_BEFORE)).any(1)
        near = paced | np.r_[paced[1:], False] | np.r_[False, paced[:-1]]
        beats["paced"] = paced
        beats.loc[near & beats.cls.isin(["PVC", "SR"]), "cls"] = "other"
    out = {"x12": x12, "beats": beats, "eligible": False, "reason": None, "ot_like": False}
    if not len(beats) or not bool(beats.eligible.iloc[0]):
        out["reason"] = "no SR template"
        return out
    pvc, sr = beats[beats.cls == "PVC"], beats[beats.cls == "SR"]
    if pvc.cluster.nunique() > 1:
        if pvc_cluster == "single":
            out["reason"] = "more than one PVC cluster"
            return out
        pvc = pvc[pvc.cluster == pvc.cluster.value_counts().idxmax()]      # dominant PVC cluster
    if sr_fallback and len(sr) < rules["min_sr"] and "p_score" in beats:
        # fallback: narrow beats of the SR-template cluster preceded by a P wave, whatever their RR
        # and whether a PVC precedes them (bigeminy leaves no beat that meets the strict SR rule)
        dom = beats.loc[beats.cls == "SR", "cluster"].mode()
        dom = dom.iloc[0] if len(dom) else beats.loc[beats.cls != "PVC", "cluster"].mode().iloc[0]
        sr = beats[(beats.cluster == dom) & (beats.cls != "PVC") & (beats.qrs_ms < rules["sr_qrs_ms_max"])
                   & (beats.p_score >= P_CORR_MIN) & ~beats.get("paced", False)]
        out["sr_fallback"] = True
    tp, ts = template(x12, pvc.fid.tolist()), template(x12, sr.fid.tolist())
    out["pvc"], out["sr"] = tp, ts
    if tp is None or tp["n"] < rules["min_pvc"]:
        out["reason"] = "no usable PVC"
    elif ts is None or ts["n"] < rules["min_sr"]:
        out["reason"] = "fewer than 2 usable SR beats"
    elif ts["qrs_ms"] >= rules["sr_qrs_ms_max"]:
        out["reason"] = "SR template QRS >= 120 ms"
    else:
        out["eligible"] = True
        out["r_inferior"] = [r_fraction(tp["T"], tp["on"], tp["off"], lead) for lead in INFERIOR]
        out["ot_like"] = all(r > 0.5 for r in out["r_inferior"] if not np.isnan(r)) and not any(
            np.isnan(out["r_inferior"]))
    return out


def from_beats(x12_raw, fs, dataset, pvc_fid, sr_fid, excluded=()):
    """Templates, eligibility and OT-like flag when the PVC and SR beats are given (Ningbo: the
    RECON-ORIGIN adjudication, docs/decisions.md item 22). Fiducials at 500 Hz in the
    pre-processed signal; `excluded` beats (paced and their neighbours) are dropped. Same template
    code, lead-quality rule and eligibility as `process`, without the beat classification."""
    if lead_quality(x12_raw, fs).any():
        return {"eligible": False, "reason": "flat or incomplete lead", "ot_like": False}
    x12 = preprocess(x12_raw, fs, dataset)
    ex = set(excluded)
    pvc_fid = [f for f in pvc_fid if f not in ex]
    sr_fid = [f for f in sr_fid if f not in ex]
    tp, ts = template(x12, pvc_fid), template(x12, sr_fid)
    out = {"x12": x12, "pvc": tp, "sr": ts, "eligible": False, "reason": None, "ot_like": False}
    if tp is None or tp["n"] < RULES["min_pvc"]:
        out["reason"] = "no usable PVC"
    elif ts is None or ts["n"] < RULES["min_sr"]:
        out["reason"] = "fewer than 2 usable SR beats"
    elif ts["qrs_ms"] >= RULES["sr_qrs_ms_max"]:
        out["reason"] = "SR template QRS >= 120 ms"
    else:
        out["eligible"] = True
        out["r_inferior"] = [r_fraction(tp["T"], tp["on"], tp["off"], lead) for lead in INFERIOR]
        out["ot_like"] = all(r > 0.5 for r in out["r_inferior"] if not np.isnan(r)) and not any(
            np.isnan(out["r_inferior"]))
    return out
