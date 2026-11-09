"""Beat detection, delineation, clustering and rule-based classification (analysis plan 4.2-4.3).

Input: one preprocessed segment, 8 leads (I, II, V1-V6) x samples, µV, 500 Hz. Output
(`classify`): one row per detected beat with the fiducial, QRS onset and duration, preceding RR,
correlation with the provisional SR template, P-wave score, cluster and class (PVC, SR, other,
unclassified).

Detection candidates, to be chosen on the INCART tuning subset: NeuroKit2 on lead II
(`nk_ii`), WFDB XQRS on lead II (`xqrs_ii`), NeuroKit2 on the VCG magnitude (`nk_vcg`). Every
detector mark is moved to the peak of the smoothed spatial slope envelope within ±60 ms, so the
fiducial does not depend on the detector. QRS bounds use the RECON-ORIGIN slope-envelope method
(`qrs_window`, ported unchanged); durations add `QRS_OFFSET_MS` (machine scale, `src.qrs_scale`). Clusters: average-linkage correlation clustering cut at
r = 0.85 on ±120 ms windows (RECON-ORIGIN decision 8). The provisional SR template is the median
beat of the largest cluster whose median beat is narrow; a segment without one is ineligible. A P wave is "present" when the segment
from 300 to 40 ms before the beat's QRS onset matches the same segment of the SR template
(multi-lead correlation, lag up to ±40 ms), i.e. a P wave at the usual PR interval. The P wave is
used by the SR rule only; it was removed from the PVC rule (docs/decisions.md item 7).
"""
import warnings

import neurokit2 as nk
import numpy as np
import pandas as pd
import wfdb.processing
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.signal import detrend

from src.beatnet import probs as beatnet_probs, threshold as beatnet_threshold
from src.config import FS, QRS_OFFSET_MS, RULES
from src.signal import vcg

MS = FS / 1000.0
HALF = int(200 * MS)                      # beat window ±200 ms (qrs_window geometry)
CLUSTER_HALF = int(120 * MS)
CLUSTER_CUT = 0.15                        # correlation distance, r >= 0.85
SNAP = int(60 * MS)                       # detector mark -> envelope peak within ±60 ms
MERGE = int(200 * MS)                     # marks closer than 200 ms are one beat
GAP = int(60 * MS)                        # envelope dips up to 60 ms are bridged
PRE = int(500 * MS)                       # pre-QRS context for the P-wave score
P_SEG = (int(300 * MS), int(40 * MS))     # P segment: 300 to 40 ms before QRS onset
P_LAG = int(40 * MS)
TEMPLATE_LAG = int(20 * MS)
P_CORR_MIN = 0.50                         # starting value; tuned on INCART (K1)
METHODS = ("nk_ii", "xqrs_ii", "nk_vcg")


def envelope(x):
    """Spatial slope envelope of x (leads x samples), 10 ms moving average."""
    return np.convolve(np.sqrt((np.gradient(x, axis=-1) ** 2).sum(0)), np.ones(5) / 5, "same")


def qrs_window(X):
    """QRS onset and offset (samples) in a ±200 ms beat window X (leads x samples), RECON-ORIGIN
    plan 5.2: the envelope peak is searched within ±80 ms of the centre; the QRS is the run above
    floor + 15% of (peak - floor), floor = the envelope median, that contains the peak, after
    bridging dips of up to 60 ms."""
    env = envelope(X)
    c = X.shape[1] // 2
    p = c - 40 + int(np.argmax(env[c - 40:c + 40]))
    floor = np.median(env)
    above = np.flatnonzero(env > floor + 0.15 * (env[p] - floor))
    runs = np.split(above, np.flatnonzero(np.diff(above) > GAP + 1) + 1)
    run = next(r for r in runs if r[0] <= p <= r[-1])
    return max(run[0] - 1, 0), min(run[-1] + 1, X.shape[1] - 1)


def _marks(x8, method):
    lead_ii = x8[1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if method == "nk_ii":
            return np.asarray(nk.ecg_peaks(lead_ii, sampling_rate=FS, method="neurokit")[1]["ECG_R_Peaks"])
        if method == "xqrs_ii":
            return np.asarray(wfdb.processing.xqrs_detect(lead_ii / 1000.0, fs=FS, verbose=False))
        if method == "nk_vcg":
            mag = np.linalg.norm(vcg(x8), axis=0)
            return np.asarray(nk.ecg_peaks(mag, sampling_rate=FS, method="neurokit")[1]["ECG_R_Peaks"])
    raise ValueError(method)


def detect(x8, method):
    """Fiducials (samples): detector marks snapped to the envelope peak, duplicates merged."""
    try:
        marks = _marks(x8, method)
    except Exception:                     # a detector failing on a flat or noisy segment: no beats
        return np.array([], dtype=int)
    env = envelope(x8)
    n = x8.shape[1]
    fid = []
    for m in np.asarray(marks, dtype=int):
        lo, hi = max(m - SNAP, 0), min(m + SNAP + 1, n)
        if hi > lo:
            fid.append(lo + int(np.argmax(env[lo:hi])))
    fid = np.unique(fid)
    keep = []
    for f in fid:
        if keep and f - keep[-1] < MERGE:
            if env[f] > env[keep[-1]]:
                keep[-1] = f
        else:
            keep.append(f)
    return np.asarray(keep, dtype=int)


def _windows(x8, fid, before, after):
    return np.stack([x8[:, f - before:f + after + 1] for f in fid]) if len(fid) else np.empty((0, 8, before + after + 1))


def _mcorr(a, b, trend=False):
    """Multi-lead correlation: per-lead centring (or linear detrending), then one correlation
    over all leads."""
    if trend:
        a, b = detrend(a, axis=-1), detrend(b, axis=-1)
    a = a - a.mean(-1, keepdims=True)
    b = b - b.mean(-1, keepdims=True)
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else 0.0


def _best_corr(seg_fn, ref, lag, trend=False):
    """Max multi-lead correlation of ref with seg_fn(shift) over shifts in [-lag, lag]."""
    best = -1.0
    for s in range(-lag, lag + 1):
        seg = seg_fn(s)
        if seg is not None:
            best = max(best, _mcorr(seg, ref, trend))
    return best


def cluster(W):
    """Cluster labels of beat windows W (beats x leads x samples), RECON-ORIGIN decision 8."""
    if len(W) < 2:
        return np.ones(len(W), dtype=int)
    X = detrend(W, axis=-1)
    X = (X - X.mean(axis=(1, 2), keepdims=True)) / (X.std(axis=(1, 2), keepdims=True) + 1e-12)
    return fcluster(linkage(X.reshape(len(X), -1), "average", metric="correlation"), CLUSTER_CUT, "distance")


def classify(x8, method, rules=RULES, p_corr_min=P_CORR_MIN):
    """Beat table of one segment; see the module docstring."""
    n = x8.shape[1]
    fid = detect(x8, method)
    fid = fid[(fid - PRE - P_LAG >= 0) & (fid + HALF < n)]
    cols = ["fid", "on", "qrs_ms", "rr_ms", "med_rr", "corr", "p_score", "cluster", "prob", "eligible", "cls"]
    if not len(fid):
        return pd.DataFrame(columns=cols)
    W = _windows(x8, fid, HALF, HALF)
    bounds = np.array([qrs_window(w) for w in W])
    on_abs = fid - HALF + bounds[:, 0]
    qrs_ms = (bounds[:, 1] - bounds[:, 0]) / MS + QRS_OFFSET_MS
    rr = np.r_[np.nan, np.diff(fid) / MS]
    lab = cluster(W[:, :, HALF - CLUSTER_HALF:HALF + CLUSTER_HALF + 1])

    # provisional SR template: largest cluster (>= 2 beats) whose median beat is narrow
    dom, tmpl = None, None
    for c in sorted(set(lab), key=lambda c: -(lab == c).sum()):
        if (lab == c).sum() < 2:
            break
        T = np.median(W[lab == c], axis=0)
        on_t, off_t = qrs_window(T)
        if (off_t - on_t) / MS + QRS_OFFSET_MS < rules["sr_qrs_ms_max"]:
            dom, tmpl, tmpl_on = c, T, on_t
            break
    out = pd.DataFrame({"fid": fid, "on": on_abs, "qrs_ms": qrs_ms, "rr_ms": rr, "cluster": lab,
                        "prob": beatnet_probs(x8, fid) if rules.get("beatnet", False) else np.nan})
    if dom is None:
        out["med_rr"], out["corr"], out["p_score"], out["cls"] = np.nan, np.nan, np.nan, "unclassified"
        out["eligible"] = False
        return out[cols]

    # correlation with the SR template (±120 ms, lag up to ±20 ms)
    ref = tmpl[:, HALF - CLUSTER_HALF:HALF + CLUSTER_HALF + 1]

    def shifted(f, s):
        a, b = f + s - CLUSTER_HALF, f + s + CLUSTER_HALF + 1
        return x8[:, a:b] if a >= 0 and b <= n else None

    corr = np.array([_best_corr(lambda s, f=f: shifted(f, s), ref, TEMPLATE_LAG) for f in fid])

    # P-wave score: segment 300-40 ms before QRS onset vs the same segment of the SR template
    sr_fid = fid[lab == dom]
    long_t = np.median(_windows(x8, sr_fid, PRE, HALF), axis=0)       # template with pre-QRS context
    t_on = PRE - HALF + tmpl_on
    p_ref = long_t[:, t_on - P_SEG[0]:t_on - P_SEG[1]]

    def p_seg(o, s):
        a, b = o + s - P_SEG[0], o + s - P_SEG[1]
        return x8[:, a:b] if a >= 0 and b <= n else None

    p_score = np.array([_best_corr(lambda s, o=o: p_seg(o, s), p_ref, P_LAG, trend=True) for o in on_abs])
    p_present = p_score >= p_corr_min

    sr_pairs = (lab[1:] == dom) & (lab[:-1] == dom)
    # sinus cycle: median SR-SR interval; without two such pairs (e.g. bigeminy) the mean RR,
    # which equals the sinus cycle when pauses are compensatory
    med_rr = np.median(rr[1:][sr_pairs]) if sr_pairs.sum() >= 2 else np.nanmean(rr)
    cls = []
    for k in range(len(fid)):
        if np.isnan(rr[k]):
            cls.append("unclassified")
            continue
        is_pvc = (rr[k] < rules["pvc_rr_ratio_max"] * med_rr
                  and rules["pvc_qrs_ms_min"] <= qrs_ms[k] <= rules["pvc_qrs_ms_max"]
                  and corr[k] < rules["pvc_corr_max"])
        is_sr = (lab[k] == dom and p_present[k] and qrs_ms[k] < rules["sr_qrs_ms_max"]
                 and abs(rr[k] / med_rr - 1) <= rules["sr_rr_tol"] and (k == 0 or cls[-1] != "PVC"))
        cls.append("PVC" if is_pvc else "SR" if is_sr else "other")
    if rules.get("propagate", False):
        # a cluster with a premature, wide, dissimilar beat is a PVC morphology: its other wide,
        # dissimilar beats are PVCs too, premature or not (late PVCs, ventricular runs)
        cls = np.array(cls, dtype=object)
        for c in set(lab[cls == "PVC"]) - {dom}:
            add = ((lab == c) & (cls == "other") & (qrs_ms >= rules["pvc_qrs_ms_min"])
                   & (qrs_ms <= rules["pvc_qrs_ms_max"]) & (corr < rules["pvc_corr_max"]))
            cls[add] = "PVC"
        cls = list(cls)
    if rules.get("beatnet", False):
        # second cue (docs/decisions.md item 14): PVC when BeatNet >= its card threshold; a beat
        # right after a PVC is then no longer SR
        cls = np.array(cls, dtype=object)
        cls[(out.prob.to_numpy() >= beatnet_threshold()) & (cls != "unclassified")] = "PVC"
        after = np.r_[False, cls[:-1] == "PVC"]
        cls[after & (cls == "SR")] = "other"
        cls = list(cls)
    out["med_rr"], out["corr"], out["p_score"], out["cls"] = med_rr, corr, p_score, cls
    out["eligible"] = True
    return out[cols]


def failed_criteria(row, rules=RULES, p_corr_min=P_CORR_MIN):
    """Which PVC criteria a classified beat fails (for tuning reports)."""
    fails = []
    if row.rr_ms >= rules["pvc_rr_ratio_max"] * row.med_rr:
        fails.append("rr")
    if row.qrs_ms < rules["pvc_qrs_ms_min"]:
        fails.append("qrs")
    if row.qrs_ms > rules["pvc_qrs_ms_max"]:
        fails.append("qrs_max")
    if row["corr"] >= rules["pvc_corr_max"]:
        fails.append("corr")
    return fails
