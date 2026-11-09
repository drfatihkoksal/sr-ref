"""Paths, signal constants, seeds, beat-rule starting values and gate thresholds.

Local data paths come from config/data_sources.local.yaml (not committed; see
config/data_sources.example.yaml). Row-level outputs go to the `artifacts` path set there,
outside the repository. Values follow prereg/ANALYSIS_PLAN_draft.md; beat rules are starting
values until they are locked at gate K1.
"""
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
DOCS = ROOT / "docs"
PLAN = ROOT / "prereg" / "ANALYSIS_PLAN.md"          # exists only after the freeze (gate K4)


def _sources():
    path = CONFIG / "data_sources.local.yaml"
    if not path.exists():
        raise FileNotFoundError(f"{path} missing: copy data_sources.example.yaml and set the paths")
    return yaml.safe_load(path.read_text())


SOURCES = _sources()
IVA = Path(SOURCES["ningbo_iva"])
IVA_RAW = IVA / "PVCVTRawECGData"
IVA_DIAGNOSIS = IVA / "download" / "Diagnosis.xlsx"
INCART = Path(SOURCES["incart"])
PTBXL = Path(SOURCES["ptbxl"])
LUDB = Path(SOURCES["ludb"])
MIMIC = Path(SOURCES["mimic_iv_ecg"])
MIMIC_IV = Path(SOURCES["mimic_iv"])
HEEDB = Path(SOURCES["heedb"])
RECON_ORIGIN = Path(SOURCES["recon_origin"])
RECON_ORIGIN_ARTIFACTS = Path(SOURCES["recon_origin_artifacts"])
BEAT_CLASSIFIER = Path(SOURCES["beat_classifier"])
ARTIFACTS = Path(SOURCES["artifacts"])

# ---- seeds and resampling ---------------------------------------------------------------
SEED = 1337
N_BOOT = 2000

# ---- signal (plan 4.1) ------------------------------------------------------------------
FS = 500                                  # Hz, every dataset after resampling
BAND = (0.5, 100.0)                       # Hz, Butterworth, zero-phase (forward-backward)
BAND_ORDER = 2                            # per pass; filtfilt doubles the effective order
NOTCH_Q = 30.0
MAINS = {"ningbo": (50.0, 100.0), "incart": (50.0,), "ptbxl": (50.0,),
         "mimic": (60.0,), "heedb": (60.0,)}
LEADS12 = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]
LEADS8 = ["I", "II", "V1", "V2", "V3", "V4", "V5", "V6"]   # independent leads (analysis space)

# Kors regression matrix, rows X, Y, Z; columns I, II, V1-V6 (Kors et al., Eur Heart J
# 1990;11:1083-1092). Axes: +X left, +Y inferior, +Z posterior. Same coefficients as
# HIDDEN-DELTA (stage2/hdsig.py), checked there against Vondrak & Penhakert, IEEE JTEHM
# 2022;10:1900208, Table 1.
KORS = np.array([[0.38, -0.07, -0.13, 0.05, -0.01, 0.14, 0.06, 0.54],
                 [-0.07, 0.93, 0.06, -0.02, -0.05, 0.06, -0.17, 0.13],
                 [0.11, -0.23, -0.43, -0.06, -0.14, -0.20, -0.11, 0.31]])

# QRS duration scale (docs/decisions.md item 5): the slope-envelope method measures the QRS
# shorter than the MIMIC-IV-ECG machine global QRS; this offset is added to every measured
# duration. `python -m src.qrs_scale` recomputes it and checks this value.
QRS_OFFSET_MS = 20.0

# ---- INCART split (plan 3.3) -----------------------------------------------------------
INCART_FS = 257
INCART_N_TUNING = 10                      # patients; the other 22 form the validation subset

# ---- beat rules (plan 4.3): candidate set from the INCART tuning subset (docs/wp1_incart_tuning.md,
# docs/decisions.md item 9); locked at gate K1. Concept starting values in brackets.
DETECTOR = "nk_vcg"
RULES = {
    "pvc_rr_ratio_max": 0.95,             # preceding RR < 0.95 x sinus RR [0.85]
    "pvc_qrs_ms_min": 115.0,              # [120]
    "pvc_qrs_ms_max": 200.0,              # plausibility bound, new: removes artefact and T-wave detections
    "pvc_corr_max": 0.85,                 # multi-lead correlation with the SR template [0.80]
    "propagate": True,                    # other wide, dissimilar beats of a PVC cluster are PVCs (new)
    "beatnet": True,                      # PVC also when the JECG BeatNet probability >= its card threshold
    "sr_qrs_ms_max": 120.0,
    "sr_rr_tol": 0.15,                    # RR within +-15% of the sinus RR
    "min_pvc": 1,
    "min_sr": 2,
}
EC57_MATCH_MS = 150.0                     # beat matching window (ANSI/AAMI EC57)

# ---- gates (plan 8) --------------------------------------------------------------------
K1_SE = 0.95
K1_PPV = 0.95
K1_KAPPA = 0.80
K3_VALID_SHARE = 0.95
K5_AUC_M1 = 0.75


def frozen():
    """SR-REF label-dependent analyses on Ningbo run only after the plan is frozen (gate K4)."""
    return PLAN.exists()
