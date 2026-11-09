"""Pre-processing shared by every dataset (analysis plan 4.1).

Resampling to 500 Hz (polyphase), zero-phase Butterworth band-pass 0.5-100 Hz, zero-phase mains
notch per dataset, the 8 independent leads, and the Kors VCG. Signals are leads x samples in µV.
"""
from fractions import Fraction

import numpy as np
from scipy.signal import butter, filtfilt, iirnotch, resample_poly, sosfiltfilt

from src.config import BAND, BAND_ORDER, FS, KORS, LEADS8, LEADS12, MAINS, NOTCH_Q

SOS = butter(BAND_ORDER, BAND, btype="bandpass", fs=FS, output="sos")
NOTCHES = {f0: iirnotch(f0, NOTCH_Q, FS) for f0 in (50.0, 60.0, 100.0)}
IDX8 = [LEADS12.index(lead) for lead in LEADS8]


def resample(x, fs):
    """Polyphase resampling of x (leads x samples) from fs to FS."""
    if fs == FS:
        return np.asarray(x, dtype=float)
    f = Fraction(FS, int(round(fs))).limit_denominator(1000)
    return resample_poly(np.asarray(x, dtype=float), f.numerator, f.denominator, axis=-1)


def preprocess(x, fs, dataset):
    """Resample, band-pass and notch x (leads x samples). Pacing screens run before this, at the
    native rate (plan 4.1)."""
    y = sosfiltfilt(SOS, resample(x, fs), axis=-1)
    for f0 in MAINS[dataset]:
        y = filtfilt(*NOTCHES[f0], y, axis=-1)
    return y


def lead8(x12):
    """The 8 independent leads I, II, V1-V6 of a 12-lead array in LEADS12 order."""
    return x12[IDX8]


def vcg(x8):
    """Kors VCG (X, Y, Z) of an 8-lead array (I, II, V1-V6)."""
    return KORS @ x8
