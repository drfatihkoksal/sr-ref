"""JECG beat classifier (BeatNet), used unchanged as the second PVC cue (docs/decisions.md item 14).

MIT-BIH-trained CNN of the JECG paper (Köksal, J Electrocardiol 2026;154486), the same weights and
input as in RECON-ORIGIN (`src/signals.py`): leads II and V1, ±200 ms around the beat, each beat
z-scored per lead, plus the segment's median beat and the beat minus that median (6 channels), and
the preceding and following RR relative to the segment's median RR. Input band 0.5-40 Hz (a 40 Hz
low-pass on the SR-REF signal). Operating threshold: the model card's (0.99), not tuned on INCART.
"""
import json

import numpy as np
import torch
import torch.nn as nn
from scipy.signal import butter, sosfiltfilt

from src.config import BEAT_CLASSIFIER, FS, LEADS8

HALF_WIN = 100                                   # samples at 500 Hz, ±200 ms
LOWPASS = butter(3, 40.0, "lowpass", fs=FS, output="sos")
IDX = [LEADS8.index("II"), LEADS8.index("V1")]
_NET = None


class BeatNet(nn.Module):
    def __init__(s, n_ctx=3):
        super().__init__()
        s.f = nn.Sequential(
            nn.Conv1d(6, 32, 7, 1, 3), nn.BatchNorm1d(32), nn.ReLU(), nn.MaxPool1d(2),
            nn.Conv1d(32, 64, 5, 1, 2), nn.BatchNorm1d(64), nn.ReLU(), nn.MaxPool1d(2),
            nn.Conv1d(64, 64, 3, 1, 1), nn.BatchNorm1d(64), nn.ReLU(), nn.AdaptiveAvgPool1d(25))
        s.h = nn.Sequential(
            nn.Flatten(), nn.Dropout(0.3),
            nn.Linear(64 * 25 + n_ctx, 128), nn.ReLU(), nn.Dropout(0.3), nn.Linear(128, 1))

    def forward(s, x, ctx):
        return s.h(torch.cat([s.f(x).flatten(1), ctx], 1)).squeeze(-1)


def threshold():
    return float(json.loads((BEAT_CLASSIFIER.parent / "model_card.json").read_text())["threshold"])


def _net():
    global _NET
    if _NET is None:
        torch.set_num_threads(1)
        _NET = BeatNet()
        _NET.load_state_dict(torch.load(BEAT_CLASSIFIER, map_location="cpu"))
        _NET.eval()
    return _NET


@torch.no_grad()
def probs(x8, fid):
    """PVC probability per beat of a segment (8 leads x samples, 500 Hz); NaN near the edges."""
    sf = sosfiltfilt(LOWPASS, x8[IDX], axis=-1)
    rr = np.diff(fid) / FS
    med = np.median(rr) if len(rr) else 1.0
    pre, post = np.r_[np.nan, rr] / med, np.r_[rr, np.nan] / med
    p = np.full(len(fid), np.nan)
    idx = np.flatnonzero((fid >= HALF_WIN) & (fid <= sf.shape[1] - HALF_WIN))
    if len(idx):
        b = np.stack([sf[:, f - HALF_WIN:f + HALF_WIN] for f in fid[idx]])
        b = (b - b.mean(-1, keepdims=True)) / (b.std(-1, keepdims=True) + 1e-6)
        t = np.median(b, axis=0)
        x = np.concatenate([b, np.repeat(t[None], len(b), 0), b - t], 1)
        a, c = np.nan_to_num(pre[idx], nan=1.0), np.nan_to_num(post[idx], nan=1.0)
        ctx = np.stack([np.clip(a, 0, 3), np.clip(c, 0, 3), np.clip(a + c, 0, 6) / 2], 1)
        p[idx] = torch.sigmoid(_net()(torch.tensor(x, dtype=torch.float32),
                                      torch.tensor(ctx, dtype=torch.float32))).numpy()
    return p
