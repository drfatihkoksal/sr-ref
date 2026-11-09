"""H4 on a HEEDB site and the transform estimator's κ (analysis plan 4.7, 5.2-5.4).

Inputs: `src.stage_a` outputs (records, 8-lead templates, matched pairs). Per pair:
- SR loop: the SR template's VCG (Kors) from QRS onset to offset, resampled to 100 points (3 x 100);
- PVC window: the PVC template's VCG from 10 ms before to 180 ms after the QRS onset (3 x 96);
- A = (X2 X1' + λI)(X1 X1' + λI)^-1, λ = κ tr(X1 X1')/3, X1 -> X2 the two recordings' SR loops;
- d(P, Q) = 1 - cosine similarity of the PVC windows, Q shifted by up to ±10 ms to the best value;
  d̄ = ½ [d(A12 P1, P2) + d(A21 P2, P1)]; secondary: normalised RMS ||P - Q|| / ||Q|| at that shift.
Arms: raw (A = I), own (the patient's SR loops), null (A from a donor's two SR loops: same sex, age
within ±5 years at the first recording, same interval class: same day, 1-365 days, > 365 days;
20 donors drawn with `SEED` among patients of the same stratum, d_null the mean over donors).
Test: d_own < d_raw and d_own < d_null, each a paired one-sided Wilcoxon signed-rank test at 0.025.
Effect: (median d_raw - median d_own) / median d_raw, 95% CI from 2000 patient resamples.

κ is chosen on HEEDB-MGH from KAPPAS by the median d_own in the primary stratum and then fixed
(plan 5.2, item 8 of Section 10); on MGH the H4 results are descriptive (development set).

HEEDB-EUH is the confirmatory test (plan 5.4): κ is read from objects/kappa.json, nothing is tuned,
the run refuses an uncommitted working tree and records the commit. H4 is supported when, in the
primary stratum, both one-sided Wilcoxon p values are < 0.025.

    python -m src.h4 mgh   -> docs/wp2_h4_mgh.md; ARTIFACTS/wp2/mgh_h4.parquet; objects/kappa.json,
                              objects/noise_A.npy (A_own of the primary stratum, both directions)
    python -m src.h4 euh   -> docs/wp2_h4_euh.md (confirmatory); ARTIFACTS/wp2/euh_h4.parquet
"""
import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from src.config import ARTIFACTS, DOCS, HEEDB, KORS, N_BOOT, SEED

KAPPAS = [0.01, 0.03, 0.1, 0.3, 1.0]
N_DONORS = 20
AGE_TOL = 5.0
LOOP_N = 100
WIN = (5, 90)                              # samples at 500 Hz: -10 to +180 ms around the QRS onset
LAG = 5                                    # ±10 ms
OUT = ARTIFACTS / "wp2"
OBJ = OUT / "objects"
SITE_DIR = {"mgh": "I0001", "euh": "I0006"}


def sr_loop(t):
    v = KORS @ t["sr"][:, t["sr_on"]:t["sr_off"] + 1]
    x = np.linspace(0, v.shape[1] - 1, LOOP_N)
    return np.stack([np.interp(x, np.arange(v.shape[1]), a) for a in v])


def pvc_vcg(t):
    return KORS @ t["pvc"], t["pvc_on"]


def window(v, on, s=0):
    return v[:, on - WIN[0] + s:on + WIN[1] + 1 + s]


def fit_a(x1, x2, kappa):
    c = x1 @ x1.T
    lam = kappa * np.trace(c) / 3
    return (x2 @ x1.T + lam * np.eye(3)) @ np.linalg.inv(c + lam * np.eye(3))


def disc(p, qv, q_on):
    """(d, nRMS) of window p against qv shifted to the best cosine within ±LAG."""
    best, best_q = -2.0, None
    for s in range(-LAG, LAG + 1):
        if q_on - WIN[0] + s < 0 or q_on + WIN[1] + 1 + s > qv.shape[1]:
            continue
        q = window(qv, q_on, s)
        c = float((p * q).sum() / np.sqrt((p * p).sum() * (q * q).sum() + 1e-12))
        if c > best:
            best, best_q = c, q
    return 1 - best, float(np.linalg.norm(p - best_q) / (np.linalg.norm(best_q) + 1e-12))


def dbar(a12, a21, e1, e2):
    """Both directions averaged; e = (PVC VCG, onset)."""
    (v1, o1), (v2, o2) = e1, e2
    d12, n12 = disc(a12 @ window(v1, o1), v2, o2)
    d21, n21 = disc(a21 @ window(v2, o2), v1, o1)
    return (d12 + d21) / 2, (n12 + n21) / 2


def demographics(site, files):
    """Sex and age (years) per file stem. Age: AgeAtAcquisition (days), else birth date to acquisition
    time; birth dates are "None" or "redacted" in part of HEEDB (redacted ones usually keep
    AgeAtAcquisition)."""
    sex_col = {"mgh": "SexDSC", "euh": "Sex"}[site]
    m = pd.read_csv(HEEDB / "ECG" / SITE_DIR[site] / "metadata" / "metadata.csv",
                    usecols=["FileName", sex_col, "DateOfBirth", "ECGAcquisitionTime", "AgeAtAcquisition"],
                    dtype=str, engine="pyarrow")
    m["SexDSC"] = m[sex_col].str.capitalize()                       # MGH "Male"/"Female", EUH "MALE"/"FEMALE"
    m["stem"] = m.FileName.str.rsplit("/", n=1).str[-1]
    m = m[m.stem.isin({f.rsplit("/", 1)[-1] for f in files})].drop_duplicates("stem").set_index("stem")
    from_dob = (pd.to_datetime(m.ECGAcquisitionTime, errors="coerce", format="mixed")
                - pd.to_datetime(m.DateOfBirth, errors="coerce", format="mixed")).dt.days / 365.25
    age = (pd.to_numeric(m.AgeAtAcquisition, errors="coerce") / 365.25).fillna(from_dob)
    return m.SexDSC.where(m.SexDSC.isin(["Male", "Female"])), age


def interval_class(days):
    return np.where(days < 1, "same day", np.where(days <= 365, "1-365 days", "> 365 days"))


def run_stratum(pr, rec, tmpl, kappa, sex, age, rng):
    """Per-patient d_raw, d_own, d_null (and nRMS) for one stratum."""
    loops = {k: sr_loop(tmpl[k]) for k in set(pr.i) | set(pr.j)}
    pv = {k: pvc_vcg(tmpl[k]) for k in loops}
    own = [(fit_a(loops[i], loops[j], kappa), fit_a(loops[j], loops[i], kappa)) for i, j in zip(pr.i, pr.j)]
    stem = rec.file.str.rsplit("/", n=1).str[-1]
    s = sex.reindex(stem.loc[pr.i]).to_numpy()
    a = age.reindex(stem.loc[pr.i]).to_numpy()
    cls = interval_class(pr.days.to_numpy())
    rows = []
    eye = np.eye(3)
    for k, (i, j) in enumerate(zip(pr.i, pr.j)):
        raw, raw_n = dbar(eye, eye, pv[i], pv[j])
        o, o_n = dbar(*own[k], pv[i], pv[j])
        ok = np.flatnonzero((s == s[k]) & (np.abs(a - a[k]) <= AGE_TOL) & (cls == cls[k]))
        ok = ok[ok != k]
        if pd.isna(s[k]) or pd.isna(a[k]) or not len(ok):
            nul, nul_n, nd = np.nan, np.nan, 0
        else:
            donors = rng.choice(ok, min(N_DONORS, len(ok)), replace=False)
            vals = [dbar(*own[q], pv[i], pv[j]) for q in donors]
            nul, nul_n, nd = np.mean([v[0] for v in vals]), np.mean([v[1] for v in vals]), len(donors)
        rows.append({"d_raw": raw, "d_own": o, "d_null": nul, "n_raw": raw_n, "n_own": o_n, "n_null": nul_n,
                     "n_donors": nd})
    return pd.DataFrame(rows), own


def summarise(df):
    out = {"n": len(df), "with_donors": int(df.d_null.notna().sum())}
    for c in ("d_raw", "d_own", "d_null", "n_raw", "n_own", "n_null"):
        out[c] = float(np.nanmedian(df[c]))
    out["p_own_raw"] = float(wilcoxon(df.d_own, df.d_raw, alternative="less").pvalue)
    d = df.dropna(subset=["d_null"])
    out["p_own_null"] = float(wilcoxon(d.d_own, d.d_null, alternative="less").pvalue) if len(d) else np.nan
    share = lambda x: (np.median(x.d_raw) - np.median(x.d_own)) / np.median(x.d_raw)
    rng = np.random.default_rng(SEED)
    boots = [share(df.iloc[rng.integers(0, len(df), len(df))]) for _ in range(N_BOOT)]
    out["share"], out["share_ci"] = float(share(df)), (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))
    return out


def clean_commit():
    """Commit of the repository; refuses when tracked files differ from it (confirmatory run)."""
    root = Path(__file__).resolve().parents[1]
    dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
                           capture_output=True, text=True, check=True).stdout.strip()
    if dirty:
        raise SystemExit("confirmatory run needs a committed working tree:\n" + dirty)
    return subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True,
                          check=True).stdout.strip()


def confirmatory(site):
    commit = clean_commit()
    kappa = json.loads((OBJ / "kappa.json").read_text())["kappa"]
    rec = pd.read_parquet(OUT / f"{site}_records.parquet")
    tmpl = np.load(OUT / f"{site}_templates.npy", allow_pickle=True)
    pairs = pd.read_parquet(OUT / f"{site}_pairs.parquet")
    sex, age = demographics(site, rec.file)
    res, frames = {}, []
    for stratum in ("primary (<= 365 days)", "same day", "any interval"):
        pr = pairs[pairs.stratum == stratum].reset_index(drop=True)
        df, _ = run_stratum(pr, rec, tmpl, kappa, sex, age, np.random.default_rng(SEED))
        df["stratum"] = stratum
        frames.append(df)
        res[stratum] = summarise(df)
    pd.concat(frames).to_parquet(OUT / f"{site}_h4.parquet")
    p = res["primary (<= 365 days)"]
    supported = p["p_own_raw"] < 0.025 and p["p_own_null"] < 0.025
    lines = [f"# WP2: H4 confirmatory test on HEEDB-{site.upper()}", "",
             f"Generated by `python -m src.h4 {site}` at commit `{commit}`; κ = {kappa} frozen on HEEDB-MGH",
             "(`objects/kappa.json`). Definitions: module docstring; analysis plan 5.2-5.4.", "",
             "| Stratum | Patients (with donors) | median d_raw | median d_own | median d_null | p own<raw | p own<null | correctable share (95% CI) |",
             "|---|---|---|---|---|---|---|---|"]
    for s, r in res.items():
        lines.append(f"| {s} | {r['n']} ({r['with_donors']}) | {r['d_raw']:.4f} | {r['d_own']:.4f} | {r['d_null']:.4f} | "
                     f"{r['p_own_raw']:.1e} | {r['p_own_null']:.1e} | {r['share']:.3f} ({r['share_ci'][0]:.3f} to {r['share_ci'][1]:.3f}) |")
    lines += ["", "Secondary measure (median normalised RMS): " + "; ".join(
        f"{s}: raw {r['n_raw']:.3f}, own {r['n_own']:.3f}, null {r['n_null']:.3f}" for s, r in res.items()) + ".",
              "", f"**H4 {'supported' if supported else 'not supported'}** (primary stratum: both one-sided "
              f"Wilcoxon p < 0.025; own < raw p = {p['p_own_raw']:.1e}, own < null p = {p['p_own_null']:.1e})."]
    (DOCS / f"wp2_h4_{site}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site", choices=["mgh", "euh"])
    site = ap.parse_args().site
    if site == "euh":
        return confirmatory(site)
    rec = pd.read_parquet(OUT / f"{site}_records.parquet")
    tmpl = np.load(OUT / f"{site}_templates.npy", allow_pickle=True)
    pairs = pd.read_parquet(OUT / f"{site}_pairs.parquet")
    sex, age = demographics(site, rec.file)
    primary = pairs[pairs.stratum.str.startswith("primary")].reset_index(drop=True)
    # κ: median d_own in the primary stratum
    k_tab = []
    for kappa in KAPPAS:
        loops = {k: sr_loop(tmpl[k]) for k in set(primary.i) | set(primary.j)}
        pv = {k: pvc_vcg(tmpl[k]) for k in loops}
        d = [dbar(fit_a(loops[i], loops[j], kappa), fit_a(loops[j], loops[i], kappa), pv[i], pv[j])[0]
             for i, j in zip(primary.i, primary.j)]
        k_tab.append((kappa, float(np.median(d))))
    kappa = min(k_tab, key=lambda t: t[1])[0]
    # shape check beyond the pre-specified grid (reported, not used for the choice)
    check = []
    for kx in (3.0, 10.0, 1e6):
        d = [dbar(fit_a(loops[i], loops[j], kx), fit_a(loops[j], loops[i], kx), pv[i], pv[j])[0]
             for i, j in zip(primary.i, primary.j)]
        check.append((kx, float(np.median(d))))
    OBJ.mkdir(parents=True, exist_ok=True)
    (OBJ / "kappa.json").write_text(json.dumps({"kappa": kappa, "grid": k_tab, "site": site}))
    res, frames = {}, []
    for stratum in ("primary (<= 365 days)", "same day", "any interval"):
        pr = pairs[pairs.stratum == stratum].reset_index(drop=True)
        df, own = run_stratum(pr, rec, tmpl, kappa, sex, age, np.random.default_rng(SEED))
        df["stratum"] = stratum
        frames.append(df)
        res[stratum] = summarise(df)
        if stratum.startswith("primary"):
            np.save(OBJ / "noise_A.npy", np.array([m for pair in own for m in pair]))
    pd.concat(frames).to_parquet(OUT / f"{site}_h4.parquet")
    lines = [f"# WP2: H4 on HEEDB-{site.upper()} (development, descriptive) and κ", "",
             f"Generated by `python -m src.h4 {site}`. Definitions: module docstring; analysis plan 5.2-5.4.", "",
             "## κ (median d_own, primary stratum)", "", "| κ | median d_own |", "|---|---|"]
    lines += [f"| {k} | {v:.4f} |" for k, v in k_tab] + ["", f"**Chosen κ = {kappa}.** Check beyond the grid (not used "
              "for the choice): " + ", ".join(f"κ {k:g}: {v:.4f}" for k, v in check) + "; κ → ∞ gives d_raw.", "",
              "## H4 (MGH is the development set; the confirmatory test is on EUH)", "",
              "| Stratum | Patients (with donors) | median d_raw | median d_own | median d_null | p own<raw | p own<null | correctable share (95% CI) |",
              "|---|---|---|---|---|---|---|---|"]
    for s, r in res.items():
        lines.append(f"| {s} | {r['n']} ({r['with_donors']}) | {r['d_raw']:.4f} | {r['d_own']:.4f} | {r['d_null']:.4f} | "
                     f"{r['p_own_raw']:.1e} | {r['p_own_null']:.1e} | {r['share']:.3f} ({r['share_ci'][0]:.3f} to {r['share_ci'][1]:.3f}) |")
    lines += ["", "Secondary measure (median normalised RMS): " + "; ".join(
        f"{s}: raw {r['n_raw']:.3f}, own {r['n_own']:.3f}, null {r['n_null']:.3f}" for s, r in res.items()) + "."]
    (DOCS / f"wp2_h4_{site}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
