"""Controls for the test-retest result of the SR-referenced criteria (post hoc, label-free;
docs/decisions.md item 32). Primary-stratum pairs of `src.stage_a`.

Scope, as in `src.criteria_reliability`, either from the first ECG only or from both ECGs. The
accuracy ceiling (1 - f/2) needs both ECGs in scope, because it compares two ECGs to which the
criterion applies.

Arms, each giving the share of pairs whose call differs (mean of the two directions from arm b on):
a. raw: each ECG with its own PVC and its own SR;
b. shared SR: PVC of one ECG, SR of the other, against the other ECG's own PVC and SR, no transform
   (this removes the between-ECG change of the SR denominator only);
c. own transform: as b, with the PVC first mapped by the 8 x 8 own-SR transform (`h4_leads.fit8`,
   κ8 = 1.0), the corrected arm reported in docs/wp2_criteria_multifocal.md;
d. donor transform: as c, with the transform of a donor patient (same sex, age within ±5 years, same
   interval class; up to 20 donors drawn with `SEED`).
e. SR only: the same PVC measured with its own SR against the other ECG's SR (the call change that
   the between-ECG change of the SR denominator causes on its own).
The effect of the patient's own SR is c against b (beyond a shared denominator) and c against d
(beyond any transform). a against b is the part of the raw discordance that a shared denominator
removes. Differences with 95% CIs from 2000 patient resamples.

    python -m src.criteria_controls mgh euh   -> docs/wp2_criteria_controls.md
"""
import sys

import numpy as np
import pandas as pd

from src.config import DOCS, N_BOOT, SEED
from src.criteria_reliability import scope
from src.h4 import AGE_TOL, N_DONORS, OUT, demographics, interval_class
from src.h4_clinical import V, V2, rfrac, tz
from src.h4_leads import fit8, sr_loop8

KAPPA8 = 1.0


def calls(rp, rs):
    """(V2TR LVOT call, TZ index LVOT call) from the R fractions of a PVC and an SR template. Plain
    bools: the arms add two comparisons, and NumPy bools add as logical OR."""
    return bool(rp[V2] / max(rs[V2], 0.01) >= 0.60), bool(tz(rp[V]) - tz(rs[V]) < 0)


def pair_rows(site):
    rec = pd.read_parquet(OUT / f"{site}_records.parquet")
    tmpl = np.load(OUT / f"{site}_templates.npy", allow_pickle=True)
    pairs = pd.read_parquet(OUT / f"{site}_pairs.parquet")
    pr = pairs[pairs.stratum.str.startswith("primary")].reset_index(drop=True)
    return rows(rec, tmpl, pr, *demographics(site, rec.file))


def rows(rec, tmpl, pr, sex, age):
    """One row per pair of `pr` (default index): scope flags, the call of each ECG and the arms."""
    stem = rec.file.str.rsplit("/", n=1).str[-1]
    s, a = sex.reindex(stem.loc[pr.i]).to_numpy(), age.reindex(stem.loc[pr.i]).to_numpy()
    cls = interval_class(pr.days.to_numpy())
    loops = {k: sr_loop8(tmpl[k]) for k in set(pr.i) | set(pr.j)}
    own = [(fit8(loops[i], loops[j], KAPPA8), fit8(loops[j], loops[i], KAPPA8)) for i, j in zip(pr.i, pr.j)]
    rng = np.random.default_rng(SEED)
    pvc_r = lambda t, A=None: rfrac(t["pvc"] if A is None else A @ t["pvc"], t["pvc_on"], t["pvc_off"])
    sr_r = {k: rfrac(tmpl[k]["sr"], tmpl[k]["sr_on"], tmpl[k]["sr_off"]) for k in loops}
    out = []
    for k, (i, j) in enumerate(zip(pr.i, pr.j)):
        t1, t2 = tmpl[i], tmpl[j]
        s1, s2 = scope(t1), scope(t2)
        r1, r2 = pvc_r(t1), pvc_r(t2)
        ref1, ref2 = calls(r1, sr_r[i]), calls(r2, sr_r[j])
        shared = (calls(r1, sr_r[j]), calls(r2, sr_r[i]))
        cor = (calls(pvc_r(t1, own[k][0]), sr_r[j]), calls(pvc_r(t2, own[k][1]), sr_r[i]))
        ok = np.flatnonzero((s == s[k]) & (np.abs(a - a[k]) <= AGE_TOL) & (cls == cls[k]))
        ok = ok[ok != k]
        don = [np.nan, np.nan]
        if not (pd.isna(s[k]) or pd.isna(a[k]) or not len(ok)):
            flips = [[], []]
            for q in rng.choice(ok, min(N_DONORS, len(ok)), replace=False):
                d1, d2 = calls(pvc_r(t1, own[q][0]), sr_r[j]), calls(pvc_r(t2, own[q][1]), sr_r[i])
                for c in (0, 1):
                    flips[c].append(((d1[c] != ref2[c]) + (d2[c] != ref1[c])) / 2)
            don = [float(np.mean(flips[0])), float(np.mean(flips[1]))]
        row = {"tzi_first": s1[0], "tzi_both": s1[0] and s2[0], "v2tr_first": s1[1], "v2tr_both": s1[1] and s2[1],
               "days": pr.days[k], "r": pr.r[k]}
        for c, name in enumerate(("v2tr", "tzi")):
            row[f"{name}_1"], row[f"{name}_2"] = ref1[c], ref2[c]
            row[f"{name}_raw"] = float(ref1[c] != ref2[c])
            row[f"{name}_shared"] = ((shared[0][c] != ref2[c]) + (shared[1][c] != ref1[c])) / 2
            row[f"{name}_sronly"] = ((shared[0][c] != ref1[c]) + (shared[1][c] != ref2[c])) / 2
            row[f"{name}_own"] = ((cor[0][c] != ref2[c]) + (cor[1][c] != ref1[c])) / 2
            row[f"{name}_donor"] = don[c]
        out.append(row)
    return pd.DataFrame(out)


def boot_diff(a, b):
    """Mean of the paired difference a - b and its percentile 95% CI (patient resamples)."""
    rng = np.random.default_rng(SEED)
    d = (a - b).to_numpy()
    boots = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(N_BOOT)]
    return d.mean(), np.percentile(boots, 2.5), np.percentile(boots, 97.5)


def diff_ci(a, b):
    m, lo, hi = boot_diff(a, b)
    return f"{m * 100:+.1f} ({lo * 100:+.1f} to {hi * 100:+.1f})"


def main():
    sites = sys.argv[1:] or ["mgh", "euh"]
    lines = ["# Controls of the criteria test-retest result (post hoc, label-free)", "",
             "Generated by `python -m src.criteria_controls`; definitions in the module docstring;",
             "`docs/decisions.md` item 32. Primary-stratum pairs. Percent of pairs whose call differs.", ""]
    for site in sites:
        d = pair_rows(site)
        d.to_parquet(OUT / f"{site}_criteria_controls.parquet")
        lines += [f"## HEEDB-{site.upper()}" if site != "mimic" else "## MIMIC-IV-ECG", "",
                  "| Criterion, scope | Pairs | a raw | ceiling 1 − f/2 | b shared SR | c own transform | d donor transform | e SR only | a − b, points (95% CI) | c − b, points (95% CI) | c − d, points (95% CI) |",
                  "|---|---|---|---|---|---|---|---|---|---|---|"]
        for c, label in (("v2tr", "V2TR"), ("tzi", "TZ index")):
            for sc, sl in (("first", "first ECG in scope"), ("both", "both ECGs in scope")):
                g = d[d[f"{c}_{sc}"]]
                gd = g.dropna(subset=[f"{c}_donor"])
                f = g[f"{c}_raw"].mean()
                lines.append(f"| {label}, {sl} | {len(g)} | {f:.1%} | ≤ {1 - f / 2:.1%} | {g[f'{c}_shared'].mean():.1%} | "
                             f"{g[f'{c}_own'].mean():.1%} | {gd[f'{c}_donor'].mean():.1%} | {g[f'{c}_sronly'].mean():.1%} | "
                             f"{diff_ci(g[f'{c}_raw'], g[f'{c}_shared'])} | "
                             f"{diff_ci(g[f'{c}_own'], g[f'{c}_shared'])} | {diff_ci(gd[f'{c}_own'], gd[f'{c}_donor'])} |")
        lines.append("")
    (DOCS / "wp2_criteria_controls.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
