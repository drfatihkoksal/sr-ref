"""Why same-day V2TR discordance differs between sources (exploratory, post hoc, label-free; outside
the registration https://osf.io/wyx28; docs/decisions.md item 36).

Prompted by the registered result: in the same-day stratum, with both ECGs in the V2TR scope, the call
differed in 38.8% of MIMIC-IV-ECG pairs against 14.4% on HEEDB-EUH. Nothing here is confirmatory and no
registered result changes. Inputs are the stored pairs and templates; no pair is rebuilt and the
registered run is not repeated. Cells with fewer than 5 pairs are not shown.

Pairs: the same-day stratum (< 24 h) and, for context, the primary stratum at 1-365 days. Looked at:
- the time between the two ECGs and the discordance by interval;
- how alike the two SR templates are (multi-lead correlation, `stage_a.match_r`), as a proxy for
  whether the electrodes were left in place;
- the distance of the calls from the threshold, the size of the change in V2, and QS-like PVCs in V2
  (R fraction r < 0.05);
- template quality: number of PVC beats, and high-frequency noise of the V2 SR template before the QRS
  (SD of the second difference, µV);
- MIMIC-IV-ECG only: same or different cart, bandwidth and filter setting (HEEDB metadata has no
  device field).
The effect of the interval mix is shown by direct standardisation: the discordance of one source at
the interval mix of another, with a bin-stratified bootstrap 95% CI (2000 resamples).

    python -m src.sameday_explore   -> docs/wp3_sameday_exploratory.md
"""
import numpy as np
import pandas as pd

from src.audit import wilson
from src.config import ARTIFACTS, DOCS, MIMIC, N_BOOT, SEED
from src.criteria_reliability import scope
from src.h4 import OUT
from src.h4_clinical import V2, rfrac
from src.stage_a import match_r

OUT3 = ARTIFACTS / "wp3"
SOURCES = {"HEEDB-MGH": (OUT, "mgh"), "HEEDB-EUH": (OUT, "euh"), "MIMIC-IV-ECG": (OUT3, "mimic")}
SAME_DAY_BINS = [("< 10 min", 0, 10), ("10-60 min", 10, 60), ("1-6 h", 60, 360), ("6-24 h", 360, 1440)]
PRIMARY_BINS = [("< 1 h", 0, 60), ("1-24 h", 60, 1440), ("1-30 days", 1440, 43200), ("31-365 days", 43200, 525700)]
MIN_CELL = 5


def cells(g, bins):
    return [g[g.minutes.between(lo, hi, inclusive="left")].differ.to_numpy(float) for _, lo, hi in bins]


def standard(g, ref, bins):
    """Discordance of g at the interval mix of ref, with a bin-stratified bootstrap 95% CI."""
    c = cells(g, bins)
    if min(len(x) for x in c) == 0:
        return "—"
    w = np.array([len(x) for x in cells(ref, bins)], float)
    w /= w.sum()
    rng = np.random.default_rng(SEED)
    boots = [sum(wi * x[rng.integers(0, len(x), len(x))].mean() for wi, x in zip(w, c)) for _ in range(N_BOOT)]
    return f"{sum(wi * x.mean() for wi, x in zip(w, c)):.1%} ({np.percentile(boots, 2.5):.0%}-{np.percentile(boots, 97.5):.0%})"


def noise(t):
    y = t["sr"][V2, :max(t["sr_on"] - 10, 0)]
    return float(np.std(np.diff(y, 2))) if len(y) >= 5 else np.nan


def table(rec, tmpl, pr):
    n_pvc = rec.n_pvc.to_numpy()
    out = []
    for i, j, days, r in zip(pr.i, pr.j, pr.days, pr.r):
        a, b = tmpl[i], tmpl[j]
        ra, rb = rfrac(a["pvc"], a["pvc_on"], a["pvc_off"])[V2], rfrac(b["pvc"], b["pvc_on"], b["pvc_off"])[V2]
        sa, sb = rfrac(a["sr"], a["sr_on"], a["sr_off"])[V2], rfrac(b["sr"], b["sr_on"], b["sr_off"])[V2]
        va, vb = ra / max(sa, 0.01), rb / max(sb, 0.01)
        out.append({"i": i, "j": j, "minutes": days * 1440, "match": r, "scope": scope(a)[1] and scope(b)[1],
                    "differ": bool((va >= 0.6) != (vb >= 0.6)), "lvot": bool(va >= 0.6), "va": va, "vb": vb, "ra": ra, "rb": rb,
                    "sa": sa, "sb": sb, "sr_r": match_r(a["sr"], b["sr"], a["sr_on"], b["sr_on"]),
                    "beats": int(min(n_pvc[i], n_pvc[j])), "noise": max(noise(a), noise(b))})
    return pd.DataFrame(out)


def rate(g):
    n = len(g)
    if n < MIN_CELL:
        return f"n < {MIN_CELL}"
    k = int(g.differ.sum())
    lo, hi = wilson(k, n)
    return f"{k / n:.1%} ({max(lo, 0.0):.0%}-{hi:.0%}), n {n}"


def med(x, d=2):
    return f"{np.median(x):.{d}f} [{np.percentile(x, 25):.{d}f}, {np.percentile(x, 75):.{d}f}]"


def cart_rows(rec, d):
    mm = pd.read_csv(MIMIC / "machine_measurements.csv", usecols=["study_id", "cart_id", "bandwidth", "filtering"], dtype=str).set_index("study_id")
    stem = rec.file.str.rsplit("/", n=1).str[-1].to_numpy()
    m = mm.reindex(stem)
    out = d.copy()
    for c in ("cart_id", "bandwidth", "filtering"):
        v = m[c].to_numpy()
        out[c] = ["unknown" if pd.isna(v[i]) or pd.isna(v[j]) else ("same" if v[i] == v[j] else "different") for i, j in zip(d.i, d.j)]
    return out


def main():
    data = {}
    for name, (root, key) in SOURCES.items():
        rec = pd.read_parquet(root / f"{key}_records.parquet")
        tmpl = np.load(root / f"{key}_templates.npy", allow_pickle=True)
        pairs = pd.read_parquet(root / f"{key}_pairs.parquet")
        same = table(rec, tmpl, pairs[pairs.stratum == "same day"])
        prim = table(rec, tmpl, pairs[pairs.stratum.str.startswith("primary")])
        data[name] = (rec, same, prim[prim.minutes >= 1440], prim)
    names = list(SOURCES)
    head = "| Measure | " + " | ".join(names) + " |"
    rule = "|---|" + "---|" * len(names)
    lines = ["# Same-day V2TR discordance between sources (exploratory, post hoc, label-free)", "",
             "Generated by `python -m src.sameday_explore`; definitions in the module docstring. Outside the registration",
             "https://osf.io/wyx28. Nothing here is confirmatory. V2TR, both ECGs in the published scope unless stated.", "",
             "## Discordance by interval between the two ECGs", "", head, rule]
    for label, lo, hi in SAME_DAY_BINS:
        lines.append(f"| Same day, {label} | " + " | ".join(rate(data[n][1][data[n][1].scope & data[n][1].minutes.between(lo, hi, inclusive='left')]) for n in names) + " |")
    lines.append("| Same day, all | " + " | ".join(rate(data[n][1][data[n][1].scope]) for n in names) + " |")
    lines.append("| 1-365 days | " + " | ".join(rate(data[n][2][data[n][2].scope]) for n in names) + " |")
    lines += ["", "All pairs, whatever the scope:", "", head, rule]
    for label, lo, hi in SAME_DAY_BINS:
        lines.append(f"| Same day, {label} | " + " | ".join(rate(data[n][1][data[n][1].minutes.between(lo, hi, inclusive='left')]) for n in names) + " |")
    lines.append("| Same day, all | " + " | ".join(rate(data[n][1]) for n in names) + " |")
    lines.append("| 1-365 days | " + " | ".join(rate(data[n][2]) for n in names) + " |")
    s = {n: data[n][1][data[n][1].scope] for n in names}
    a = {n: data[n][1] for n in names}
    p = {n: data[n][3][data[n][3].scope] for n in names}
    lines += ["", "## Primary stratum, both ECGs in scope: share of pairs and discordance by interval", "", head, rule]
    for label, lo, hi in PRIMARY_BINS:
        bin_ = {n: p[n][p[n].minutes.between(lo, hi, inclusive="left")] for n in names}
        lines.append(f"| {label} | " + " | ".join(f"{len(bin_[n]) / len(p[n]):.0%} of pairs; {rate(bin_[n])}" for n in names) + " |")
    lines.append("| All | " + " | ".join(rate(p[n]) for n in names) + " |")
    M, E, X = names
    lines += ["", "## Direct standardisation to another source's interval mix", "",
              "| Pairs | MGH crude | EUH crude | MIMIC crude | MIMIC at the MGH mix | MGH at the MIMIC mix | EUH at the MIMIC mix |",
              "|---|---|---|---|---|---|---|"]
    for label, src, bins in (("Same day, in scope", s, SAME_DAY_BINS), ("Same day, all pairs", a, SAME_DAY_BINS),
                             ("Primary, in scope", p, PRIMARY_BINS), ("Primary, all pairs", {n: data[n][3] for n in names}, PRIMARY_BINS)):
        lines.append(f"| {label} | {src[M].differ.mean():.1%} | {src[E].differ.mean():.1%} | {src[X].differ.mean():.1%} | "
                     f"{standard(src[X], src[M], bins)} | {standard(src[M], src[X], bins)} | {standard(src[E], src[X], bins)} |")
    row = lambda label, f, src=s: f"| {label} | " + " | ".join(f(src[n]) if len(src[n]) >= MIN_CELL else f"n < {MIN_CELL}" for n in names) + " |"
    near = lambda g: g.va.between(0.4, 0.9) | g.vb.between(0.4, 0.9)
    qs1 = lambda g: (g.ra < 0.05) != (g.rb < 0.05)
    lines += ["", "## Same-day pairs in scope: what differs between the sources", "", head, rule,
              row("Pairs", lambda g: str(len(g))),
              row("Minutes between ECGs, median [IQR]", lambda g: med(g.minutes, 0)),
              row("Under 60 minutes apart", lambda g: f"{(g.minutes < 60).mean():.0%}"),
              row("SR templates alike, r, median [IQR]", lambda g: med(g.sr_r, 3)),
              row("PVC match r, median [IQR]", lambda g: med(g.match, 3)),
              row("LVOT call on the first ECG", lambda g: f"{g.lvot.mean():.0%}"),
              row("A V2TR value within 0.40-0.90", lambda g: f"{near(g).mean():.0%}"),
              row("PVC r in V2, first ECG, median [IQR]", lambda g: med(g.ra, 2)),
              row("SR r in V2, first ECG, median [IQR]", lambda g: med(g.sa, 2)),
              row("|change| of PVC r in V2, median [IQR]", lambda g: med((g.ra - g.rb).abs(), 3)),
              row("|change| of SR r in V2, median [IQR]", lambda g: med((g.sa - g.sb).abs(), 3)),
              row("QS-like PVC in V2 on exactly one ECG", lambda g: f"{qs1(g).mean():.0%}"),
              row("Fewest PVC beats per template, median [IQR]", lambda g: med(g.beats, 0)),
              row("Single-beat PVC template on an ECG", lambda g: f"{(g.beats == 1).mean():.0%}"),
              row("V2 SR template noise, larger of the two, µV", lambda g: med(g.noise.dropna(), 1)),
              "", "## Same-day discordance by features (all same-day pairs, more numbers than in scope)", "", head, rule]
    pooled = np.nanmedian(pd.concat([a[n].sr_r for n in names]))
    noise_cut = np.nanmedian(pd.concat([a[n].noise for n in names]))
    lines += [f"| SR templates alike (r ≥ pooled median {pooled:.3f}) | " + " | ".join(rate(a[n][a[n].sr_r >= pooled]) for n in names) + " |",
              "| SR templates less alike | " + " | ".join(rate(a[n][a[n].sr_r < pooled]) for n in names) + " |",
              f"| V2 noise below the pooled median ({noise_cut:.1f} µV) | " + " | ".join(rate(a[n][a[n].noise < noise_cut]) for n in names) + " |",
              "| V2 noise at or above it | " + " | ".join(rate(a[n][a[n].noise >= noise_cut]) for n in names) + " |",
              "| At least 2 PVC beats on both ECGs | " + " | ".join(rate(a[n][a[n].beats >= 2]) for n in names) + " |",
              "| A single-beat PVC template | " + " | ".join(rate(a[n][a[n].beats == 1]) for n in names) + " |",
              "| PVC match r ≥ 0.95 | " + " | ".join(rate(a[n][a[n].match >= 0.95]) for n in names) + " |",
              "| PVC match r < 0.95 | " + " | ".join(rate(a[n][a[n].match < 0.95]) for n in names) + " |"]
    m = cart_rows(data["MIMIC-IV-ECG"][0], a["MIMIC-IV-ECG"])
    lines += ["", "## MIMIC-IV-ECG: same-day pairs by acquisition setting", "", "| Setting | Same | Different | Unknown |", "|---|---|---|---|"]
    for c, label in (("cart_id", "Cart"), ("bandwidth", "Bandwidth"), ("filtering", "Filter")):
        lines.append(f"| {label} | " + " | ".join(rate(m[m[c] == v]) for v in ("same", "different", "unknown")) + " |")
    for label, mask in (("cart, < 60 min apart", m.minutes < 60), ("cart, ≥ 60 min apart", m.minutes >= 60)):
        lines.append(f"| By {label} | " + " | ".join(rate(m[mask & (m.cart_id == v)]) for v in ("same", "different", "unknown")) + " |")
    (DOCS / "wp3_sameday_exploratory.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
