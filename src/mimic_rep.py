"""MIMIC-IV-ECG replication of the criteria test-retest result and of H4 (preregistered on OSF;
docs/decisions.md items 33, 34).

The feasibility step reads only the machine statements, the acquisition times and the MIMIC-IV
`patients` table and writes counts. The run step reads the waveforms of the sample once, after
registration; it refuses code or constants that differ from prereg/mimic_manifest.json, an
uncommitted working tree, and a second start. Before registration, `dryrun` checks the pair and
analysis code on HEEDB-EUH, and `smoke` checks the record pipeline on MIMIC records of excluded
subjects that were already read for K1-T (counts only, no pair, no criterion).

Sequence (one-sided level 0.025 each; testing stops at the first hypothesis not met): R1 V2TR
discordance f above 0.10 (lower Wilson 95% limit); R2 H4 (`src.h4`); R3 f above 0.18; R4 own 8 x 8
transform below shared SR (upper limit of the 95% CI of the paired difference below 0). Primary
stratum; for f and R4, pairs with both ECGs in the published V2TR scope.

Candidates: the K1-T statement rules of `src.audit` (INC_MIMIC, EXC_MIMIC), all records of a subject.
Prior contact: subjects among the first 400 of the K1-T candidate order (284 processed for K1-T, 400
for the yield check) and subjects of the 500 ECGs of the QRS scale (`src.qrs_scale`).
Projection: subjects with at least two candidate records times the yield of HEEDB-MGH and HEEDB-EUH
(primary pairs and V2 transition ratio pairs with both ECGs in scope, per patient with at least two
500-Hz candidate records). Power: exact binomial probability that the lower Wilson 95% bound exceeds
each threshold (0.10: a single-ECG accuracy of at least 95%, as implied by sensitivity 95% and
specificity 100%, needs a discordance of at most 10%; 0.18: 91% accuracy, 1 - f/2); for the transform
arms, a one-sided 0.025 normal approximation from the HEEDB paired SDs.

    python -m src.mimic_rep feasibility            -> docs/wp3_mimic_feasibility.md
    python -m src.mimic_rep dryrun                 -> docs/wp3_mimic_dryrun_euh.md
    python -m src.mimic_rep smoke
    python -m src.mimic_rep lock                   -> prereg/mimic_manifest.json (prints its SHA-256)
    python -m src.mimic_rep run --registered URL   -> docs/wp3_mimic_replication.md; ARTIFACTS/wp3
"""
import argparse
import hashlib
import json
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom, norm

from src.audit import EXC_MIMIC, INC_MIMIC, mimic_candidates, wilson
from src.config import ARTIFACTS, BEAT_CLASSIFIER, DOCS, MIMIC, MIMIC_IV, SEED
from src.criteria_controls import boot_diff, diff_ci, rows
from src.criteria_reliability import kappa as cohen_kappa
from src.h4 import OBJ, OUT, clean_commit, demographics as heedb_demographics, run_stratum, summarise
from src.qrs_scale import mimic_sample
from src.records import load, process
from src.signal import IDX8
from src.stage_a import pairs

ROOT = Path(__file__).resolve().parents[1]
OUT3 = ARTIFACTS / "wp3"
MANIFEST = ROOT / "prereg" / "mimic_manifest.json"
PRIMARY = "primary (<= 365 days)"
STRATA = (PRIMARY, "same day", "any interval")
CODE = ("audit", "beatnet", "beats", "config", "criteria_controls", "criteria_reliability", "h4", "h4_clinical",
        "h4_leads", "mimic_rep", "qrs_scale", "records", "signal", "stage_a")

N_PRIOR = 400           # MIMIC candidates read before registration, in the K1-T order
F0S = (0.10, 0.18)      # discordance above which a single ECG cannot reach 95% / 91% accuracy (1 - f/2)
FS = (0.20, 0.22, 0.24)


def candidates():
    cols = [f"report_{i}" for i in range(18)]
    mm = pd.read_csv(MIMIC / "machine_measurements.csv", usecols=["subject_id", "study_id"] + cols, dtype=str)
    rep = mm[cols].fillna("").agg(" ".join, axis=1).str.lower()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ok = rep.str.contains(INC_MIMIC, regex=True) & ~rep.str.contains(EXC_MIMIC, regex=True)
    rl = pd.read_csv(MIMIC / "record_list.csv", dtype=str)[["study_id", "path", "ecg_time"]]
    df = mm.loc[ok, ["subject_id", "study_id"]].merge(rl, on="study_id")
    df["time"] = pd.to_datetime(df.ecg_time)
    return df


def prior_subjects():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        k1t = set(mimic_candidates().patient.iloc[:N_PRIOR])
        qrs = set(mimic_sample().study_id.astype(str))
    rl = pd.read_csv(MIMIC / "record_list.csv", usecols=["subject_id", "study_id"], dtype=str)
    return k1t, set(rl.subject_id[rl.study_id.isin(qrs)])


def within_year(g):
    t = np.sort(g.to_numpy())
    return len(t) >= 2 and (np.diff(t) <= np.timedelta64(365, "D")).any()


def power(n, f, f0):
    return float(sum(binom.pmf(k, n, f) for k in range(1, n + 1) if wilson(k, n)[0] > f0))


def yields():
    out = {}
    for site in ("mgh", "euh"):
        rec = pd.read_parquet(OUT / f"{site}_records.parquet", columns=["patient"])
        pr = pd.read_parquet(OUT / f"{site}_pairs.parquet", columns=["stratum"])
        c = pd.read_parquet(OUT / f"{site}_criteria_controls.parquet")
        both = c[c.v2tr_both]
        diff = both.v2tr_own - both.v2tr_shared
        out[site] = {"patients": rec.patient.nunique(), "primary": int(pr.stratum.str.startswith("primary").sum()),
                     "both": len(both), "effect": float(diff.mean()), "sd": float(diff.std(ddof=1))}
    return out


def feasibility():
    c = candidates()
    n_rec = c.groupby("subject_id").size()
    two = c[c.subject_id.isin(n_rec[n_rec >= 2].index)]
    near = two.groupby("subject_id").time.apply(within_year)
    k1t, qrs = prior_subjects()
    prior = k1t | qrs
    keep = set(n_rec[n_rec >= 2].index) - prior
    keep_near = set(near[near].index) - prior
    pat = pd.read_csv(MIMIC_IV / "hosp" / "patients.csv.gz", usecols=["subject_id"], dtype=str)
    y = yields()
    z = norm.ppf(0.975)
    lines = ["# MIMIC-IV-ECG replication: metadata feasibility (counts only)", "",
             "Generated by `python -m src.mimic_rep feasibility`; definitions in the module docstring. No waveform",
             "was read, no pair was built and no criterion was measured.", "",
             "| Step | Count |", "|---|---|",
             f"| Records with a PVC statement and no excluded statement | {len(c)} |",
             f"| Subjects with such a record | {len(n_rec)} |",
             f"| Subjects with at least 2 such records | {int((n_rec >= 2).sum())} |",
             f"| Of these, at least 2 records within 365 days | {int(near.sum())} |",
             f"| Subjects with prior waveform contact: K1-T and yield check | {len(k1t)} |",
             f"| Subjects with prior waveform contact: QRS scale | {len(qrs)} |",
             f"| Subjects with at least 2 records, no prior contact | {len(keep)} |",
             f"| Of these, at least 2 records within 365 days | {len(keep_near)} |",
             f"| Of these, in the MIMIC-IV `patients` table (sex and age for donors) | {len(keep_near & set(pat.subject_id))} |", "",
             "## Projected yield (HEEDB per-patient yields applied to subjects with at least 2 records, no prior contact)", "",
             "| Yield from | HEEDB patients | Primary pairs | V2TR both in scope | Projected primary pairs | Projected V2TR pairs |",
             "|---|---|---|---|---|---|"]
    proj = {}
    for s, v in y.items():
        proj[s] = round(len(keep) * v["both"] / v["patients"])
        lines.append(f"| {s.upper()} | {v['patients']} | {v['primary']} | {v['both']} | "
                     f"{round(len(keep) * v['primary'] / v['patients'])} | {proj[s]} |")
    for f0 in F0S:
        lines += ["", f"## Power: lower Wilson 95% bound of the V2TR discordance > {f0:.2f}", "",
                  "| V2TR pairs | " + " | ".join(f"true f {f:.2f}" for f in FS) + " |",
                  "|---|" + "---|" * len(FS)]
        for n in sorted({*proj.values(), 200, 300, 500, 1000}):
            lines.append(f"| {n} | " + " | ".join(f"{power(n, f, f0):.2f}" for f in FS) + " |")
    lines += ["", "## Power of the transform comparison (own transform vs shared SR, V2TR, one-sided 0.025)", "",
              "| HEEDB source of effect and SD | Effect, points | Paired SD, points | Power at the projected pairs (MGH yield / EUH yield) |",
              "|---|---|---|---|"]
    for s, v in y.items():
        pw = [norm.cdf(abs(v["effect"]) / (v["sd"] / np.sqrt(n)) - z) for n in proj.values()]
        lines.append(f"| {s.upper()} | {v['effect'] * 100:+.1f} | {v['sd'] * 100:.1f} | " + " / ".join(f"{p:.2f}" for p in pw) + " |")
    (DOCS / "wp3_mimic_feasibility.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def manifest():
    files = {f"src/{m}.py": sha(ROOT / "src" / f"{m}.py") for m in CODE}
    files["objects/kappa.json"] = sha(OBJ / "kappa.json")
    files["beat_classifier"] = sha(BEAT_CLASSIFIER)
    files["beat_classifier_card"] = sha(BEAT_CLASSIFIER.parent / "model_card.json")
    return files


def lock():
    MANIFEST.write_text(json.dumps({"files": manifest()}, indent=1, sort_keys=True) + "\n")
    print(f"SHA-256 of {MANIFEST.relative_to(ROOT)}: {sha(MANIFEST)}")


def verify():
    locked, now = json.loads(MANIFEST.read_text())["files"], manifest()
    bad = sorted(k for k in set(locked) | set(now) if locked.get(k) != now.get(k))
    if bad:
        raise SystemExit("code or constants differ from the lock: " + ", ".join(bad))


def sample():
    """Candidate records of subjects with no prior waveform contact and at least two candidate records."""
    c = candidates().rename(columns={"subject_id": "patient"})
    k1t, qrs = prior_subjects()
    c = c[~c.patient.isin(k1t | qrs)]
    c = c[c.groupby("patient").path.transform("size") >= 2].reset_index(drop=True)
    c["file"] = [str(MIMIC / p) for p in c.path]
    return c


def _process(f):
    try:
        x, fs = load(f)
    except Exception:
        return {"reason": "unreadable"}, None
    if x is None or np.isnan(x).any():
        return {"reason": "missing lead or samples"}, None
    if fs != 500:
        return {"reason": "not 500 Hz"}, None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = process(x, fs, "mimic")
    row = {"reason": r["reason"] or ("OT-like" if r["ot_like"] else "eligible, not OT-like"),
           "n_pvc": r["pvc"]["n"] if r.get("pvc") else 0, "n_sr": r["sr"]["n"] if r.get("sr") else 0}
    if not (r["eligible"] and r["ot_like"]):
        return row, None
    return row, {"pvc": r["pvc"]["T"][IDX8].astype(np.float32), "sr": r["sr"]["T"][IDX8].astype(np.float32),
                 "pvc_on": r["pvc"]["on"], "pvc_off": r["pvc"]["off"], "sr_on": r["sr"]["on"], "sr_off": r["sr"]["off"]}


def records(c):
    """Record table and templates for candidates c (columns patient, time, file)."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    with ProcessPoolExecutor(max_workers=max(os.cpu_count() - 4, 1)) as ex:
        res = list(ex.map(_process, c.file, chunksize=50))
    rec = pd.concat([c[["patient", "time", "file"]].reset_index(drop=True), pd.DataFrame([r for r, _ in res])], axis=1)
    return rec, {i: t for i, (_, t) in enumerate(res) if t is not None}


def demographics(rec):
    """Sex and age (years) per file stem; age = anchor_age + year of the ECG - anchor_year."""
    pat = pd.read_csv(MIMIC_IV / "hosp" / "patients.csv.gz", usecols=["subject_id", "gender", "anchor_age", "anchor_year"],
                      dtype={"subject_id": str}).set_index("subject_id")
    p = pat.reindex(rec.patient)
    stem = rec.file.str.rsplit("/", n=1).str[-1].to_numpy()
    return (pd.Series(p.gender.to_numpy(), index=stem),
            pd.Series(p.anchor_age.to_numpy() + rec.time.dt.year.to_numpy() - p.anchor_year.to_numpy(), index=stem))


def desc(name, g, c):
    if not len(g):
        return f"| {name} | 0 | — | — | — |"
    x, y = g[f"{c}_1"].to_numpy(bool), g[f"{c}_2"].to_numpy(bool)
    k, n = int((x != y).sum()), len(g)
    lo, hi = wilson(k, n)
    return f"| {name} | {n} | {k / n:.1%} ({lo:.1%} to {hi:.1%}) | {cohen_kappa(x, y):.2f} | ≤ {1 - k / n / 2:.1%} |"


def analyse(rec, tmpl, pr, sex, age):
    """Decisions and descriptive tables as report lines, and the per-pair table of the primary stratum."""
    kappa = json.loads((OBJ / "kappa.json").read_text())["kappa"]
    tab = {s: pr[pr.stratum == s].reset_index(drop=True) for s in STRATA}
    d = rows(rec, tmpl, tab[PRIMARY], sex, age)
    both = d[d.v2tr_both]
    k, n = int(both.v2tr_raw.sum()), len(both)
    lo, hi = wilson(k, n) if n else (0.0, 1.0)
    h = {s: summarise(run_stratum(t, rec, tmpl, kappa, sex, age, np.random.default_rng(SEED))[0]) for s, t in tab.items() if len(t)}
    p = h[PRIMARY]
    m, b_lo, b_hi = boot_diff(both.v2tr_own, both.v2tr_shared) if n else (np.nan, np.nan, np.nan)
    seq = [("R1", f"V2TR discordance above 0.10: lower Wilson limit {lo:.3f}", lo > 0.10),
           ("R2", f"H4: own < raw p = {p['p_own_raw']:.1e}, own < null p = {p['p_own_null']:.1e}",
            p["p_own_raw"] < 0.025 and p["p_own_null"] < 0.025),
           ("R3", f"V2TR discordance above 0.18: lower Wilson limit {lo:.3f}", lo > 0.18),
           ("R4", f"Own transform minus shared SR {m * 100:+.1f} points (95% CI {b_lo * 100:+.1f} to {b_hi * 100:+.1f}); "
                  "upper limit below 0", bool(b_hi < 0))]
    lines = [f"Primary stratum: {len(d)} pairs; V2TR with both ECGs in scope: {n} pairs; discordance "
             f"{k / max(n, 1):.1%} (Wilson 95% CI {lo:.1%} to {hi:.1%}); accuracy ceiling ≤ {1 - k / max(n, 1) / 2:.1%}.", "",
             "| Hypothesis | Criterion | Status |", "|---|---|---|"]
    going = True
    for name, text, met in seq:
        lines.append(f"| {name} | {text} | {('met' if met else 'not met') if going else 'not tested: sequence stopped (descriptive)'} |")
        going = going and met
    conclusion = ("the V2TR call differed in more than 10% of pairs; a single-ECG accuracy of at least 95% cannot hold in "
                  "this population if the origin is fixed" if lo > 0.10 else
                  "the HEEDB finding did not replicate in this dataset" if hi <= 0.10 else
                  "inconclusive; the estimate and its CI are given")
    lines += ["", f"R1 conclusion fixed in advance: {conclusion}.", "", "## H4 (R2)", "",
              "| Stratum | Patients (with donors) | median d_raw | median d_own | median d_null | p own<raw | p own<null | correctable share (95% CI) |",
              "|---|---|---|---|---|---|---|---|"]
    for s, r in h.items():
        lines.append(f"| {s} | {r['n']} ({r['with_donors']}) | {r['d_raw']:.4f} | {r['d_own']:.4f} | {r['d_null']:.4f} | "
                     f"{r['p_own_raw']:.1e} | {r['p_own_null']:.1e} | {r['share']:.3f} ({r['share_ci'][0]:.3f} to {r['share_ci'][1]:.3f}) |")
    extra = {s: rows(rec, tmpl, tab[s], sex, age) if len(tab[s]) else d.iloc[:0] for s in ("same day", "any interval")}
    for c, label in (("v2tr", "V2TR (≥ 0.60)"), ("tzi", "TZ index (< 0)")):
        g = d[d[f"{c}_both"]]
        lines += ["", f"## {label}: calls (descriptive except R1, R3)", "",
                  "| Subset | Pairs | Call differs (Wilson 95% CI) | κ | Ceiling |", "|---|---|---|---|---|",
                  desc("Primary, both ECGs in scope", g, c), desc("Primary, first ECG in scope", d[d[f"{c}_first"]], c),
                  desc("Both in scope, 1-365 days", g[g.days >= 1], c), desc("Both in scope, match r ≥ 0.95", g[g.r >= 0.95], c)]
        for s, t in extra.items():
            lines.append(desc(f"{s.capitalize()} stratum, both in scope", t[t[f"{c}_both"]], c))
        gd = g.dropna(subset=[f"{c}_donor"])
        lines += ["", "| Arm (primary, both in scope) | Percent of pairs |", "|---|---|"]
        for arm, nm in (("raw", "a raw"), ("shared", "b shared SR"), ("own", "c own 8 x 8 transform"), ("sronly", "e SR only")):
            lines.append(f"| {nm} | {g[f'{c}_{arm}'].mean():.1%} |" if len(g) else f"| {nm} | — |")
        lines.append(f"| d donor transform (n {len(gd)}) | {gd[f'{c}_donor'].mean():.1%} |" if len(gd) else "| d donor transform | — |")
        if len(g):
            lines += ["", f"a − b {diff_ci(g[f'{c}_raw'], g[f'{c}_shared'])} points; c − b {diff_ci(g[f'{c}_own'], g[f'{c}_shared'])} points"
                      + (f"; c − d {diff_ci(gd[f'{c}_own'], gd[f'{c}_donor'])} points." if len(gd) else ".")]
    return lines, d


def dryrun():
    """Same pair and analysis code on HEEDB-EUH, where the results are known."""
    rec = pd.read_parquet(OUT / "euh_records.parquet")
    tmpl = np.load(OUT / "euh_templates.npy", allow_pickle=True)
    stored = pd.read_parquet(OUT / "euh_pairs.parquet")
    pr = pairs(rec, tmpl)
    cols = ["stratum", "i", "j"]
    key = lambda t: t[cols].sort_values(cols).reset_index(drop=True)
    if not key(pr).equals(key(stored)):
        raise SystemExit("pairs differ from the stored HEEDB-EUH pairs")
    lines, _ = analyse(rec, tmpl, pr, *heedb_demographics("euh", rec.file))
    head = ["# MIMIC-IV-ECG replication code: dry run on HEEDB-EUH", "",
            "Generated by `python -m src.mimic_rep dryrun`. Pairs re-derived and equal to the stored EUH pairs. This is a",
            "code check on known results, not a test of any hypothesis.", ""]
    (DOCS / "wp3_mimic_dryrun_euh.md").write_text("\n".join(head + lines) + "\n")
    print("\n".join(head + lines))


def smoke():
    """Record pipeline on MIMIC records already read for K1-T (their subjects are excluded); counts only."""
    c = mimic_candidates().iloc[:30]
    rl = pd.read_csv(MIMIC / "record_list.csv", usecols=["study_id", "ecg_time"], dtype=str)
    c = c.merge(rl, on="study_id")
    c["time"] = pd.to_datetime(c.ecg_time)
    rec, _ = records(c)
    sex, age = demographics(rec)
    print(rec.reason.value_counts().to_string())
    print(f"sex known {sex.notna().mean():.0%}; age known {age.notna().mean():.0%}, all 18-100: {bool(age.between(18, 100).all())}")


def run(url):
    if not url.startswith("https://osf.io/"):
        raise SystemExit("--registered needs the OSF registration URL")
    verify()
    commit = clean_commit()
    OUT3.mkdir(parents=True, exist_ok=True)
    marker = OUT3 / "run_started.json"
    if marker.exists():
        raise SystemExit("the registered analysis has already started; it runs once (plan: report any coding error)")
    marker.write_text(json.dumps({"url": url, "commit": commit, "manifest": sha(MANIFEST), "started": datetime.now().isoformat()}))
    c = sample()
    rec, tmpl = records(c)
    rec.to_parquet(OUT3 / "mimic_records.parquet")
    np.save(OUT3 / "mimic_templates.npy", np.array([tmpl.get(i) for i in range(len(rec))], dtype=object), allow_pickle=True)
    pr = pairs(rec, tmpl)
    pr.to_parquet(OUT3 / "mimic_pairs.parquet")
    lines, d = analyse(rec, tmpl, pr, *demographics(rec))
    d.to_parquet(OUT3 / "mimic_criteria.parquet")
    ot = rec[rec.reason == "OT-like"].groupby("patient").size()
    head = ["# MIMIC-IV-ECG replication (preregistered)", "",
            f"Registration: {url}. Generated by `python -m src.mimic_rep run` at commit `{commit}`; manifest SHA-256 "
            f"`{sha(MANIFEST)}`.", "",
            f"Records {len(rec)} of {rec.patient.nunique()} subjects with at least two candidate records and no prior waveform "
            f"contact; subjects with at least two eligible OT-like records: {int((ot >= 2).sum())}.", "",
            "| Outcome | Records |", "|---|---|"] + [f"| {k} | {v} |" for k, v in rec.reason.value_counts().items()] + [
            "", "| Stratum | Pairs |", "|---|---|"] + [f"| {s} | {int((pr.stratum == s).sum()) if len(pr) else 0} |" for s in STRATA] + [
            "", "## Decisions", ""]
    (DOCS / "wp3_mimic_replication.md").write_text("\n".join(head + lines) + "\n")
    print("\n".join(head + lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["feasibility", "dryrun", "smoke", "lock", "run"])
    ap.add_argument("--registered", default="")
    a = ap.parse_args()
    if a.step == "run":
        run(a.registered)
    else:
        {"feasibility": feasibility, "dryrun": dryrun, "smoke": smoke, "lock": lock}[a.step]()


if __name__ == "__main__":
    main()
