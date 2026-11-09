"""Table 1 and the single file of the numbers used in the paper. Aggregates only.

Inputs are the stored row-level outputs and, for H4, the generated reports. Nothing is estimated anew: the
registered MIMIC-IV-ECG run is not repeated. Sources: criteria (`src.criteria_controls`, `src.mimic_rep`),
intervals (`src.sameday_explore`, exploratory), H4 (`docs/wp2_h4_*.md`, `docs/wp3_mimic_replication.md`),
reader measurement (`src.v2_readers`, exploratory), published sinus rule (`src.betensky_rule`, exploratory), sex
(`src.sex_discordance`, exploratory).

Flow counts group the record outcomes of the cohort step. Table 1 describes the primary-stratum pairs, and those with both ECGs in the V2TR scope. Age and sex are
those at the first ECG of the pair. QRS durations are the first ECG's templates on the machine scale.

    python -m src.paper_numbers   -> results/paper_numbers.json; docs/wp3_table1.md
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from src.audit import wilson
from src.beats import MS
from src.betensky_rule import SOURCES as RULE_SOURCES, source_rows
from src.config import ARTIFACTS, DOCS, N_BOOT, QRS_OFFSET_MS, SEED
from src.criteria_controls import boot_diff
from src.criteria_reliability import kappa, scope
from src.h4 import OUT, demographics as heedb_demographics
from src.mimic_rep import demographics as mimic_demographics
from src.sameday_explore import PRIMARY_BINS, SAME_DAY_BINS, cells, table as interval_table
from src.sex_discordance import compute as sex_numbers
from src.v2_readers import WHO, ecg_table, pair_table, weighted

ROOT = Path(__file__).resolve().parents[1]
OUT3 = ARTIFACTS / "wp3"
SITES = {"mgh": (OUT, "HEEDB-MGH"), "euh": (OUT, "HEEDB-Emory"), "mimic": (OUT3, "MIMIC-IV-ECG")}
PRIMARY = "primary (<= 365 days)"
BIN_KEYS = ("lt1h", "h1to24", "d1to30", "d31to365")
GROUPS = {"signal": ("flat or incomplete lead", "missing lead or samples", "unreadable", "not 500 Hz"),
          "no_sr": ("no SR template", "fewer than 2 usable SR beats", "SR template QRS >= 120 ms")}
H4_ROW = re.compile(r"^\| primary \(<= 365 days\) \| (\d+) \((\d+)\) \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \| ([\de.+-]+) \| "
                    r"([\de.+-]+) \| ([\d.]+) \(([\d.-]+) to ([\d.-]+)\) \|$", re.M)


def prop(x):
    x = np.asarray(x, float)
    k, n = int(x.sum()), len(x)
    if not n:
        return {"k": 0, "n": 0}
    lo, hi = wilson(k, n)
    return {"k": k, "n": n, "f": k / n, "lo": max(float(lo), 0.0), "hi": float(hi)}


def diff(a, b):
    return dict(zip(("m", "lo", "hi"), (float(v) for v in boot_diff(a, b))))


def criteria(d, c):
    out = {}
    for sc, g in (("both", d[d[f"{c}_both"]]), ("first", d[d[f"{c}_first"]]), ("match95", d[d[f"{c}_both"] & (d.r >= 0.95)])):
        p = prop(g[f"{c}_raw"])
        p["kappa"] = float(kappa(g[f"{c}_1"].to_numpy(bool), g[f"{c}_2"].to_numpy(bool)))
        p["ceiling"] = 1 - p["f"] / 2
        out[sc] = p
    return out


def arms(d, c):
    g = d[d[f"{c}_both"]]
    gd = g.dropna(subset=[f"{c}_donor"])
    return {"n": len(g), "n_donor": len(gd), **{k: float(g[f"{c}_{k}"].mean()) for k in ("raw", "shared", "own", "sronly")},
            "donor": float(gd[f"{c}_donor"].mean()), "raw_minus_shared": diff(g[f"{c}_raw"], g[f"{c}_shared"]),
            "own_minus_shared": diff(g[f"{c}_own"], g[f"{c}_shared"]), "own_minus_donor": diff(gd[f"{c}_own"], gd[f"{c}_donor"])}


def standard(g, ref, bins):
    c = cells(g, bins)
    w = np.array([len(x) for x in cells(ref, bins)], float)
    w /= w.sum()
    rng = np.random.default_rng(SEED)
    boots = [sum(wi * x[rng.integers(0, len(x), len(x))].mean() for wi, x in zip(w, c)) for _ in range(N_BOOT)]
    return {"f": float(sum(wi * x.mean() for wi, x in zip(w, c))), "lo": float(np.percentile(boots, 2.5)), "hi": float(np.percentile(boots, 97.5))}


def h4(path):
    m = H4_ROW.search(path.read_text(encoding="utf-8"))
    if m is None:
        raise SystemExit(f"no primary-stratum H4 row in {path}")
    g = m.groups()
    return {"n": int(g[0]), "with_donors": int(g[1]), "d_raw": float(g[2]), "d_own": float(g[3]), "d_null": float(g[4]),
            "p_own_raw": float(g[5]), "p_own_null": float(g[6]), "share": float(g[7]), "lo": float(g[8]), "hi": float(g[9])}


def readers():
    R = OUT3 / "readers"
    e = ecg_table(R)
    pop = json.loads((R / "sample.json").read_text())
    nd, nc = pop["discordant"], pop["concordant"]
    out = {"population": pop, "ecgs": len(e), "pairs": int(e.pair.nunique())}
    for w in WHO:
        t = pair_table(e, w)
        d, c = prop(t[t.stratum == "discordant"].differ), prop(t[t.stratum == "concordant"].differ)
        ok = e[e[f"{w}_ok"]]
        out[w] = {"ecgs_evaluable": len(ok), "pairs_evaluable": len(t), "discordant": d, "concordant": c,
                  "weighted": float(weighted(t, pop)[0]), "weighted_wilson_lo": (nd * d["lo"] + nc * c["lo"]) / (nd + nc),
                  "weighted_wilson_hi": (nd * d["hi"] + nc * c["hi"]) / (nd + nc),
                  "kappa_vs_auto": float(kappa(ok[f"{w}_call"].to_numpy(bool), ok.auto_call.to_numpy(bool)))}
    ok = e[e.reader1_ok & e.reader2_ok]
    out["kappa_reader1_reader2"] = float(kappa(ok.reader1_call.to_numpy(bool), ok.reader2_call.to_numpy(bool)))
    return out


def characteristics(rec, tmpl, pr, sex, age):
    stem = rec.file.str.rsplit("/", n=1).str[-1].to_numpy()
    s, a = sex.reindex(stem[pr.i]).to_numpy(), age.reindex(stem[pr.i]).to_numpy(dtype=float)
    n_pvc = rec.n_pvc.to_numpy()
    rows = [{"age": a[k], "female": None if pd.isna(s[k]) else str(s[k]).upper().startswith("F"), "days": d, "pvc_beats": int(min(n_pvc[i], n_pvc[j])),
             "pvc_qrs": (tmpl[i]["pvc_off"] - tmpl[i]["pvc_on"]) / MS + QRS_OFFSET_MS, "sr_qrs": (tmpl[i]["sr_off"] - tmpl[i]["sr_on"]) / MS + QRS_OFFSET_MS,
             "scope_v2tr": scope(tmpl[i])[1] and scope(tmpl[j])[1]} for k, (i, j, d) in enumerate(zip(pr.i, pr.j, pr.days))]
    return pd.DataFrame(rows)


def summary(d):
    q = lambda x: {"median": float(np.nanmedian(x)), "q1": float(np.nanpercentile(x, 25)), "q3": float(np.nanpercentile(x, 75))}
    known = d.female.notna()
    return {"pairs": len(d), "age": q(d.age), "age_known": float(d.age.notna().mean()), "female": float(d.female[known].astype(bool).mean()), "sex_known": float(known.mean()),
            "days": q(d.days), "same_day": float((d.days < 1).mean()), "pvc_beats": q(d.pvc_beats), "single_beat": float((d.pvc_beats == 1).mean()),
            "pvc_qrs_ms": q(d.pvc_qrs), "sr_qrs_ms": q(d.sr_qrs)}


def main():
    numbers = {k: {} for k in ("table1", "pairs", "flow", "discordance", "arms", "interval", "sameday", "h4", "published_rule")}
    frames, prims, sames = {}, {}, {}
    for key, (root, name) in SITES.items():
        rec = pd.read_parquet(root / f"{key}_records.parquet")
        tmpl = np.load(root / f"{key}_templates.npy", allow_pickle=True)
        pairs = pd.read_parquet(root / f"{key}_pairs.parquet")
        pr = pairs[pairs.stratum == PRIMARY].reset_index(drop=True)
        sex, age = mimic_demographics(rec) if key == "mimic" else heedb_demographics(key, rec.file)
        c = characteristics(rec, tmpl, pr, sex, age)
        frames[key] = c
        numbers["table1"][key] = {"primary": summary(c), "v2tr_scope": summary(c[c.scope_v2tr])}
        why = rec.reason.fillna("unreadable")
        ot = rec[why == "OT-like"].groupby("patient").size()
        numbers["flow"][key] = {"records": len(rec), "patients": int(rec.patient.nunique()),
                                "signal": int(why.isin(GROUPS["signal"]).sum()), "no_sr": int(why.isin(GROUPS["no_sr"]).sum()),
                                "no_pvc": int((why == "no usable PVC").sum()), "multiform": int((why == "more than one PVC cluster").sum()),
                                "not_ot": int((why == "eligible, not OT-like").sum()), "ot_like": int((why == "OT-like").sum()),
                                "patients_2ot": int((ot >= 2).sum())}
        numbers["pairs"][key] = {"records": len(rec), "subjects": int(rec.patient.nunique()),
                                 **{s.split(" ")[0]: int((pairs.stratum == s).sum()) for s in (PRIMARY, "same day", "any interval")}}
        d = pd.read_parquet(root / ("mimic_criteria.parquet" if key == "mimic" else f"{key}_criteria_controls.parquet"))
        numbers["discordance"][key] = {c_: criteria(d, c_) for c_ in ("v2tr", "tzi")}
        numbers["arms"][key] = {c_: arms(d, c_) for c_ in ("v2tr", "tzi")}
        prim = interval_table(rec, tmpl, pr)
        same = interval_table(rec, tmpl, pairs[pairs.stratum == "same day"])
        prims[key], sames[key] = prim, same
        p = prim[prim.scope]
        numbers["interval"][key] = {k: {**prop(g.differ), "share": len(g) / len(p)} for k, g in
                                    ((b, p[p.minutes.between(lo, hi, inclusive="left")]) for b, (_, lo, hi) in zip(BIN_KEYS, PRIMARY_BINS))}
        s = same[same.scope]
        numbers["sameday"][key] = {"in_scope": prop(s.differ), "all_pairs": prop(same.differ), "minutes_median": float(s.minutes.median()), "under_60_min": float((s.minutes < 60).mean())}
    px = {k: v[v.scope] for k, v in prims.items()}
    sx = {k: v[v.scope] for k, v in sames.items()}
    numbers["standardised"] = {"mimic_at_mgh_primary_scope": standard(px["mimic"], px["mgh"], PRIMARY_BINS),
                               "mimic_at_mgh_sameday_scope": standard(sx["mimic"], sx["mgh"], SAME_DAY_BINS)}
    numbers["h4"] = {"mgh": h4(DOCS / "wp2_h4_mgh.md"), "euh": h4(DOCS / "wp2_h4_euh.md"), "mimic": h4(DOCS / "wp3_mimic_replication.md")}
    numbers["readers"] = readers()
    for (name, (root, key, dataset)) in RULE_SOURCES.items():
        d, same_share, lost = source_rows(root, key, dataset)
        numbers["published_rule"][key] = {**prop(d.rule), "template_f": float(d.template.mean()), "difference": diff(d.rule, d.template),
                                          "ecg_changed": float(d.ecg_changed.mean()), "templates_reproduced": float(same_share), "lost": lost}
    numbers["sex"] = sex_numbers()[0]
    (ROOT / "results").mkdir(exist_ok=True)
    clean = lambda o: {k: clean(v) for k, v in o.items()} if isinstance(o, dict) else (round(float(o), 6) if isinstance(o, (float, np.floating)) else (int(o) if isinstance(o, np.integer) else o))
    (ROOT / "results" / "paper_numbers.json").write_text(json.dumps(clean(numbers), indent=1, sort_keys=True) + "\n")

    fmt = lambda q, d=0: f"{q['median']:.{d}f} [{q['q1']:.{d}f}, {q['q3']:.{d}f}]"
    lines = ["# Table 1: pairs by source", "", "Generated by `python -m src.paper_numbers`. Medians [IQR]. Age and sex at the first ECG of the pair.", ""]
    for block, title in (("primary", "Primary-stratum pairs"), ("v2tr_scope", "Pairs with both ECGs in the V2TR scope")):
        t = {k: numbers["table1"][k][block] for k in SITES}
        rows = [("Pairs (one per patient)", lambda v: f"{v['pairs']:,}"), ("Age, years", lambda v: fmt(v["age"]) + (f" (known {v['age_known']:.0%})" if v["age_known"] < 0.995 else "")),
                ("Female", lambda v: f"{v['female']:.1%}" + (f" (known {v['sex_known']:.0%})" if v["sex_known"] < 0.995 else "")),
                ("Interval between ECGs, days", lambda v: fmt(v["days"], 1)), ("Same day (< 24 h)", lambda v: f"{v['same_day']:.1%}"),
                ("Fewest PVC beats per template", lambda v: fmt(v["pvc_beats"])), ("Single-beat PVC template on an ECG", lambda v: f"{v['single_beat']:.1%}"),
                ("PVC QRS duration, ms", lambda v: fmt(v["pvc_qrs_ms"])), ("SR QRS duration, ms", lambda v: fmt(v["sr_qrs_ms"]))]
        lines += [f"## {title}", "", "| | " + " | ".join(n for _, n in SITES.values()) + " |", "|---|---|---|---|"]
        lines += [f"| {label} | " + " | ".join(f(t[k]) for k in SITES) + " |" for label, f in rows] + [""]
    (DOCS / "wp3_table1.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
