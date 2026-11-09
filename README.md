# SR-REF

Code for a study of how often sinus rhythm (SR)-referenced electrocardiographic criteria for outflow
tract premature ventricular complexes (PVCs) change their call between 2 electrocardiograms (ECGs) of
the same patient. The criteria are the V2 transition ratio (V2TR) and the transition zone (TZ) index.
Pairs of ECGs with a matched outflow tract-like PVC come from 2 sites of the Harvard-Emory ECG Database
(HEEDB; Massachusetts General Hospital and Emory) and from MIMIC-IV-ECG. The study also tests whether a
transform estimated from the 2 ECGs' sinus beats corrects the PVC for electrode position.

The MIMIC-IV-ECG replication was preregistered on the Open Science Framework
(https://doi.org/10.17605/OSF.IO/WYX28) and run once. The code and constants of that run are locked by
`prereg/mimic_manifest.json`; `python -m src.mimic_rep run` refuses to start if any locked file differs.
Every definition, change and result is dated in `docs/decisions.md`.

## Data

No data are in this repository: no signals, records, identifiers or row-level results. Paths are set in
`config/data_sources.local.yaml` (copy `config/data_sources.example.yaml`). Row-level outputs go to the
`artifacts` path set there.

| Source | Access | Used for |
|---|---|---|
| HEEDB 5.0 (doi:10.60508/765d-f248) | Brain Data Science Platform, data use agreement | pairs at MGH (development) and Emory |
| MIMIC-IV-ECG 1.0 (doi:10.13026/4nqg-sb35) | PhysioNet, open access | preregistered replication; QRS duration scale; template audit |
| MIMIC-IV 3.1 (doi:10.13026/kpb9-mt58), `hosp/patients` | PhysioNet, credentialed | sex and age for donor matching |
| St Petersburg INCART database (doi:10.13026/C2V88N) | PhysioNet, open access | detector tuning and validation |
| LUDB 1.0.1 | PhysioNet, open access | secondary QRS duration check |

The second PVC cue is a published beat classifier (Köksal, J Electrocardiol 2026;99:154486), used
with its weights and threshold unchanged (`beat_classifier` in the configuration).

## Setup

Python 3.13 with the packages in `requirements.txt`.

## Reproducing the results

Commands run from the repository root. Each writes the report named after the arrow; reports hold
aggregates only.

| Step | Command | Output |
|---|---|---|
| QRS duration scale | `python -m src.qrs_scale` | `docs/wp1_qrs_scale.md` |
| Detector tuning and validation (INCART) | `python -m src.incart tuning`; `python -m src.incart validation --confirm-locked` (once) | `docs/wp1_incart_tuning.md`, `docs/wp1_incart_validation.md` |
| Template audit (blinded readers) | `python -m src.audit sample`, `pages`, `score` | `docs/wp1_audit_sample.md`, `docs/wp1_audit_result.md` |
| Cohorts and matched pairs (HEEDB) | `python -m src.stage_a cohort --site mgh` (about 8 min on 28 processes), `--site euh` | `docs/wp2_cohort_mgh.md`, `docs/wp2_cohort_euh.md` |
| Placement test (H4) | `python -m src.h4 mgh`, `python -m src.h4 euh` | `docs/wp2_h4_mgh.md`, `docs/wp2_h4_euh.md` |
| Reliability of the criteria (HEEDB) | `python -m src.criteria_reliability`, `python -m src.criteria_controls mgh euh` | `docs/wp2_criteria_reliability.md`, `docs/wp2_criteria_controls.md` |
| MIMIC-IV-ECG replication (preregistered) | `python -m src.mimic_rep feasibility`, `dryrun`, `lock`, `run --registered URL` | `docs/wp3_mimic_feasibility.md`, `docs/wp3_mimic_dryrun_euh.md`, `docs/wp3_mimic_replication.md` |
| Time between ECGs (exploratory) | `python -m src.sameday_explore` | `docs/wp3_sameday_exploratory.md` |
| Human measurement of V2 (exploratory) | `python -m src.v2_readers sample`, `pages`, `score` | `docs/wp3_readers_sample.md`, `docs/wp3_readers_result.md` |
| Published sinus measurement (exploratory) | `python -m src.betensky_rule` | `docs/wp3_betensky_rule.md` |
| Discordance by sex (exploratory) | `python -m src.sex_discordance` | `docs/wp3_sex_discordance.md` |
| Table 1 and all reported numbers | `python -m src.paper_numbers` | `docs/wp3_table1.md`, `results/paper_numbers.json` |
| Figures 1-4 | `python -m src.figures` | `results/figures/` |

Table 2 and the supplemental tables are read from `results/paper_numbers.json` and the reports above.
The registered run (`src.mimic_rep run`) refuses a second start.

`src/ningbo.py`, `src/objects.py` and `src/power.py` belong to a separate analysis of ablation-confirmed
origins, which has not been run.

## Licence

MIT (`LICENSE`).
