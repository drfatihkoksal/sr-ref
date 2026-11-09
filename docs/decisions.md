# Decisions

Dated record of approved and frozen definitions. Times are UTC+03:00.

## 2026-11-08: concept

1. **Concept v0.2 proposals approved** (author). Ningbo amplitudes are scaled by one
   dataset-level factor (no patient-level SR scaling, which would pass SR information into M1).
   The power analysis plans on a half effect.
2. **A1: H6 rewritten** (author). RECON-ORIGIN H5 found no AUC gain from concatenating the sinus
   template at any montage. H6 now compares normalisation (ΔLL(M3 − M1)) with low-dimensional
   concatenation (ΔLL(M1+SR − M1), new model M1+SR with an SR PCA basis, k = 3). Exploratory.
3. **A2: power baseline** (author). Effects transfer on the d′ scale (Nakasone Δd′ = 0.368);
   grid of baseline AUC 0.86, 0.90, 0.945; planning scenario baseline 0.945, half effect.

## 2026-11-08: analysis plan draft

4. **Plan items 1–24 confirmed** (author): all defaults of `prereg/ANALYSIS_PLAN_draft.md`
   Section 10, including one-sided α = 0.025 for H1–H3, the new gate K5 (AUC(M1) ≥ 0.75),
   HEEDB-EUH as the confirmatory H4 set, and one INCART re-validation after a revision. Still
   open before the freeze: the audit reader's name (item 3), learned variant specified or dropped
   (item 9), registration (item 24).

## 2026-11-08: WP1 detector

5. **QRS duration scale: constant offset on the machine scale** (author). The slope-envelope QRS
   method keeps its 15% threshold; a constant offset is added to every measured duration
   (`config.QRS_OFFSET_MS` = 20 ms). Reference: MIMIC-IV-ECG machine global QRS, 500 random sinus
   ECGs (bias −20 ms, IQR −26 to −16; MAE after correction 7 ms); LUDB expert QRS is a secondary
   check (`python -m src.qrs_scale` -> `docs/wp1_qrs_scale.md`). The author first chose LUDB
   (+44 ms) and switched after the LUDB convention proved longer than the machine one (median
   111 ms in records without conduction or structural statements), which made 37% instead of 18%
   of the INCART tuning PVCs fall in ineligible windows.
6. **K1 scope: eligible windows** (author). PVC sensitivity and PPV are computed over windows with
   a valid SR template; the share of reference PVCs in ineligible windows is reported beside them.

7. **P wave removed from the PVC rule** (author). On the INCART tuning subset the P criterion
   ("no P wave 120–300 ms before the QRS") did not separate PVCs (PVC sensitivity ≈ 0.73 with it,
   ≈ 0.90 without; PPV unchanged). The PVC rule is prematurity + width + dissimilarity to the SR
   template. The P wave stays in the SR rule ("preceded by a P wave").
8. **K1 fallback** (author). One more rule-refinement round on the tuning subset; if PVC
   sensitivity still cannot reach 0.95, K1 becomes PPV ≥ 0.95 and sensitivity ≥ 0.85, recorded as
   a deviation from the concept. Reason: a false PVC contaminates the PVC template, a missed PVC
   only lowers the number of eligible records.

9. **Candidate PVC rule set** (from the INCART tuning subset, 2026-11-08; locked only at K1). The
   cluster-propagation round did not lift sensitivity above ≈ 0.90 (PPV ≈ 0.94), so the K1
   fallback of item 8 applies (PPV ≥ 0.95 and Se ≥ 0.85, point estimates). Most remaining false
   positives were detections without a reference beat with a very wide measured QRS (median
   188 ms vs 138 ms in true PVCs), so a plausibility bound QRS ≤ 200 ms was added; a correlation
   floor was rejected (true PVCs of opposite polarity have negative correlation). Grid over QRS
   110–120 / 200–220 ms, RR 0.90–0.95, correlation 0.85–0.90 with propagation: a flat optimum
   (Se 0.865–0.88, PPV 0.963–0.968). Chosen, largest margin: detector `nk_vcg`; PVC = preceding
   RR < 0.95 × sinus RR, QRS 115–200 ms, correlation with the SR template < 0.85, plus
   propagation within a PVC cluster. Tuning subset: V Se 0.869 (record-bootstrap 95% CI 0.732 to
   0.943), V PPV 0.967 (0.934 to 0.983); 17.8% of reference PVCs in ineligible windows. The wide
   CI means the validation subset may fall below Se 0.85.

10. **Rules locked** (author, 2026-11-08 20:51 UTC+03:00): item 9 as in `src/config.py` (`DETECTOR`, `RULES`,
    `QRS_OFFSET_MS` = 20). Combined SHA-256 prefix of `src/{config,beats,signal,incart}.py` at the
    lock: `fe45df9385bd724d`. Next: the INCART validation subset, run 1 of at most 2.

11. **K1 moved to the template level** (author, 2026-11-08, after validation run 1 failed K1;
    a deviation from the concept and the draft plan, reported as such).
    - The locked rules (item 10) are kept. The beat-level INCART result is reported as measured:
      V Se 0.788, PPV 0.943, below both K1 thresholds. The second validation run is not used.
    - Reason: every later analysis uses one median PVC template and one median SR template per
      record, not single beats; a missed PVC lowers the yield of eligible records, a false one can
      contaminate a template. The quantity to gate is the template, and it is checked on data the
      detector will run on.
    - **New gate K1-T (template audit)**, fixed before any audit packet is built:
      - Sample: 150 eligible records (pipeline passed: ≥ 1 PVC, ≥ 2 SR, one PVC cluster), 75 from
        HEEDB-MGH and 75 from MIMIC-IV-ECG, drawn with `SEED` among eligible records with an
        OT-like PVC. Records are shuffled; packets show the 12-lead strip with the classified
        beats marked and the PVC and SR templates, without dataset, age or sex.
      - Two independent readers judge each record: PVC template is a ventricular ectopic beat
        (yes/no); SR template is a conducted sinus beat (yes/no). Disagreements are settled by
        consensus; the consensus is the reference.
      - Pass: PVC-template PPV ≥ 0.95 and SR-template PPV ≥ 0.95 (point estimates), each with a
        Wilson 95% CI. Inter-reader agreement is reported (percent agreement and Cohen κ); κ is
        not gated, because with prevalence near 1 it is unstable.
      - Fail: stop and consult before Stage A.
    - The beat-level audit of the draft plan (300 detections) is replaced by K1-T.

12. **K1-T readers and sample** (author, 2026-11-08). Reader 1: Fatih Köksal; reader 2: Fatih
    Levent (plan item 3). The sample was drawn by `python -m src.audit sample` and is fixed
    (`docs/wp1_audit_sample.md`; HEEDB-MGH 1671 and MIMIC 314 candidates processed for 75 + 75
    OT-like eligible records). Reader pages built by `python -m src.audit pages` in
    `ARTIFACTS/wp1/audit/` (each reader in an own random order). The pages contain HEEDB waveforms
    (BDSP DUA): they are read on this machine only, by readers with BDSP HEEDB access. MIMIC-IV-ECG
    is open access (ODbL 1.0, checked on its PhysioNet page 2026-11-08); only MIMIC-IV clinical
    data are credentialed. The author confirmed on 2026-11-08 that reader 2 has the required access.

13. **Lead-quality exclusion implemented** (2026-11-08, raised by the author while reading the audit
    pages). Plan 4.5 excludes records with noisy or missing leads, but `records.process` did not
    check it. Added: a record is excluded when any lead is flat (peak-to-peak < 50 µV) or holds a
    constant run > 0.5 s (lead off, or sequential 2.5-s lead groups with padding). Effects:
    - INCART: no flat leads; 4 validation records (none in tuning) have a constant run. Without
      them the validation result is unchanged (V Se 0.786, PPV 0.943; those 4 records: Se 0.849),
      so incomplete leads do not explain the low sensitivity. The misses concentrate in a few
      records (44% of false negatives in 3 records, 60% in 5). This look at the validation
      results is allowed because the second validation run is not used (item 11).
    - K1-T sample: 6 of 150 cases fail the rule (4 HEEDB-MGH, 2 MIMIC).

14. **Second PVC cue: JECG BeatNet** (author, 2026-11-08; deviation from the concept's purely
    rule-based detector). The MIT-BIH-trained beat classifier of the JECG paper is used unchanged
    (weights and the model card's threshold 0.99; input as in RECON-ORIGIN `signals.py`, leads II
    and V1, 40 Hz low-pass): a beat is a PVC when the locked rules say so or its BeatNet
    probability is ≥ 0.99 (`src/beatnet.py`). INCART tuning subset, eligible windows: rules alone
    Se 0.869 / PPV 0.967; BeatNet alone (≥ 0.99) 0.551 / 0.987; rules OR BeatNet ≥ 0.99
    **0.940 (95% CI 0.837 to 0.984) / 0.964 (0.930 to 0.982)**; with ≥ 0.9, 0.949 / 0.960. The card
    threshold was kept so that no new parameter is tuned on INCART. BeatNet never saw INCART.
15. **Detector v2 locked; the last INCART validation run is used** (author, 2026-11-08 21:18 UTC+03:00). This
    reverses the "not used" part of item 11. Declared: before this run, the run-1 validation
    results were examined once in aggregate for item 13 (false negatives concentrate in a few
    records); that look did not shape item 14, which rests on the tuning subset. Combined SHA-256
    prefix of `src/{config,beats,beatnet,signal,incart}.py` at the lock: `c4bdf6f11c50bb06`. K1-T stays the
    gate (item 11); the beat-level result is reported.

16. **Validation run 2 result** (2026-11-08): V Se 0.927 (0.889 to 0.959),
    V PPV 0.931 (0.858 to 0.975); run 1 was 0.788 / 0.943. Sensitivity
    rose, PPV fell slightly; beat-level K1 not met. No validation run is left.

17. **K1-T sample redrawn with detector v2 and the lead-quality rule** (2026-11-08). The first sample
    (rules only, before item 13) is archived, not deleted, in `ARTIFACTS/wp1/audit_superseded_v1_rules_only/`.
    New sample: `python -m src.audit sample` (HEEDB-MGH 1630 and MIMIC 284 candidates processed for
    75 + 75; `docs/wp1_audit_sample.md`). Reader pages store answers under a key tied to the sample
    (hash of key.csv), so answers given on the first pages do not carry over; reader orders are kept
    in `orders.json` for later rebuilds.
    Label-free yield on the same candidate records, detector v1 -> v2: OT-like eligible HEEDB-MGH
    63 -> 69 of 471, MIMIC 92 -> 109 of 399; "no usable PVC" 92 -> 57 and 116 -> 59; "more than one
    PVC cluster" 24 -> 34 and 27 -> 43. The lead-quality rule removes 36% of 500-Hz HEEDB-MGH
    records with PVC statements (170 of 471; MIMIC 8 of 399), which will shrink the HEEDB cohort of
    Stage A.

18. **K1-T passed** (2026-11-08; `python -m src.audit score` -> `docs/wp1_audit_result.md`). Both
    forms came from the redrawn pages (answers after the 21:22 rebuild) and were complete.
    - PVC template ventricular: 144/150 = 0.960 (Wilson 95% CI 0.915 to 0.982); agreement 98.7%,
      Cohen κ 0.79; 2 disagreements (K077, K128) settled by consensus as "no" (Fatih Köksal and
      Fatih Levent). SR template sinus: 150/150 (0.975 to 1.000).
    - All 6 non-ventricular PVC templates are from HEEDB-MGH (69/75 = 0.92; MIMIC 75/75). Reader
      notes: atrial premature beats (2), aberrant conduction (1), a T wave and a P wave detected as
      beats (1 each), no note (1). The lower CI bound (0.915) is below 0.95; the gate was set on
      the point estimate.
    - Reading took 15 and 12 min (median 4 and 3 s per case).
    - For Stage A: HEEDB carries the H4 test and the frozen objects. The cross-recording PVC match
      (r ≥ 0.85, plan 5.1) should screen out most spurious templates in H4; the template and PCA
      bases fitted on HEEDB-MGH may carry about 8% non-ventricular PVC templates.

19. **κ = 1.0** (HEEDB-MGH, primary stratum, rule of plan 5.2; `docs/wp2_h4_mgh.md`). Median d_own fell
    over the grid {0.01, 0.03, 0.1, 0.3, 1} (0.0256 to 0.0178), so the choice sits at the grid's edge;
    a check beyond the grid, not used for the choice, shows a plateau (κ 3: 0.0178) and a rise (κ 10:
    0.0193; κ → ∞ gives d_raw 0.0212). The grid was not changed.
20. **H4 on HEEDB-MGH (development, descriptive).** Primary stratum, 6104 patients (5807 with null
    donors): median d_raw 0.0212, d_own 0.0178, d_null 0.0319; own < raw p = 3e-53; correctable share
    0.163 (95% CI 0.135 to 0.194). Same day: 0.069 (0.031 to 0.100); any interval: 0.157. The null
    arm raises the discrepancy above raw, so "own < null" is easy; the informative contrast is own <
    raw. The smaller same-day share fits placement varying more between sessions than within a day.
    Ages: AgeAtAcquisition, else birth date (birth date is "None" or "redacted" in part of HEEDB);
    with only birth dates, 18% of patients had no donor.
21. **Stage A objects built on HEEDB-MGH** (`python -m src.objects` -> `docs/wp2_objects.md`;
    `ARTIFACTS/wp2/objects/` with `hashes.json`): 18,210 patients (first eligible OT-like recording of
    the patients processed by the cohort step, i.e. with ≥ 2 PVC records at 500 Hz); noise model
    12,208 matrices (||A − I||_F median 0.176); bases explain 0.867 (raw PVC, k 6), 0.860 (normalised
    PVC, k 6) and 0.731 (SR, k 3). Formal freeze at K2. Still to build: the Ningbo scale factor
    (needs the Ningbo pipeline, label-free).

22. **Ningbo beat classes from the RECON-ORIGIN adjudication** (author, 2026-11-08; deviation from
    "the same detector on all datasets"). The SR-REF detector found a valid PVC and SR pair in only
    176/334 Ningbo records: the single-PVC-cluster rule (79 lost) and the strict SR rule, which
    leaves no SR beat in bigeminy (58 lost), do not suit short EP-lab strips with frequent PVCs.
    Primary: PVC and SR beats as adjudicated in RECON-ORIGIN (blinded to label, sub-site and sex;
    κ 0.97), templates and eligibility by SR-REF code (`records.from_beats`), paced beats and their
    neighbours dropped, AF/flutter-SR records excluded. Sensitivity: the SR-REF detector with the
    dominant PVC cluster and an SR fallback (narrow beats of the SR cluster with a P wave, whatever
    their RR), `python -m src.ningbo --source auto`: 303 eligible. Both reproduce the RECON-ORIGIN
    templates (median r 1.000). Ningbo scale factor (plan 5.6 item 6): 0.0758 µV per raw unit
    (primary route; auto 0.0748), in `objects/ningbo_scale.json`.
23. **K3 not met** (2026-11-08): 310/334 = 92.8% valid (threshold 95%). Excluded: SR template QRS
    ≥ 120 ms 12, fewer than 2 usable SR beats 6, adjudication G0b 3, AF/flutter SR 2, rejected 1.
    As the plan says, the analysis uses the valid subset and the excluded records are reported by
    origin after the freeze. No record was flagged by the template agreement check (all r ≥ 0.90),
    so the adjudication stands as the label-blind visual QC.

24. **HEEDB-EUH cohort** (2026-11-08, label-free; `docs/wp2_cohort_euh.md`). EUH file names carry no
    patient or acquisition time; both come from the metadata (BDSPPatientID, ECGAcquisitionTime);
    sex from its `Sex` column. 21,349 PVC records of 6,915 patients with ≥ 2; 1,329 patients with ≥ 2
    eligible OT-like records; matched pairs: primary 788, same day 311, any interval 943. The
    confirmatory H4 run (`python -m src.h4 euh`) uses κ from MGH, refuses an uncommitted tree and
    records the commit.

25. **H4 supported on HEEDB-EUH** (confirmatory, 2026-11-08, commit `0d6bbed`; `docs/wp2_h4_euh.md`).
    Primary stratum, 788 patients (786 with donors): median d_raw 0.0272, d_own 0.0215, d_null
    0.0373; own < raw p = 2.7e-12, own < null p = 5.2e-84; correctable share 0.210 (95% CI 0.131 to
    0.283). Same day: 0.088 (−0.071 to 0.214); any interval: 0.209. Replicates MGH (0.163).
26. **K2 met** (2026-11-08): K1-T passed (item 18), H4 done (item 25), code committed (`0d6bbed`),
    Stage A objects frozen with these SHA-256 content hashes (`ARTIFACTS/wp2/objects/hashes.json`):
      - `kappa.json`: `6143e79214d30d4bd68b817cf3d0144520a0897f56e0bdb4ea6024a2785ccfa4`
      - `noise_A.npy`: `32888791815c48c27ee9b4080d83dd346744a0b3ac3aa9efadfa2ef52bdb2b1b`
      - `sr_template.npy`: `d65572817099a9c832e8132b230cb89a9cad00c125ca703a88d8b9b4105f93f4`
      - `sr_loop_rms.json`: `020c10db95d77254e0e3ee78e48ceba7cf09b3364e3afb6abe8c3d4319c9e569`
      - `basis_raw_pvc.npz`: `51a2b67bc5287ee1953a8c7cc1724b4f65a1cfd2e716dcff25a7d84522109036`
      - `basis_norm_pvc.npz`: `1a7081313cfd0e89e40cd9e83f6ba15aa947dabce5c3cef91c6e0a78e2a7eff3`
      - `basis_sr.npz`: `0375289e9c6cdc0f002c1442744ab4adce0c6e5e21afbeb460810c70f6f078c2`
      - `ningbo_scale.json`: `b665d2036523dc82fcb62c89b72fe84d4a66be161a7089046a2c9a11787b3c0d`
    Any change to them from here on needs the author's approval and a dated note here.

27. **H1–H3 decision rule: Nadeau–Bengio; power** (author, 2026-11-08; `python -m src.power` ->
    `docs/wp2_power.md`). The power simulation (n 310, 239/71, plan's nested CV, 500 sets per cell)
    showed that the draft rule, "upper 95% bootstrap CI bound of ΔLL < 0", has a type I error of
    0.04–0.06 (nominal 0.025): the patient bootstrap does not refit the models. The label permutation of plan 6.6 was first
    proposed as the new primary test and accepted; a check before implementing showed it invalid for H1 (type I 0.345 in 200 sets: permuting the
    label makes both models uninformative, so it tests "no association", not "equal
    information"). The author then chose the Nadeau–Bengio corrected resampled t-test (100
    per-fold ΔLL values, df 99, one-sided p < 0.025): type I 0.02–0.03 in all no-effect cells. The
    bootstrap CI stays as the estimate; the label permutation is kept for "M3 against chance"
    only. Power with Nadeau–Bengio, planning scenario (baseline 0.945, half effect): 0.16, 0.29,
    0.42 at ρ 0.6, 0.8, 0.9 -> recorded as "H1 may be inconclusive" (plan 6.11); full effect:
    0.50, 0.74, 0.93.

28. **H4 in clinical units** (post hoc, label-free, defined after the EUH result; `python -m
    src.h4_clinical` -> `docs/wp2_h4_clinical.md`). Primary-stratum pairs, the own-SR VCG transform
    applied to the dipolar part of the 8-lead PVC (L' = K⁺AKL + (I − K⁺K)L). Findings:
    - The VCG correction does not stabilise lead-level PVC criteria: graded precordial transition
      differs between the two ECGs in 28.2% raw vs 28.4% corrected (MGH; EUH 29.6% vs 29.3%);
      median |change| 0.06 vs 0.07 leads. Only large shifts fall: transition change ≥ 1 lead 6.7% ->
      5.5% (MGH), 9.0% -> 6.2% (EUH).
    - Published SR-referenced criteria, each ECG with its own SR, give a different RVOT/LVOT call
      between two ECGs of the same patient with a matched PVC morphology (r ≥ 0.85, ≤ 365 days):
      V2 transition ratio 12.8% (MGH) and 14.2% (EUH); TZ index 11.4% and 10.9%.
    - A first version counted a corrected change when either direction changed, which works against
      the correction; the reported version averages both directions, like the continuous measures.
    Reading: the 3 x 3 VCG transform captures the global orientation, not chest-electrode-specific
    displacement, which drives the precordial transition.

29. **8 x 8 lead-space transform** (plan 5.2 sensitivity estimator, run after item 28; descriptive;
    `python -m src.h4_leads` -> `docs/wp2_h4_leads.md`). κ8 = 1.0 by the plan's rule on MGH (edge of
    the grid again; check κ8 = 3: 0.0288 vs 0.0293, κ8 = 10: 0.0311). Correctable share in lead space
    0.232 (MGH) and 0.297 (EUH), above the VCG version (0.163, 0.210). Clinical units, raw -> 8 x 8
    corrected: transition change ≥ 1 lead 6.7% -> 4.4% (MGH) and 9.0% -> 5.8% (EUH); graded transition
    differs 28.2% -> 25.7% (p 2e-10) and 29.6% -> 27.5% (p 0.06); median |change| 0.06 -> 0.07-0.08
    leads (small shifts grow slightly, large shifts shrink). The lead-space transform removes about a
    third of the large transition shifts; most between-ECG change remains.

30. **Test-retest reliability of the published criteria** (post hoc, label-free, prompted by the
    author's question whether reliability limits the claimed accuracy; `python -m
    src.criteria_reliability` -> `docs/wp2_criteria_reliability.md`). Primary-stratum pairs, each ECG
    on its own templates. V2 transition ratio (≥ 0.60) within its published scope (LBBB, inferior
    axis, V3 transition; scope from the first ECG): call differs in 23.2% (MGH, 2088 pairs) and
    23.3% (EUH, 270), κ 0.54 and 0.53; same day 14.9% (both sites), 1-365 days 30.6% and 27.8%; with
    a stricter match (r ≥ 0.95) 15.8% and 16.4%. Implied single-ECG accuracy ceiling when the origin
    is fixed: ≤ 88.4% (in scope), ≤ 84.7% (different days). Published claims: sensitivity 95%,
    specificity 100% (retrospective, n 40) and 91% accuracy (prospective, n 21) (Betensky 2011, full
    text); pooled specificity below 0.90 in the 2025 network meta-analysis (He et al., Heart). TZ index
    (< 0) in scope (LBBB, inferior axis): 15.6% / 15.5%, κ 0.68; same day 10.6% / 11.6%; 1-365 days
    20.5% / 17.4%. Caveats: automated template measurement (Betensky took the largest sinus R and S
    over 10 s), unselected HEEDB PVC patients rather than ablation candidates, and a morphology
    match that may admit different foci.

31. **Multifocality and the V2TR test-retest result** (post hoc, label-free, prompted by the author's
    question; `python -m src.criteria_reliability multifocal` -> `docs/wp2_criteria_multifocal.md`).
    V2TR pairs in scope. Multifocality does not explain the call changes: patients without multiform
    PVCs in any ECG 22.3% (with 25.6%, MGH); patients with ≥ 3 ECGs whose PVCs all match 24.3% (MGH;
    EUH 25.0%). Mapping one ECG's PVC into the other's lead space with the 8 x 8 own-SR transform
    lowers the rate from 23.2% to 14.0% (MGH) and 23.3% to 13.1% (EUH), which a different focus would
    not allow. What changes is the PVC numerator, not the SR denominator (SR-driven in 11%): in the
    differing pairs the PVC is QS or near-QS in V2 (r < 0.05) in one ECG in 59% (MGH) and 64% (EUH),
    with r_V2 median 0.03 vs 0.27 (absolute change median 0.19 and 0.27; agreeing pairs 0.03 and
    0.05). A first guess (an unstable SR denominator) was wrong.
    **Corrected 2026-11-08 (item 32):** the figures 14.0% and 13.1% came from a counting error. The
    correct values are 20.2% and 18.1%, and the clause "which a different focus would not allow" is
    withdrawn. The checks on multiform PVCs and on patients with ≥ 3 ECGs, and the reading that the
    PVC numerator drives the V2TR change, stand.

32. **Controls of the test-retest result; counting error corrected** (post hoc, label-free, 2026-11-08;
    `python -m src.criteria_controls` -> `docs/wp2_criteria_controls.md`).
    - *Error.* The two-direction mean of the V2TR arms added NumPy booleans, which add as logical OR, so
      a pair flipping in both directions counted 0.5 instead of 1. It was found because the raw rate
      exceeded the sum of the two single-factor arms, which is impossible. It affected
      `criteria_reliability.multifocal` (item 31) and the first run of `criteria_controls`. The TZ index
      arms and the graded-transition counts of `h4_clinical` and `h4_leads` (items 28, 29) use Python
      bools and are unaffected. Raw rates (item 30) are unaffected.
    - *Scope on both ECGs.* V2TR call differs in 22.3% (MGH, 1631 pairs) and 24.2% (EUH, 211), against
      23.2% and 23.3% with scope from the first ECG; ceilings ≤ 88.9% and ≤ 87.9%. TZ index 15.3% and
      14.9%; ceilings ≤ 92.3% and ≤ 92.6%. Selection on the first ECG does not explain the result.
    - *V2TR arms, both ECGs in scope (MGH; EUH).* Raw 22.3; 24.2. Shared SR 21.7; 25.6 (raw minus shared
      +0.6, 95% CI −0.4 to +1.5; −1.4, −3.8 to +1.2). SR only 5.9; 6.2. Own 8 x 8 transform 19.7; 18.5.
      Donor transform 24.7; 28.2. Own minus shared SR −2.0 points (−3.2 to −0.6) and −7.1 (−11.1 to −3.1).
      Own minus donor −4.9 (−6.1 to −3.7) and −9.7 (−13.1 to −6.2).
    - *TZ index arms, both ECGs in scope.* Raw 15.3; 14.9. Shared SR 7.9; 8.4 (raw minus shared +7.4,
      +6.3 to +8.5; +6.4, +3.4 to +9.5). SR only 13.2; 14.0. Own transform 7.2; 7.6 (own minus shared
      −0.7, −1.4 to −0.1; −0.8, −2.6 to +1.1).
    Reading: the V2TR change comes from the PVC, and sharing the SR does not lower it. The TZ index
    change comes largely from the SR, because the graded SR transition moves between ECGs. The own-SR
    8 x 8 transform removes a modest part of the V2TR discordance (12% MGH, 24% EUH, against raw, both
    ECGs in scope), and most of it remains. A foreign transform makes the call less stable.

33. **MIMIC-IV-ECG replication: draft registration** (2026-11-08; form text kept locally; not approved,
    not frozen, not registered). Metadata-only feasibility (`python -m src.mimic_rep feasibility` ->
    `docs/wp3_mimic_feasibility.md`): 5128 subjects with ≥ 2 PVC-statement records and no prior waveform
    contact; projected 541-611 primary pairs and 145-164 V2TR pairs with both ECGs in scope. Prior contact,
    all excluded in the draft: QRS scale (497 subjects), K1-T and yield check (first 400 candidates of one
    fixed order). Proposed sequence: R1 discordance > 0.10 (power 0.94-0.99), R2 H4 replication, R3
    discordance > 0.18 (power 0.23-0.26 at f 0.22), R4 own transform below shared SR. Waveforms are not
    read and no pair is built until the author approves the text and the code is locked.

34. **MIMIC-IV-ECG replication approved and locked** (author, 2026-11-09). The author approved the
    registration text with the primary threshold 0.10, exclusion of all subjects with prior waveform
    contact, and registration on OSF. `python -m src.mimic_rep run --registered URL` reads the waveforms of
    the sample once. It refuses code or constants that differ from `prereg/mimic_manifest.json` (SHA-256
    `56365530a6408135c674c816042b4a71414678392cc4e7f3d608d1ec004a1611`), an uncommitted working tree and a second start. Checks before the lock: `dryrun` on HEEDB-EUH
    re-derived the stored pairs and reproduced the known results (788 primary pairs, 211 V2TR pairs, 24.2%;
    H4 p 2.7e-12 and 5.2e-84, share 0.210; `docs/wp3_mimic_dryrun_euh.md`); `smoke` ran the record pipeline
    on 30 already-read records of excluded subjects (outcome counts and sex and age coverage only). Same-day
    calls in the replication come from the same-day stratum pairs, as in the registration text; the
    development tables of items 30-32 took same-day pairs from the primary stratum.

35. **MIMIC-IV-ECG replication: registered and run once** (2026-11-09; OSF DOI 10.17605/OSF.IO/WYX28,
    https://osf.io/wyx28; commit `5783bfe`; `docs/wp3_mimic_replication.md`). 17,195 records of 5,128
    subjects; 1,169 subjects with at least two eligible OT-like records; 720 primary pairs; 185 V2TR pairs
    with both ECGs in scope. All four hypotheses met: R1 and R3, V2TR discordance 36.2% (Wilson 95% CI 29.6%
    to 43.4%; ceiling ≤ 81.9%); R2, H4 own < raw p 1.5e-11, own < null p 1.1e-69, correctable share 0.121
    (0.070 to 0.194); R4, own transform minus shared SR −4.6 points (−8.9 to −0.3). The estimate is above the
    registered prediction (0.22 to 0.24), and the same-day stratum (38.8%, 67 pairs) does not fall as on
    HEEDB-EUH (14.4%, 97 pairs); the cause is not known. Descriptive: TZ index 19.8% (16.3% to 23.9%); raw
    minus shared SR +0.8 points for V2TR and +10.7 for the TZ index, as on HEEDB; the own transform adds
    nothing for the TZ index (+0.1 points). The run may not be repeated; anything further on this sample is
    exploratory and is labelled so.

36. **Same-day V2TR discordance between sources** (exploratory, post hoc, label-free, outside the
    registration; 2026-11-09; `python -m src.sameday_explore` -> `docs/wp3_sameday_exploratory.md`). Prompted by
    item 35. The registered results do not change. V2TR, both ECGs in scope.
    - *Interval mix.* Same-day pairs lie a median of 356 minutes apart in MIMIC-IV-ECG, against 23 (MGH) and 21
      (EUH); under 60 minutes apart in 28%, against 54% and 57%. At the MGH interval mix the MIMIC same-day
      rate falls from 38.8% to 27.4% (95% CI 18% to 39%; MGH 14.3%) and the primary-stratum rate from 36.2% to
      31.2% (24% to 38%; MGH 22.3%). The interval mix explains about half of the same-day gap to MGH and
      about a third of the primary-stratum gap.
    - *By interval, primary stratum (MGH; EUH; MIMIC).* Under 1 h 8.7%; 14.0%; 12.5% (16 MIMIC pairs). 1-24 h
      20.9%; 16.1%; 47.4%. 1-30 days 32.5%; 36.7%; 45.2%. 31-365 days 30.6%; 24.3%; 27.5%.
    - *Remaining difference, not separable.* In same-day pairs MIMIC has less alike SR templates (median r
      0.974, against 0.991 and 0.980), about twice the V2 template noise (5.1 µV, against 2.9 and 2.3), more
      pairs with a QS-like V2 PVC on one ECG only (37%, against 13% and 23%) and more values near the threshold
      (40%, against 26% and 27%). Noise does not predict discordance within HEEDB, and almost all MIMIC pairs
      lie above the pooled noise median, so noise cannot be tested within MIMIC. MIMIC same-day pairs, all scopes:
      same cart 10.8% (102 pairs), different cart 21.4% (173); at ≥ 60 minutes apart 14.3% (49) and 22.6% (164),
      with overlapping CIs. HEEDB metadata have no device field.
    Reading: report the discordance by interval, not as one number. The causes of the MIMIC excess at 1 h to
    30 days were not established.

37. **Human measurement of V2 R and S** (exploratory, outside the registration; author, 2026-11-09; protocol
    `prereg/V2_READER_VALIDATION.md`). Purpose: to check whether the automated V2TR discordance holds on human
    measurement. The author chose 100 ECGs: 25 automated-discordant and 25 automated-concordant pairs of the 185
    MIMIC-IV-ECG primary-stratum pairs with both ECGs in the V2TR scope. Readers F. Köksal and F. Levent, each ECG a
    separate item in a separately shuffled order, blinded to pairing, stratum, interval, age, sex and automated values.
    The page and the score step were tested on synthetic data only. No real page was opened during development.
    SHA-256: protocol `f627f2befdda5481d0430f2355da67aa9c9430abc28b7083e158b81d87bb251c`; code `a686ec0040c2e888591d221d1128f8674a709ed5ac7daa3cbcfd185a79b6d951`; page
    `2b065bb8ee685d4caff18bbd00c76423ce8f8c7f0b478ada46d0b042b521f9c5`; sample key `64c2bdb6176c0f0e8e9c6083e86e0f3127c985833219524a5c7756e22d455ef4`.

38. **Human measurement of V2 R and S: result** (exploratory, outside the registration; 2026-11-09; item 37;
    `docs/wp3_readers_result.md`). Joint reading: the call differed in 23 of 23 evaluable automated-discordant pairs and in
    0 of 25 automated-concordant pairs; κ of the LVOT call against automated 0.96; reader against reader 1.00. Weighted
    human discordance 36.2%. The protocol's bootstrap CI has zero width because both stratum rates lie at a boundary;
    combined Wilson limits 31.0% to 44.7%. Reading of protocol section 8: supported at the 0.18 threshold. On the same
    templates, readers marked the same deflections as the pipeline, so the changed calls are not an artefact of the
    automated R and S extraction. Part of the agreement is structural (click snapping). Template construction and caliper
    measurement on raw beats were not tested. The two forms carry consecutive timestamps (21:45 to 21:57 and 21:58 to
    22:05 UTC; median 6 and 4 s per ECG); the author confirmed on 2026-11-09 that F. Levent completed the second form
    without seeing the first reader's answers.

39. **Published sinus measurement rule** (exploratory, post hoc, label-free, outside the registration; 2026-11-09;
    `python -m src.betensky_rule` -> `docs/wp3_betensky_rule.md`). V2TR with the sinus r from the largest R and the
    largest S over 10 s, as in Betensky 2011, in place of the median SR template; the PVC template is unchanged. Primary
    stratum, both ECGs in scope. Discordance 22.3% to 22.4% (MGH; difference +0.2 points, 95% CI −0.7 to +1.1), 24.2% to
    26.1% (Emory; +1.9, −0.5 to +4.7), 36.2% to 35.7% (MIMIC; −0.5, −2.7 to +1.6). Records were processed again and the
    stored SR templates were reproduced in 100%. The discordance does not depend on this choice of sinus measurement.
    The PVC measurement was not varied.

40. **Discordance by sex** (post hoc, exploratory, label-free, outside the registration; 2026-11-09; `python -m
    src.sex_discordance` -> `docs/wp3_sex_discordance.md`). Sex as recorded in the source, at the first ECG of the pair.
    V2TR, women against men: MGH 22.5% and 22.0% (difference +0.4 points, 95% CI −3.5 to +4.3); Emory 20.2% and 29.3% (−9.2,
    −20.7 to +2.7); MIMIC-IV-ECG 33.3% and 38.5% (−5.1, −19.3 to +8.7). TZ index: 16.4% and 14.5%, 16.0% and 13.9%, 18.8% and
    20.6%. Pooled, adjusted for source and interval class: V2TR −1.3 points (−4.8 to +2.4), TZ index +1.5 (−0.7 to +3.6). No
    difference by sex was detected; the CIs at Emory and in MIMIC are wide. H4 was not stratified by sex.

41. **Public code repository** (author, 2026-11-09). The public branch `main` starts from a single commit of the analysis
    code, reports and aggregate results. The development history stays in the local branch `wp1-wp2-stage-a`; commit
    identifiers cited above (for example `0d6bbed`, the Emory run, and `5783bfe`, the locked MIMIC-IV-ECG run) refer to
    it. The locked MIMIC-IV-ECG code is unchanged and can be checked against `prereg/mimic_manifest.json`. The concept
    note and the unfrozen analysis plan draft are not public. Figures 1-4 are drawn by `src/figures.py`. Release v1.0.0 is
    archived at Zenodo (version DOI 10.5281/zenodo.24085929; concept DOI 10.5281/zenodo.24085928).

## Open

- Before the freeze (K4): Ningbo outer folds (6.5), plan items 9 (learned variant) and 24
  (registration), plan v1.0 with the deviations folded in (record the power result in Section 0).
- Plan items 9 (learned variant) and 24 (registration).
