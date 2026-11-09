# Blinded human measurement of V2 R and S in MIMIC-IV-ECG pairs

Protocol, written 2026-11-09 before any reader has seen the sample. Exploratory. It is outside the
OSF registration 10.17605/OSF.IO/WYX28 and changes none of its results. Code: `src/v2_readers.py`,
`src/v2_readers.html`. Record: `docs/decisions.md` item 37.

## 1. Question

When the automated pipeline finds that the V2 transition ratio (V2TR) call differs between two ECGs
of one patient, does a human reader measuring R and S in V2 find the same? In the pairs where the
pipeline finds no change, does the reader agree?

## 2. Known before this protocol

The automated V2TR discordance in this sample is 36.2% (67 of 185 pairs). The commonest changed call
on HEEDB is a QS-like against an rS PVC in V2. No human has measured V2 R and S in these pairs.

## 3. Sample

- Population: MIMIC-IV-ECG primary-stratum pairs with both ECGs in the published V2TR scope, from the
  registered run (185 pairs).
- Strata from the automated calls: discordant (67 pairs) and concordant (118 pairs).
- Sample: 25 pairs from each stratum, simple random sampling without replacement, seed 1337. The
  author chose this size: 50 ECGs per stratum, 50 pairs and 100 ECGs in total.

## 4. What the readers see

- Each ECG is a separate item, in an order shuffled separately for each reader.
- An item shows lead V2 only: the PVC template and the SR template from 60 ms before to 300 ms after
  the QRS onset, drawn from the template baseline (mean of the 20 ms before the QRS onset).
- These are the same templates that the pipeline measures. When an ECG has one PVC, its template is
  that beat.
- The readers do not see the dataset, the pairing, the stratum, the interval between ECGs, age, sex,
  the automated values, the QRS window or any ratio.

## 5. What the readers do

For the PVC and for the SR of each ECG, each reader:

- marks the R wave, the largest positive deflection of the QRS, by clicking above the baseline;
- marks the S wave, the deepest negative deflection of the QRS, by clicking below the baseline. In a
  QS complex the nadir is marked as S and "no R" is chosen;
- chooses "no R" or "no S" when that wave is absent, or when a deflection is judged to be noise or
  baseline movement and not a true wave;
- chooses "not measurable" when the template cannot be read.

A click snaps to the local extremum within 8 ms. Amplitude is measured from the baseline.
Readers: F. Köksal and F. Levent. They read independently and do not discuss cases before both forms
are returned. There is no consensus step.

## 6. Derived values

- r = R / (R + S), and 0 if both are absent.
- V2TR = r of the PVC / max(r of the SR, 0.01). LVOT call: V2TR ≥ 0.60.
- QS-like PVC: r of the PVC below 0.05.
- These are computed for each reader and for the joint reading, which uses the mean of the two
  readers' R and the mean of their S for each template.
- An ECG is evaluable for a reader when both templates are measured. A pair is evaluable when both
  ECGs are. The joint reading needs both readers. Numbers not evaluable are reported.

## 7. Estimates

No test is confirmatory.

- E1. Share of pairs whose call differs on human measurement, within the automated-discordant and
  within the automated-concordant stratum, for the joint reading and for each reader (Wilson 95% CI).
- E2. Human discordance in the 185 pairs, weighting the two strata by their population sizes:
  (67 × p_discordant + 118 × p_concordant) / 185. Percentile 95% CI from 2000 bootstrap resamples
  within strata, seed 1337. Reported beside the automated 36.2%.
- E3. Reader against automated, per ECG: Cohen's κ of the LVOT call; κ of QS-like PVC; for r of the
  PVC, the mean difference with 95% limits of agreement and the Spearman correlation.
- E4. Reader against reader: the same measures, and κ of the pair-level changed call.
- E5. Among automated-discordant pairs with a QS-like PVC on exactly one ECG: the share where the
  joint reading also finds a QS-like PVC on that ECG and not on the other.

## 8. Reading fixed in advance

Applied to E2, for the joint reading:

- Lower 95% limit above 0.18: human measurement supports the automated finding at the 0.18 threshold.
- Lower limit above 0.10 and at or below 0.18: human measurement supports it at the 0.10 threshold only.
- Upper limit at or below 0.10: human measurement does not support the automated finding, and the
  study report says that the discordance largely reflects the automated measurement.
- Otherwise: inconclusive.

## 9. Order of work

1. This protocol is fixed.
2. The sample and the two reader pages are built. Nobody looks at them before reading starts.
3. The readers return one form each.
4. The score step runs once. If a coding error is found, it is corrected and reported.

## 10. Limits

The readers measure the median templates, so this checks the R and S measurement and not the template
construction. Only V2TR is examined, not the TZ index. There are two readers and one data source.
With 25 pairs per stratum the CI per stratum is wide.
