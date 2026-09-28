# Messy Lab · Machine answers

Computed by `questions/answers.py` against `db/messy_lab.duckdb`. Each answer sits next to the query that produced it. These still need a hand check against the raw files before they go into the answer key.

Shared helper used by Q18–Q20 and Q24–Q25:

```sql
CREATE OR REPLACE TEMP MACRO treatment_effect(sym) AS TABLE
WITH expr AS (
    SELECT s.*, m.value
    FROM sample s JOIN measurement m USING (sample_id) JOIN gene g USING (hgnc_id)
    WHERE g.symbol = sym AND s.include_in_analysis),
sig AS (   -- what each sample received, leaving out the solvent
    SELECT sample_id, string_agg(compound_id, '+' ORDER BY compound_id) AS compounds,
           max(exposure_h) AS exposure_h, max(exposure_type) AS exposure_type
    FROM sample_treatment WHERE compound_id NOT IN ('ETHANOL', 'DMSO') GROUP BY 1),
tre AS (
    SELECT e.study_id, e.batch, e.cell_line_id, e.derived_line_id, e.genetic_change, e.genetic_change_target,
           sig.compounds, sig.exposure_h, sig.exposure_type, avg(e.value) AS t, count(*) AS n_treated,
           list(e.sample_id ORDER BY e.sample_id) AS treated_ids
    FROM expr e JOIN sig USING (sample_id) WHERE e.role = 'treated' GROUP BY ALL),
ctl AS (   -- controls: solvent only or nothing added, same run, same cells, same genetic change
    SELECT study_id, batch, cell_line_id, derived_line_id, genetic_change, genetic_change_target,
           avg(value) AS c, count(*) AS n_control, list(sample_id ORDER BY sample_id) AS control_ids
    FROM expr WHERE role IN ('vehicle_control', 'untreated_control') GROUP BY ALL)
SELECT tre.*, ctl.n_control, ctl.control_ids, tre.t - ctl.c AS log2_change
FROM tre LEFT JOIN ctl
  ON  ctl.study_id = tre.study_id
  AND ctl.batch IS NOT DISTINCT FROM tre.batch
  AND ctl.cell_line_id = tre.cell_line_id
  AND ctl.derived_line_id IS NOT DISTINCT FROM tre.derived_line_id
  AND ctl.genetic_change IS NOT DISTINCT FROM tre.genetic_change
  AND ctl.genetic_change_target IS NOT DISTINCT FROM tre.genetic_change_target;
```

## Q01

*(CEO)* How many MCF-7 samples were exposed to 4-OHT, across all five studies? Count resistant MCF-7 lines grown in 4-OHT. Don't count samples whose tamoxifen form isn't stated.

**Answer:** `{"total": 13, "per_study": {"GSE111151": 1, "GSE117942": 2, "GSE26298": 4, "GSE4025": 6}}`

```sql
SELECT st.geo_accession, count(DISTINCT s.sample_id) AS samples
FROM sample_treatment t JOIN sample s USING (sample_id) JOIN study st USING (study_id)
WHERE t.compound_id = 'OHT' AND s.cell_line_id = 'CVCL_0031'
GROUP BY ALL ORDER BY 1
```

## Q02

*(Head of data)* How many control samples does each study have? A control is a sample given no drug, hormone or growth factor (only the solvent, or nothing at all), including 0-hour baselines and parental lines, whatever was done to its genes.

**Answer:** `{"GSE111151": 4, "GSE117942": 14, "GSE21618": 31, "GSE26298": 8, "GSE4025": 5}`

```sql
SELECT st.geo_accession, count(*) AS controls
FROM sample s JOIN study st USING (study_id)
WHERE s.role IN ('vehicle_control', 'untreated_control', 'parental')
GROUP BY ALL ORDER BY 1
```

## Q03

*(Head of data)* Which cell lines are in the data, and how many samples does each have? Count tamoxifen-resistant derivatives under the line they were made from.

**Answer:** `{"MCF-7": 188, "CAMA-1": 14, "HCC1500": 14, "BT-474": 13, "T-47D": 13, "EFM-19": 10, "MDA-MB-330": 10, "ZR-75-1": 3}`

```sql
SELECT cl.preferred_name AS cell_line, count(*) AS samples
FROM sample s JOIN cell_line cl ON cl.cell_line_id = s.cell_line_id
GROUP BY ALL ORDER BY samples DESC, cell_line
```

## Q04

*(Scientist)* Which compounds were added to cells, and in how many samples each? Count the solvents too. Keep tamoxifen of an unstated form separate from 4-OHT.

**Answer:** `{"17beta-estradiol": 76, "tamoxifen (form not stated)": 49, "heregulin-beta1 (176-246)": 42, "4-hydroxytamoxifen": 31, "GDC-0810 (brilanestrant)": 14, "dimethyl sulfoxide": 14, "fulvestrant": 14, "ethanol": 13, "GDC-0927": 6, "GNE-274": 6}`

```sql
SELECT c.preferred_name AS compound, count(DISTINCT t.sample_id) AS samples
FROM sample_treatment t JOIN compound c USING (compound_id)
GROUP BY ALL ORDER BY samples DESC, compound
```

## Q05

*(CEO)* How many distinct tamoxifen-resistant cell lines do we have, and which lab made each?

**Answer:** `{"total": 13, "per_study": {"GSE111151": ["BT-474 Tam1 (Helsinki)", "BT-474 Tam2 (Helsinki)", "MCF-7 Tam1 (Helsinki)", "T-47D Tam1 (Helsinki)", "T-47D Tam2 (Helsinki)", "ZR-75-1 Tam1 (Helsinki)", "ZR-75-1 Tam2 (Helsinki)"], "GSE21618": ["MCF-7 TamR clone #1 (RIKEN)", "MCF-7 TamR clone #2 (RIKEN)", "MCF-7 TamR clone #3 (RIKEN)", "MCF-7 TamR clone #4 (RIKEN)", "MCF-7 TamR clone #5 (RIKEN)", "MCF-7 TamR clone #6 (RIKEN)"]}}`

```sql
-- resistant lines actually measured; RIKEN's time-course clone is one of its six numbered clones
-- (cell_line.notes), so it is not counted as a seventh line
SELECT st.geo_accession, cl.cell_line_id, cl.preferred_name
FROM cell_line cl JOIN study st ON st.study_id = cl.derived_by_study
WHERE cl.is_derived
  AND cl.cell_line_id IN (SELECT derived_line_id FROM sample)
  AND cl.cell_line_id <> 'S2_TAMR_SELECTED'
ORDER BY 1, 2
```

## Q06

*(Head of data)* How many genes can we compare across all five studies?

**Answer:** `11781`

```sql
SELECT count(*) AS genes
FROM (SELECT hgnc_id FROM feature_map WHERE chosen_for_gene
      GROUP BY 1 HAVING count(DISTINCT study_id) = 5)
```

## Q07

*(Scientist)* Which samples are the right controls for the RAR-alpha knockdown samples in GSE26298?

**Answer:** `["GSM645718", "GSM645719"]`

```sql
SELECT c.sample_id
FROM sample k JOIN sample c ON c.study_id = k.study_id AND c.batch = k.batch
WHERE k.study_id = 'S3' AND k.genetic_change = 'knockdown'
  AND k.genetic_change_target = (SELECT hgnc_id FROM gene WHERE symbol = 'RARA')
  AND c.genetic_change = 'non_targeting_control' AND c.role = 'vehicle_control'
GROUP BY 1 ORDER BY 1
```

## Q08

*(Scientist)* How many samples in GSE21618 never had anything added to them before their RNA was taken?

**Answer:** `31`

```sql
SELECT count(*) AS never_treated
FROM sample s
WHERE s.study_id = 'S2'
  AND NOT EXISTS (SELECT 1 FROM sample_treatment t WHERE t.sample_id = s.sample_id)
```

## Q09

*(Scientist)* What concentration of 4-OHT did each study use?

**Answer:** `{"GSE111151": "1000", "GSE117942": "1000", "GSE21618": "not 4-OHT: tamoxifen of unstated form, dose not stated", "GSE26298": "100, 500", "GSE4025": "10"}`

```sql
SELECT st.geo_accession,
       CASE WHEN count(t.sample_id) FILTER (WHERE t.compound_id = 'OHT') > 0
            THEN string_agg(DISTINCT CAST(CAST(t.dose_nM AS INTEGER) AS VARCHAR), ', ' ORDER BY CAST(CAST(t.dose_nM AS INTEGER) AS VARCHAR))
                 FILTER (WHERE t.compound_id = 'OHT')
            WHEN count(t.sample_id) FILTER (WHERE t.compound_id = 'TAM_UNSPECIFIED') > 0
            THEN 'not 4-OHT: tamoxifen of unstated form, dose not stated'
            ELSE 'no 4-OHT' END AS oht_dose_nM
FROM study st LEFT JOIN sample s USING (study_id) LEFT JOIN sample_treatment t USING (sample_id)
GROUP BY 1 ORDER BY 1
```

## Q10

*(Scientist)* How long were cells exposed to 4-OHT before RNA was collected, in each study that used 4-OHT?

**Answer:** `{"GSE111151": "chronic (months)", "GSE117942": "24 h", "GSE26298": "not stated", "GSE4025": "24 h"}`

```sql
SELECT st.geo_accession,
       CASE WHEN bool_and(t.exposure_type = 'chronic') THEN 'chronic (months)'
            WHEN count(t.exposure_h) = 0 THEN 'not stated'
            ELSE string_agg(DISTINCT CAST(CAST(t.exposure_h AS INTEGER) AS VARCHAR), ', ') || ' h' END AS exposure,
       any_value(t.exposure_note) AS note
FROM sample_treatment t JOIN sample s USING (sample_id) JOIN study st USING (study_id)
WHERE t.compound_id = 'OHT'
GROUP BY 1 ORDER BY 1
```

## Q11

*(Head of data)* Which treated samples have no control (solvent-only or untreated) with the same cells in the same experimental run?

**Answer:** `["GSM3315688", "GSM3315689", "GSM3315690", "GSM3315691", "GSM3315692", "GSM3315693", "GSM3315694", "GSM3315695"]`

```sql
SELECT t.sample_id
FROM sample t
WHERE t.role = 'treated' AND NOT EXISTS (
    SELECT 1 FROM sample c
    WHERE c.study_id = t.study_id
      AND c.batch IS NOT DISTINCT FROM t.batch
      AND c.cell_line_id = t.cell_line_id
      AND (c.derived_line_id IS NOT DISTINCT FROM t.derived_line_id OR c.role = 'parental')
      AND c.genetic_change IS NOT DISTINCT FROM t.genetic_change
      AND c.genetic_change_target IS NOT DISTINCT FROM t.genetic_change_target
      AND c.role IN ('vehicle_control', 'untreated_control', 'parental'))
ORDER BY 1
```

## Q12

*(Scientist)* Which studies removed estrogen from the growth medium before treating the cells, and which didn't?

**Answer:** `{"deprived": ["GSE117942", "GSE21618", "GSE26298", "GSE4025"], "not_deprived": ["GSE111151"]}`

```sql
SELECT st.geo_accession, string_agg(DISTINCT s.hormone_deprived, ',') AS hormone_deprived
FROM sample s JOIN study st USING (study_id) GROUP BY 1 ORDER BY 1
```

## Q13

*(Head of data)* Which samples have labels that contradict each other?

**Answer:** `["GSM539725", "GSM539726", "GSM539727", "GSM539728"]`

```sql
SELECT DISTINCT entity_id AS sample_id FROM issue
WHERE check_name = 'fields_contradict' AND entity_type = 'sample' ORDER BY 1
```

## Q14

*(Head of data)* Is any dose written at an impossible level anywhere in the data? If so, where, and what should it be?

**Answer:** `[{"study": "GSE117942", "as_written": "1 M", "corrected_nM": 1000}]`

```sql
SELECT st.geo_accession, i.observed AS as_written, i.message
FROM issue i JOIN study st USING (study_id)
WHERE i.check_name = 'dose_above_100uM'
```

## Q15

*(Head of data)* In which studies are the deposited expression values not what GEO's own processing description says they are?

**Answer:** `["GSE117942", "GSE4025"]`

```sql
-- mismatches against GEO's own processing text (the S3 flag compares GEO with the paper instead)
SELECT DISTINCT st.geo_accession FROM issue i JOIN study st USING (study_id)
WHERE i.check_name = 'value_scale_mismatch' AND i.expected NOT ILIKE '%paper%' ORDER BY 1
```

## Q16

*(Head of data)* How many gene names have been corrupted into dates by Excel, and in which study?

**Answer:** `{"count": 28, "study": ["GSE117942"]}`

```sql
SELECT st.geo_accession, count(*) AS names FROM issue i JOIN study st USING (study_id)
WHERE i.check_name = 'date_like_gene_name' GROUP BY 1
```

## Q17

*(Scientist)* Which samples carry a replicate label that is repeated or wrong?

**Answer:** `["GSM645720", "GSM645721"]`

```sql
SELECT DISTINCT entity_id AS sample_id FROM issue
WHERE check_name = 'duplicate_titles' AND message LIKE '%replicate number is repeated%' ORDER BY 1
```

## Q18

*(Scientist)* Is GREB1 induced by 24 hours of 4-OHT in ordinary MCF-7, compared with its own control, in every study that tested this?

**Answer:** `{"induced_everywhere": "no", "log2_change": {"GSE117942": -2.52, "GSE4025": -0.07}}`

```sql
SELECT st.geo_accession, round(e.log2_change, 2) AS log2_change, e.n_treated, e.n_control, e.treated_ids, e.control_ids
FROM treatment_effect('GREB1') e JOIN study st USING (study_id)
JOIN (SELECT DISTINCT cell_line_id, derived_line_id, genetic_change FROM sample s WHERE s.cell_line_id = 'CVCL_0031' AND s.derived_line_id IS NULL AND s.genetic_change IN ('none', 'empty_vector', 'non_targeting_control')) o
  ON o.cell_line_id = e.cell_line_id AND o.derived_line_id IS NOT DISTINCT FROM e.derived_line_id
 AND o.genetic_change = e.genetic_change
WHERE e.compounds = 'OHT' AND e.exposure_h = 24 AND e.log2_change IS NOT NULL
ORDER BY 1
```

## Q19

*(Scientist)* In GSE117942, which compound raised GREB1 in every cell line, compared with DMSO from the same run?

**Answer:** `["E2"]`

```sql
-- compounds that raised GREB1 in every GSE117942 cell line, using only same-run DMSO controls;
-- a compound with any line lacking a same-run control can't qualify
WITH x AS (SELECT compounds, cell_line_id, log2_change FROM treatment_effect('GREB1') WHERE study_id = 'S5'),
lines AS (SELECT count(DISTINCT cell_line_id) AS n FROM sample WHERE study_id = 'S5')
SELECT compounds, count(*) FILTER (WHERE log2_change > 0) AS lines_up, count(*) AS lines_tested,
       count(*) FILTER (WHERE log2_change IS NULL) AS lines_without_same_run_control, (SELECT n FROM lines) AS lines_in_study
FROM x GROUP BY 1 ORDER BY 1
```

## Q20

*(Scientist)* Does 24 hours of estradiol (E2) raise GREB1 in ordinary MCF-7 in every study that tested it, compared within each experimental run?

**Answer:** `{"answer": "yes", "log2_change": {"GSE117942 S5_SAM24314": 2.85, "GSE21618 S2_wt1": 0.48, "GSE21618 S2_wt3": 1.33, "GSE21618 S2_wt4": 0.29, "GSE4025": 1.77}}`

```sql
SELECT st.geo_accession, e.batch, round(e.log2_change, 2) AS log2_change
FROM treatment_effect('GREB1') e JOIN study st USING (study_id)
JOIN (SELECT DISTINCT cell_line_id, derived_line_id, genetic_change FROM sample s WHERE s.cell_line_id = 'CVCL_0031' AND s.derived_line_id IS NULL AND s.genetic_change IN ('none', 'empty_vector', 'non_targeting_control')) o
  ON o.cell_line_id = e.cell_line_id AND o.derived_line_id IS NOT DISTINCT FROM e.derived_line_id
 AND o.genetic_change = e.genetic_change
WHERE e.compounds = 'E2' AND e.exposure_h = 24 AND e.log2_change IS NOT NULL
ORDER BY 1, 2
```

## Q21

*(Scientist)* In GSE26298, how much lower is ESR1 in the ESR1-knockdown samples than in their own experiment's control?

**Answer:** `-2.82`

```sql
WITH x AS (
    SELECT s.batch, s.genetic_change, s.genetic_change_target, m.value
    FROM sample s JOIN measurement m USING (sample_id) JOIN gene g USING (hgnc_id)
    WHERE g.symbol = 'ESR1' AND s.study_id = 'S3' AND s.role = 'vehicle_control'),
kd AS (SELECT batch, avg(value) AS v FROM x WHERE genetic_change = 'knockdown'
       AND genetic_change_target = (SELECT hgnc_id FROM gene WHERE symbol = 'ESR1') GROUP BY 1),
ctl AS (SELECT batch, avg(value) AS v FROM x WHERE genetic_change = 'non_targeting_control' GROUP BY 1)
SELECT kd.batch, round(kd.v - ctl.v, 2) AS log2_change FROM kd JOIN ctl USING (batch)
```

## Q22

*(Scientist)* In GSE111151, which tamoxifen-resistant line shows the largest drop in ESR1 compared with its parental line?

**Answer:** `"ZR-75-1 Tam2"`

```sql
WITH x AS (
    SELECT s.sample_id, s.cell_line_id, s.derived_line_id, m.value
    FROM sample s JOIN measurement m USING (sample_id) JOIN gene g USING (hgnc_id)
    WHERE g.symbol = 'ESR1' AND s.study_id = 'S4')
SELECT cl.preferred_name AS resistant_line, round(r.value - p.value, 2) AS log2_change_vs_parent
FROM x r JOIN x p ON p.cell_line_id = r.cell_line_id AND p.derived_line_id IS NULL
JOIN cell_line cl ON cl.cell_line_id = r.derived_line_id
WHERE r.derived_line_id IS NOT NULL
ORDER BY log2_change_vs_parent
```

## Q23

*(CEO)* Which studies can't simply be pooled with the others to study 4-OHT response, and why?

**Answer:** `{"GSE111151": ["4-OHT only as months-long selection, not an acute treatment", "cells not hormone-deprived, unlike the other studies", "resistant lines in drug at harvest vs parental lines without it; no replicates"], "GSE21618": ["tamoxifen form (and dose) not stated, so it may not be 4-OHT"], "GSE26298": ["4-OHT exposure time not stated"], "all studies": ["array and RNA-seq values are on different scales; only within-study changes can be combined"]}`

```sql
-- study-level facts that block pooling for a 4-OHT question
SELECT st.geo_accession,
       bool_or(t.compound_id = 'TAM_UNSPECIFIED')                       AS tamoxifen_form_not_stated,
       bool_or(t.compound_id = 'OHT' AND t.exposure_type = 'acute' AND t.exposure_h IS NULL) AS oht_exposure_not_stated,
       bool_or(t.compound_id = 'OHT' AND t.exposure_type = 'chronic')    AS oht_chronic_only,
       bool_or(s.hormone_deprived = 'no')                                AS not_hormone_deprived,
       bool_or(s.in_drug_at_harvest = 'yes' AND s.role = 'treated' AND t.exposure_type = 'chronic') AS drug_present_vs_absent,
       count(DISTINCT s.batch) FILTER (WHERE s.batch_basis = 'inferred')  AS inferred_batches,
       any_value(st.technology)                                          AS technology
FROM study st JOIN sample s USING (study_id) LEFT JOIN sample_treatment t USING (sample_id)
GROUP BY 1 ORDER BY 1
```

## Q24

*(Head of data)* Which studies measured ordinary MCF-7 given 4-OHT for a stated 24 hours with a matched control in the same run, so that their responses can be combined?

**Answer:** `["GSE117942", "GSE4025"]`

```sql
SELECT DISTINCT st.geo_accession
FROM treatment_effect('GREB1') e JOIN study st USING (study_id)
JOIN (SELECT DISTINCT cell_line_id, derived_line_id, genetic_change FROM sample s WHERE s.cell_line_id = 'CVCL_0031' AND s.derived_line_id IS NULL AND s.genetic_change IN ('none', 'empty_vector', 'non_targeting_control')) o
  ON o.cell_line_id = e.cell_line_id AND o.derived_line_id IS NOT DISTINCT FROM e.derived_line_id
 AND o.genetic_change = e.genetic_change
WHERE e.compounds = 'OHT' AND e.exposure_type = 'acute' AND e.exposure_h = 24 AND e.n_control > 0
ORDER BY 1
```

## Q25

*(CEO)* Combining only the studies that can be combined, what is the average change in GREB1 after 24 hours of 4-OHT in ordinary MCF-7?

**Answer:** `-1.29`

```sql
WITH per_study AS (
    SELECT e.study_id, avg(e.log2_change) AS change
    FROM treatment_effect('GREB1') e
    JOIN (SELECT DISTINCT cell_line_id, derived_line_id, genetic_change FROM sample s WHERE s.cell_line_id = 'CVCL_0031' AND s.derived_line_id IS NULL AND s.genetic_change IN ('none', 'empty_vector', 'non_targeting_control')) o
      ON o.cell_line_id = e.cell_line_id AND o.derived_line_id IS NOT DISTINCT FROM e.derived_line_id
     AND o.genetic_change = e.genetic_change
    WHERE e.compounds = 'OHT' AND e.exposure_type = 'acute' AND e.exposure_h = 24 AND e.n_control > 0
    GROUP BY 1)
SELECT round(avg(change), 2) AS mean_log2_change, count(*) AS studies,
       string_agg(study_id || ': ' || CAST(round(change, 2) AS VARCHAR), '; ' ORDER BY study_id) AS per_study
FROM per_study
```
