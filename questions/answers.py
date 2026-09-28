"""
Messy Lab, step 4: answer the 25 questions with code against the clean database.

Every answer is computed by one SQL query against db/messy_lab.duckdb, and the query is saved next
to the answer in questions/answers.json (and shown in questions/answers.md) so anyone can re-run it.
A few questions share a helper, `treatment_effect(gene)`, defined once in SETUP and saved with the
answers: for every group of treated samples it gives the average log2 change against controls with the
same cells, same genetic change and same experimental run.

These are the machine answers. Each is then checked by hand against the raw files before it goes
into the answer key.

Run from the Messy_Lab folder, after pipeline/build.py:   python3 questions/answers.py
"""
import json
import os
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parent.parent
DB = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db")) / "messy_lab.duckdb"

ORDINARY_MCF7 = ("s.cell_line_id = 'CVCL_0031' AND s.derived_line_id IS NULL "
                 "AND s.genetic_change IN ('none', 'empty_vector', 'non_targeting_control')")

SETUP = """
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
"""

r2 = lambda x: None if x is None else round(float(x), 2)

Q = []   # (id, sql, how to turn the rows into the answer)

Q.append(("Q01", f"""
SELECT st.geo_accession, count(DISTINCT s.sample_id) AS samples
FROM sample_treatment t JOIN sample s USING (sample_id) JOIN study st USING (study_id)
WHERE t.compound_id = 'OHT' AND s.cell_line_id = 'CVCL_0031'
GROUP BY ALL ORDER BY 1""",
    lambda d: {"total": int(d.samples.sum()), "per_study": dict(zip(d.geo_accession, d.samples.astype(int)))}))

Q.append(("Q02", """
SELECT st.geo_accession, count(*) AS controls
FROM sample s JOIN study st USING (study_id)
WHERE s.role IN ('vehicle_control', 'untreated_control', 'parental')
GROUP BY ALL ORDER BY 1""",
    lambda d: dict(zip(d.geo_accession, d.controls.astype(int)))))

Q.append(("Q03", """
SELECT cl.preferred_name AS cell_line, count(*) AS samples
FROM sample s JOIN cell_line cl ON cl.cell_line_id = s.cell_line_id
GROUP BY ALL ORDER BY samples DESC, cell_line""",
    lambda d: dict(zip(d.cell_line, d.samples.astype(int)))))

Q.append(("Q04", """
SELECT c.preferred_name AS compound, count(DISTINCT t.sample_id) AS samples
FROM sample_treatment t JOIN compound c USING (compound_id)
GROUP BY ALL ORDER BY samples DESC, compound""",
    lambda d: dict(zip(d.compound, d.samples.astype(int)))))

Q.append(("Q05", """
-- resistant lines actually measured; RIKEN's time-course clone is one of its six numbered clones
-- (cell_line.notes), so it is not counted as a seventh line
SELECT st.geo_accession, cl.cell_line_id, cl.preferred_name
FROM cell_line cl JOIN study st ON st.study_id = cl.derived_by_study
WHERE cl.is_derived
  AND cl.cell_line_id IN (SELECT derived_line_id FROM sample)
  AND cl.cell_line_id <> 'S2_TAMR_SELECTED'
ORDER BY 1, 2""",
    lambda d: {"total": len(d), "per_study": {g: list(x.preferred_name) for g, x in d.groupby("geo_accession")}}))

Q.append(("Q06", """
SELECT count(*) AS genes
FROM (SELECT hgnc_id FROM feature_map WHERE chosen_for_gene
      GROUP BY 1 HAVING count(DISTINCT study_id) = 5)""",
    lambda d: int(d.genes[0])))

Q.append(("Q07", """
SELECT c.sample_id
FROM sample k JOIN sample c ON c.study_id = k.study_id AND c.batch = k.batch
WHERE k.study_id = 'S3' AND k.genetic_change = 'knockdown'
  AND k.genetic_change_target = (SELECT hgnc_id FROM gene WHERE symbol = 'RARA')
  AND c.genetic_change = 'non_targeting_control' AND c.role = 'vehicle_control'
GROUP BY 1 ORDER BY 1""",
    lambda d: list(d.sample_id)))

Q.append(("Q08", """
SELECT count(*) AS never_treated
FROM sample s
WHERE s.study_id = 'S2'
  AND NOT EXISTS (SELECT 1 FROM sample_treatment t WHERE t.sample_id = s.sample_id)""",
    lambda d: int(d.never_treated[0])))

Q.append(("Q09", """
SELECT st.geo_accession,
       CASE WHEN count(t.sample_id) FILTER (WHERE t.compound_id = 'OHT') > 0
            THEN string_agg(DISTINCT CAST(CAST(t.dose_nM AS INTEGER) AS VARCHAR), ', ' ORDER BY CAST(CAST(t.dose_nM AS INTEGER) AS VARCHAR))
                 FILTER (WHERE t.compound_id = 'OHT')
            WHEN count(t.sample_id) FILTER (WHERE t.compound_id = 'TAM_UNSPECIFIED') > 0
            THEN 'not 4-OHT: tamoxifen of unstated form, dose not stated'
            ELSE 'no 4-OHT' END AS oht_dose_nM
FROM study st LEFT JOIN sample s USING (study_id) LEFT JOIN sample_treatment t USING (sample_id)
GROUP BY 1 ORDER BY 1""",
    lambda d: dict(zip(d.geo_accession, d.oht_dose_nM))))

Q.append(("Q10", """
SELECT st.geo_accession,
       CASE WHEN bool_and(t.exposure_type = 'chronic') THEN 'chronic (months)'
            WHEN count(t.exposure_h) = 0 THEN 'not stated'
            ELSE string_agg(DISTINCT CAST(CAST(t.exposure_h AS INTEGER) AS VARCHAR), ', ') || ' h' END AS exposure,
       any_value(t.exposure_note) AS note
FROM sample_treatment t JOIN sample s USING (sample_id) JOIN study st USING (study_id)
WHERE t.compound_id = 'OHT'
GROUP BY 1 ORDER BY 1""",
    lambda d: dict(zip(d.geo_accession, d.exposure))))

Q.append(("Q11", """
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
ORDER BY 1""",
    lambda d: list(d.sample_id)))

Q.append(("Q12", """
SELECT st.geo_accession, string_agg(DISTINCT s.hormone_deprived, ',') AS hormone_deprived
FROM sample s JOIN study st USING (study_id) GROUP BY 1 ORDER BY 1""",
    lambda d: {"deprived": list(d[d.hormone_deprived == "yes"].geo_accession),
               "not_deprived": list(d[d.hormone_deprived == "no"].geo_accession)}))

Q.append(("Q13", """
SELECT DISTINCT entity_id AS sample_id FROM issue
WHERE check_name = 'fields_contradict' AND entity_type = 'sample' ORDER BY 1""",
    lambda d: list(d.sample_id)))

Q.append(("Q14", """
SELECT st.geo_accession, i.observed AS as_written, i.message
FROM issue i JOIN study st USING (study_id)
WHERE i.check_name = 'dose_above_100uM'""",
    lambda d: [{"study": g, "as_written": w.split(" ")[0] + " " + w.split(" ")[1], "corrected_nM": 1000}
               for g, w in zip(d.geo_accession, d.as_written)] or "none"))

Q.append(("Q15", """
-- mismatches against GEO's own processing text (the S3 flag compares GEO with the paper instead)
SELECT DISTINCT st.geo_accession FROM issue i JOIN study st USING (study_id)
WHERE i.check_name = 'value_scale_mismatch' AND i.expected NOT ILIKE '%paper%' ORDER BY 1""",
    lambda d: list(d.geo_accession)))

Q.append(("Q16", """
SELECT st.geo_accession, count(*) AS names FROM issue i JOIN study st USING (study_id)
WHERE i.check_name = 'date_like_gene_name' GROUP BY 1""",
    lambda d: {"count": int(d.names.sum()), "study": list(d.geo_accession)}))

Q.append(("Q17", """
SELECT DISTINCT entity_id AS sample_id FROM issue
WHERE check_name = 'duplicate_titles' AND message LIKE '%replicate number is repeated%' ORDER BY 1""",
    lambda d: list(d.sample_id)))

Q.append(("Q18", f"""
SELECT st.geo_accession, round(e.log2_change, 2) AS log2_change, e.n_treated, e.n_control, e.treated_ids, e.control_ids
FROM treatment_effect('GREB1') e JOIN study st USING (study_id)
JOIN (SELECT DISTINCT cell_line_id, derived_line_id, genetic_change FROM sample s WHERE {ORDINARY_MCF7}) o
  ON o.cell_line_id = e.cell_line_id AND o.derived_line_id IS NOT DISTINCT FROM e.derived_line_id
 AND o.genetic_change = e.genetic_change
WHERE e.compounds = 'OHT' AND e.exposure_h = 24 AND e.log2_change IS NOT NULL
ORDER BY 1""",
    lambda d: {"induced_everywhere": "yes" if (d.log2_change > 0).all() else "no",
               "log2_change": dict(zip(d.geo_accession, d.log2_change.map(r2)))}))

Q.append(("Q19", """
-- compounds that raised GREB1 in every GSE117942 cell line, using only same-run DMSO controls;
-- a compound with any line lacking a same-run control can't qualify
WITH x AS (SELECT compounds, cell_line_id, log2_change FROM treatment_effect('GREB1') WHERE study_id = 'S5'),
lines AS (SELECT count(DISTINCT cell_line_id) AS n FROM sample WHERE study_id = 'S5')
SELECT compounds, count(*) FILTER (WHERE log2_change > 0) AS lines_up, count(*) AS lines_tested,
       count(*) FILTER (WHERE log2_change IS NULL) AS lines_without_same_run_control, (SELECT n FROM lines) AS lines_in_study
FROM x GROUP BY 1 ORDER BY 1""",
    lambda d: list(d[(d.lines_up == d.lines_in_study)].compounds)))

Q.append(("Q20", f"""
SELECT st.geo_accession, e.batch, round(e.log2_change, 2) AS log2_change
FROM treatment_effect('GREB1') e JOIN study st USING (study_id)
JOIN (SELECT DISTINCT cell_line_id, derived_line_id, genetic_change FROM sample s WHERE {ORDINARY_MCF7}) o
  ON o.cell_line_id = e.cell_line_id AND o.derived_line_id IS NOT DISTINCT FROM e.derived_line_id
 AND o.genetic_change = e.genetic_change
WHERE e.compounds = 'E2' AND e.exposure_h = 24 AND e.log2_change IS NOT NULL
ORDER BY 1, 2""",
    lambda d: {"answer": "yes" if (d.log2_change > 0).all() else "no",
               "log2_change": {f"{g} {b or ''}".strip(): r2(v) for g, b, v in zip(d.geo_accession, d.batch, d.log2_change)}}))

Q.append(("Q21", """
WITH x AS (
    SELECT s.batch, s.genetic_change, s.genetic_change_target, m.value
    FROM sample s JOIN measurement m USING (sample_id) JOIN gene g USING (hgnc_id)
    WHERE g.symbol = 'ESR1' AND s.study_id = 'S3' AND s.role = 'vehicle_control'),
kd AS (SELECT batch, avg(value) AS v FROM x WHERE genetic_change = 'knockdown'
       AND genetic_change_target = (SELECT hgnc_id FROM gene WHERE symbol = 'ESR1') GROUP BY 1),
ctl AS (SELECT batch, avg(value) AS v FROM x WHERE genetic_change = 'non_targeting_control' GROUP BY 1)
SELECT kd.batch, round(kd.v - ctl.v, 2) AS log2_change FROM kd JOIN ctl USING (batch)""",
    lambda d: r2(d.log2_change[0])))

Q.append(("Q22", """
WITH x AS (
    SELECT s.sample_id, s.cell_line_id, s.derived_line_id, m.value
    FROM sample s JOIN measurement m USING (sample_id) JOIN gene g USING (hgnc_id)
    WHERE g.symbol = 'ESR1' AND s.study_id = 'S4')
SELECT cl.preferred_name AS resistant_line, round(r.value - p.value, 2) AS log2_change_vs_parent
FROM x r JOIN x p ON p.cell_line_id = r.cell_line_id AND p.derived_line_id IS NULL
JOIN cell_line cl ON cl.cell_line_id = r.derived_line_id
WHERE r.derived_line_id IS NOT NULL
ORDER BY log2_change_vs_parent""",
    lambda d: d.resistant_line[0].replace(" (Helsinki)", "")))

Q.append(("Q23", """
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
GROUP BY 1 ORDER BY 1""",
    lambda d: {g: [why for flag, why in [
        (r.tamoxifen_form_not_stated, "tamoxifen form (and dose) not stated, so it may not be 4-OHT"),
        (r.oht_exposure_not_stated, "4-OHT exposure time not stated"),
        (r.oht_chronic_only, "4-OHT only as months-long selection, not an acute treatment"),
        (r.not_hormone_deprived, "cells not hormone-deprived, unlike the other studies"),
        (r.drug_present_vs_absent, "resistant lines in drug at harvest vs parental lines without it; no replicates")]
        if flag] for g, r in zip(d.geo_accession, d.itertuples())
        if any([r.tamoxifen_form_not_stated, r.oht_exposure_not_stated, r.oht_chronic_only, r.not_hormone_deprived])}
    | {"all studies": ["array and RNA-seq values are on different scales; only within-study changes can be combined"]}))

Q.append(("Q24", f"""
SELECT DISTINCT st.geo_accession
FROM treatment_effect('GREB1') e JOIN study st USING (study_id)
JOIN (SELECT DISTINCT cell_line_id, derived_line_id, genetic_change FROM sample s WHERE {ORDINARY_MCF7}) o
  ON o.cell_line_id = e.cell_line_id AND o.derived_line_id IS NOT DISTINCT FROM e.derived_line_id
 AND o.genetic_change = e.genetic_change
WHERE e.compounds = 'OHT' AND e.exposure_type = 'acute' AND e.exposure_h = 24 AND e.n_control > 0
ORDER BY 1""",
    lambda d: list(d.geo_accession)))

Q.append(("Q25", f"""
WITH per_study AS (
    SELECT e.study_id, avg(e.log2_change) AS change
    FROM treatment_effect('GREB1') e
    JOIN (SELECT DISTINCT cell_line_id, derived_line_id, genetic_change FROM sample s WHERE {ORDINARY_MCF7}) o
      ON o.cell_line_id = e.cell_line_id AND o.derived_line_id IS NOT DISTINCT FROM e.derived_line_id
     AND o.genetic_change = e.genetic_change
    WHERE e.compounds = 'OHT' AND e.exposure_type = 'acute' AND e.exposure_h = 24 AND e.n_control > 0
    GROUP BY 1)
SELECT round(avg(change), 2) AS mean_log2_change, count(*) AS studies,
       string_agg(study_id || ': ' || CAST(round(change, 2) AS VARCHAR), '; ' ORDER BY study_id) AS per_study
FROM per_study""",
    lambda d: r2(d.mean_log2_change[0])))


def main():
    con = duckdb.connect(str(DB), read_only=True)
    con.execute(SETUP)
    texts = {}
    block = None
    for line in (ROOT / "questions" / "questions.md").read_text().splitlines():
        if line.startswith("**Q"):
            block = line[2:5]
            texts[block] = line.split("** ", 1)[1]
        elif block and line.strip() and not line.startswith("**Answer"):
            texts[block] += " " + line.strip()
        elif not line.strip():
            block = None
    out = []
    for qid, sql, fmt in Q:
        rows = con.execute(sql).df()
        ans = fmt(rows)
        out.append({"id": qid, "question": texts.get(qid, ""), "answer": ans, "query": sql.strip(),
                    "rows": json.loads(rows.to_json(orient="records", default_handler=str))})
    (ROOT / "questions" / "answers.json").write_text(json.dumps({"setup": SETUP.strip(), "answers": out}, indent=2, default=str))

    md = ["# Messy Lab · Machine answers", "",
          "Computed by `questions/answers.py` against `db/messy_lab.duckdb`. Each answer sits next to the query "
          "that produced it. These still need a hand check against the raw files before they go into the answer key.", "",
          "Shared helper used by Q18–Q20 and Q24–Q25:", "", "```sql", SETUP.strip(), "```", ""]
    for o in out:
        md += [f"## {o['id']}", "", o["question"], "", f"**Answer:** `{json.dumps(o['answer'], default=str)}`", "",
               "```sql", o["query"], "```", ""]
    (ROOT / "questions" / "answers.md").write_text("\n".join(md))
    for o in out:
        print(f"{o['id']}: {json.dumps(o['answer'], default=str)}")


if __name__ == "__main__":
    main()
