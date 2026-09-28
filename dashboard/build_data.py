"""
Messy Lab dashboard, step 1 of 2: pull everything the dashboard shows out of db/messy_lab.duckdb
into dashboard/dashboard_data.json. Step 2 (build_page.py) puts that data into the page.

Run from the Messy_Lab folder, after pipeline/build.py:   python3 dashboard/build_data.py
"""
import json
import os
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parent.parent
DB = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db")) / "messy_lab.duckdb"
con = duckdb.connect(str(DB), read_only=True)
rows = lambda sql: [dict(zip([c[0] for c in con.description], r)) for r in con.execute(sql).fetchall()]
SETUP = (ROOT / "questions/answers.py").read_text().split('SETUP = """')[1].split('"""')[0]
con.execute(SETUP)   # treatment_effect(gene): the same matched-control rule the answer key uses

_fp = Path(DB).parent / "build_fingerprint.txt"
out = {"fingerprint": next((l.split()[1] for l in _fp.read_text().splitlines() if l.startswith("fingerprint")), "") if _fp.exists() else ""}

out["studies"] = rows("""
SELECT st.study_id, st.geo_accession, st.first_author, st.year, st.technology, st.lab_institution,
       count(*) AS samples,
       count(*) FILTER (WHERE s.role = 'treated') AS treated,
       count(*) FILTER (WHERE s.role IN ('vehicle_control', 'untreated_control')) AS controls,
       count(*) FILTER (WHERE s.role = 'parental') AS parental
FROM study st JOIN sample s USING (study_id) GROUP BY ALL ORDER BY st.study_id""")

out["cell_lines"] = rows("""
SELECT c.preferred_name AS cell_line, count(*) AS samples,
       count(*) FILTER (WHERE s.derived_line_id IS NOT NULL) AS derived,
       string_agg(DISTINCT st.geo_accession, ', ' ORDER BY st.geo_accession) AS studies
FROM sample s JOIN cell_line c USING (cell_line_id) JOIN study st USING (study_id)
GROUP BY ALL ORDER BY samples DESC, cell_line""")

out["treatments"] = rows("""
SELECT c.preferred_name AS compound, c.is_vehicle, t.compound_id,
       t.dose_nM, t.dose_status,
       CASE WHEN t.exposure_type = 'chronic' THEN 'chronic'
            WHEN t.exposure_h IS NULL THEN 'not stated'
            ELSE CAST(CAST(t.exposure_h AS INTEGER) AS VARCHAR) END AS time,
       count(DISTINCT t.sample_id) AS samples,
       string_agg(DISTINCT st.geo_accession, ', ' ORDER BY st.geo_accession) AS studies
FROM sample_treatment t JOIN compound c USING (compound_id) JOIN sample s USING (sample_id) JOIN study st USING (study_id)
GROUP BY ALL ORDER BY c.is_vehicle, compound, t.dose_nM NULLS LAST""")

out["issues"] = rows("""
SELECT st.geo_accession, i.resolution, count(*) AS n
FROM issue i LEFT JOIN study st USING (study_id) GROUP BY ALL ORDER BY 1, 2""")
out["open_issue_types"] = rows("""
SELECT st.geo_accession, i.check_name, any_value(i.inventory_id) AS inventory_id, any_value(i.impact) AS impact,
       any_value(i.message) AS example, count(*) AS n
FROM issue i LEFT JOIN study st USING (study_id) WHERE i.resolution = 'flagged'
GROUP BY ALL ORDER BY 1, n DESC""")
out["totals"] = rows("""
SELECT (SELECT count(*) FROM sample) AS samples, (SELECT count(*) FROM study) AS studies,
       (SELECT count(DISTINCT cell_line_id) FROM sample) AS cell_lines,
       (SELECT count(*) FROM compound) AS compounds,
       (SELECT count(*) FROM (SELECT hgnc_id FROM feature_map WHERE chosen_for_gene GROUP BY 1 HAVING count(DISTINCT study_id) = 5)) AS genes_all5,
       (SELECT count(*) FROM issue) AS issues,
       (SELECT count(*) FROM issue WHERE resolution = 'flagged') AS open_issues,
       (SELECT count(*) FROM measurement) AS values""")[0]

# Gene lookup: 4-OHT (or unstated tamoxifen) vs matched control, ordinary MCF-7, one comparison per study.
g = {r["study_id"] + "|" + str(r["batch"]) + "|" + str(r["genetic_change"]) + "|" + r["compounds"] + "|" + str(r["exposure_h"]): r
     for r in rows("SELECT * FROM treatment_effect('GREB1') WHERE cell_line_id = 'CVCL_0031' AND derived_line_id IS NULL")}
COMPARISONS = [
    ("GSE4025",   "S1|None|empty_vector|OHT|24.0",          "4-OHT 10 nM, 24 h, vs ethanol", None),
    ("GSE21618",  "S2|S2_wt2|none|TAM_UNSPECIFIED|24.0",    "tamoxifen 24 h vs 0 h baseline", "Form of tamoxifen never stated, so this may not be 4-OHT (M03). One treated sample."),
    ("GSE26298",  "S3|S3_A|non_targeting_control|OHT|None",  "4-OHT 100 and 500 nM vs ethanol", "Exposure time not stated (M25)."),
    ("GSE111151", None,                                     "no matched comparison", "Resistant lines were grown in 4-OHT for months; parental lines never saw it. Drug and resistance can't be separated (M35)."),
    ("GSE117942", "S5|S5_SAM24314|none|OHT|24.0",           "4-OHT 1 µM, 24 h, vs DMSO", "Dose written as \"1 M\" in GEO; corrected to 1 µM (M13)."),
]
comps, sql_parts = [], []
for i, (gse, key, label, note) in enumerate(COMPARISONS):
    r = g.get(key) if key else None
    comps.append({"study": gse, "label": label, "note": note,
                  "n_treated": r["n_treated"] if r else None, "n_control": r["n_control"] if r else None,
                  "treated": r["treated_ids"] if r else [], "control": r["control_ids"] if r else []})
    if r:
        t = ",".join(f"'{x}'" for x in r["treated_ids"]); c = ",".join(f"'{x}'" for x in r["control_ids"])
        sql_parts.append(f"round(avg(value) FILTER (WHERE sample_id IN ({t})) - avg(value) FILTER (WHERE sample_id IN ({c})), 2) AS c{i}")
    else:
        sql_parts.append(f"NULL AS c{i}")
missing = [c["study"] for c, (_, k, _, _) in zip(comps, COMPARISONS) if k and not c["treated"]]
assert not missing, f"comparison not found: {missing}"
ids = ",".join(f"'{x}'" for c in comps for x in c["treated"] + c["control"])
genes = con.execute(f"""
SELECT g.symbol, g.name, {', '.join(sql_parts)}
FROM measurement m JOIN gene g USING (hgnc_id)
WHERE m.sample_id IN ({ids})
GROUP BY ALL HAVING count(*) > 0 ORDER BY g.symbol""").fetchall()
out["comparisons"] = comps
out["genes"] = [[s, n, *[None if v is None else float(v) for v in vals]] for s, n, *vals in genes]

path = ROOT / "dashboard/dashboard_data.json"
path.write_text(json.dumps(out, default=str, separators=(",", ":")))
print(f"{path.relative_to(ROOT)}: {len(out['genes']):,} genes, {path.stat().st_size/1e6:.1f} MB")
for s in ("GREB1", "PGR", "ESR1", "TFF1"):
    print(s, next((x[2:] for x in out["genes"] if x[0] == s), None))
