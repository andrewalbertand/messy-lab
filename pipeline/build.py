"""
Messy Lab: one command builds everything.

    python3 pipeline/build.py

Runs every step in order, from the raw downloads to the finished database:
  0. check_reference.py  the reference lists load and every spelling appears in its raw file
  1. loaders.py          raw files -> db/staging.duckdb (raw files proven unchanged)
  2. standardize.py      sample labels -> clean samples and treatments, with provenance
  3. genes.py            every study's row IDs -> HGNC gene IDs
  4. values.py           every value -> one scale per technology
  5. checks.py           every problem found -> the issue table
  6. this file           everything -> db/messy_lab.duckdb, built fresh from pipeline/schema.sql
  7. coverage.py        db/coverage.md: what was loaded, what was found, how all 75 inventory rows are handled

The finished database is deleted and rebuilt on every run, so it only ever reflects the raw files
and the code. The run ends by printing a fingerprint of every table's contents: delete the
database, build again, and the fingerprint is identical.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

import duckdb
import yaml

ROOT = Path(__file__).resolve().parent.parent
DB_DIR = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db"))
STAGING = DB_DIR / "staging.duckdb"
FINAL = DB_DIR / "messy_lab.duckdb"
REF = ROOT / "pipeline" / "reference"
STEPS = ["reference/check_reference.py", "loaders.py", "standardize.py", "genes.py", "values.py", "checks.py"]
TABLES = ["source_file", "study", "compound", "compound_alias", "cell_line", "cell_line_alias", "gene",
          "feature_map", "sample", "sample_treatment", "measurement", "provenance", "issue"]


def run_steps():
    for step in STEPS:
        t = time.time()
        print(f"\n=== {step} " + "=" * (80 - len(step)))
        subprocess.run([sys.executable, str(ROOT / "pipeline" / step)], cwd=ROOT, check=True)
        print(f"--- {step} done in {time.time() - t:.1f}s")


def assemble():
    for old in (FINAL, FINAL.with_name(FINAL.name + ".wal")):
        if old.exists():
            old.unlink()
    con = duckdb.connect(str(FINAL))
    con.execute((ROOT / "pipeline" / "schema.sql").read_text())
    con.execute(f"ATTACH '{STAGING}' AS stg (READ_ONLY)")

    con.execute("""INSERT INTO source_file
                   SELECT source_file_id, path, study_id, file_kind, sha256, size_bytes, origin_url,
                          DATE '2026-09-26', NULL FROM stg.stg_source_file ORDER BY source_file_id""")

    facts = yaml.safe_load(open(REF / "study.yaml"))
    for sid, f in facts.items():
        geo = dict(con.execute("SELECT field, value FROM stg.stg_series WHERE study_id=? AND n=1", [sid]).fetchall())
        rel = [v for (v,) in con.execute("SELECT value FROM stg.stg_series WHERE study_id=? AND field='Series_relation' "
                                         "AND value LIKE 'SubSeries of%'", [sid]).fetchall()]
        proc = " | ".join(v for (v,) in con.execute(
            "SELECT DISTINCT value FROM stg.stg_sample_field WHERE study_id=? AND field LIKE 'data_processing%' ORDER BY value", [sid]).fetchall())
        n = con.execute("SELECT count(DISTINCT sample_id) FROM stg.stg_sample_field WHERE study_id=?", [sid]).fetchone()[0]
        con.execute("INSERT INTO study VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", [
            sid, geo["Series_geo_accession"], geo["Series_title"], f["paper_title"], f["first_author"], f["year"],
            f["journal"], f["pmid"], f["lab_institution"], f["technology"], geo["Series_platform_id"],
            f["platform_description"], proc, f["raw_value_scale"], n,
            rel[0].split(":")[-1].strip() if rel else None, f["notes"]])

    csv = lambda t: f"read_csv('{REF / (t + '.csv')}', header=true, nullstr='')"
    con.execute(f"INSERT INTO compound SELECT * FROM {csv('compound')} ORDER BY compound_id")
    con.execute(f"INSERT INTO compound_alias SELECT * FROM {csv('compound_alias')} ORDER BY alias, seen_in_study")
    con.execute(f"INSERT INTO cell_line SELECT * FROM {csv('cell_line')} WHERE NOT is_derived ORDER BY cell_line_id")
    con.execute(f"INSERT INTO cell_line SELECT * FROM {csv('cell_line')} WHERE is_derived ORDER BY cell_line_id")
    con.execute(f"INSERT INTO cell_line_alias SELECT * FROM {csv('cell_line_alias')} ORDER BY alias, seen_in_study")

    con.execute("INSERT INTO gene SELECT * FROM stg.std_gene ORDER BY hgnc_id")
    con.execute("INSERT INTO feature_map SELECT * FROM stg.std_feature_map ORDER BY study_id, feature_id")
    con.execute("INSERT INTO sample SELECT * FROM stg.std_sample ORDER BY sample_id")
    con.execute("INSERT INTO sample_treatment SELECT * FROM stg.std_sample_treatment ORDER BY sample_id, compound_id")
    con.execute("INSERT INTO measurement SELECT * FROM stg.std_measurement ORDER BY sample_id, hgnc_id")
    con.execute("""INSERT INTO provenance (entity_type, entity_id, field, clean_value, source_file_id, source_location, raw_text, rule)
                   SELECT * FROM stg.std_provenance
                   ORDER BY entity_type, entity_id, field, rule, source_location, clean_value""")
    con.execute("""INSERT INTO issue (check_name, inventory_id, impact, study_id, entity_type, entity_id, field, observed,
                                      expected, resolution, message)
                   SELECT check_name, inventory_id, impact, study_id, entity_type, entity_id, field, observed, expected,
                          resolution, message FROM stg.std_issue
                   ORDER BY check_name, study_id, entity_type, entity_id, field, observed""")
    con.execute("DETACH stg")
    return con


def fingerprint(con):
    """Order-independent hash of every table's contents (the issue timestamp is left out)."""
    out = {}
    for t in TABLES:
        cols = "* EXCLUDE (created_at)" if t == "issue" else "*"
        n, h = con.execute(f"SELECT count(*), md5(coalesce(string_agg(md5(CAST(x AS VARCHAR)), '' ORDER BY md5(CAST(x AS VARCHAR))), '')) "
                           f"FROM (SELECT {cols} FROM {t}) x").fetchone()
        out[t] = (n, h)
    overall = __import__("hashlib").md5("".join(h for _, h in out.values()).encode()).hexdigest()
    return out, overall


def main():
    start = time.time()
    run_steps()
    print("\n=== assemble db/messy_lab.duckdb " + "=" * 50)
    con = assemble()
    tables, overall = fingerprint(con)
    con.close()
    lines = [f"{t:18} {n:>12,} rows   {h}" for t, (n, h) in tables.items()]
    report = "\n".join(lines) + f"\n\nfingerprint {overall}\n"
    (ROOT / "db" / "build_fingerprint.txt").write_text(
        "# Contents of db/messy_lab.duckdb from the last build. Rebuilding from the same raw files and code\n"
        "# gives the same fingerprint.\n" + report)
    print(report)
    print("\n=== coverage.py " + "=" * 67)
    subprocess.run([sys.executable, str(ROOT / "pipeline" / "coverage.py")], cwd=ROOT, check=True)
    print(f"Built {FINAL} in {time.time() - start:.0f}s.")


if __name__ == "__main__":
    main()
