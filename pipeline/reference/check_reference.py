"""
Checks the reference lists before anything is built on them. Run from the Messy_Lab folder:

    python3 pipeline/reference/check_reference.py

1. Loads every reference CSV into the real schema, so the database's own rules check them.
2. Confirms every alias said to be seen in a study really appears in that study's raw files.
3. Confirms vocabulary.yaml and schema.sql list exactly the same allowed values.
"""
import csv, gzip, io, re, sys, tarfile
from pathlib import Path
import duckdb, yaml

ROOT = Path(".")
REF = ROOT / "pipeline/reference"
RAW = ROOT / "studies/data/raw"
problems = []

# 1. Load into the schema
con = duckdb.connect(":memory:")
con.execute((ROOT / "pipeline/schema.sql").read_text())
loads = [("compound", ""), ("compound_alias", ""),
         # cell_line points to itself, so base lines go in before the lines derived from them
         ("cell_line", " WHERE NOT is_derived"), ("cell_line", " WHERE is_derived"), ("cell_line_alias", "")]
for table, where in loads:
    try:
        con.execute(f"INSERT INTO {table} SELECT * FROM read_csv('{REF / (table + '.csv')}', header=true, nullstr=''){where}")
    except Exception as e:
        problems.append(f"{table}: schema rejected it: {e}")
for t in ["compound", "compound_alias", "cell_line", "cell_line_alias"]:
    print(f"{t}: {con.execute(f'select count(*) from {t}').fetchone()[0]} rows loaded")

# 2. Every alias appears in its study's raw files
def text_of(study):
    files = {"S1": ["paper_1_GSE4025_series_matrix.txt.gz"], "S2": ["paper_2_GSE21618_series_matrix.txt.gz"],
             "S3": ["paper_3_GSE26298_series_matrix.txt.gz"], "S4": ["paper_4_GSE111151_series_matrix.txt.gz"],
             "S5": ["paper_5_GSE117942_series_matrix.txt.gz"]}[study]
    out = ""
    for f in files:
        with gzip.open(RAW / f, "rt", errors="replace") as h:
            for line in h:
                if line.startswith("!series_matrix_table_begin"): break
                out += line
    if study == "S5":
        with gzip.open(RAW / "paper_5_GSE117942_RNA-seq_counts.tsv.gz", "rt") as h:
            out += h.readline() + h.readline()
    return out
texts = {s: text_of(s) for s in ["S1", "S2", "S3", "S4", "S5"]}
checked = 0
for name in ["compound_alias.csv", "cell_line_alias.csv"]:
    for r in csv.DictReader(open(REF / name)):
        if "paper only" in r.get("seen_in_field", ""):
            continue
        checked += 1
        if r["alias"] not in texts[r["seen_in_study"]]:
            problems.append(f"{name}: '{r['alias']}' not found in {r['seen_in_study']} raw files")
print(f"aliases checked against raw files: {checked}")

# 3. Vocabulary matches schema
vocab = yaml.safe_load(open(REF / "vocabulary.yaml"))["allowed_values"]
schema = (ROOT / "pipeline/schema.sql").read_text()
tables = dict(re.findall(r"CREATE TABLE (\w+) \((.*?)\n\);", schema, re.S))
for key, values in vocab.items():
    table, col = key.split(".")
    m = re.search(rf"\b{col}\b[^\n]*?\n?[^\n]*?IN \(([^)]*)\)", tables.get(table, ""), re.S)
    if not m:
        problems.append(f"vocabulary: no CHECK list found for {key}")
        continue
    in_schema = set(re.findall(r"'([^']+)'", m.group(1)))
    in_yaml = {str(k).lower() if isinstance(k, bool) else str(k) for k in values}
    in_yaml = {("yes" if v == "true" else "no" if v == "false" else v) for v in in_yaml}
    if in_schema != in_yaml:
        problems.append(f"vocabulary: {key} differs. schema {sorted(in_schema)} vs yaml {sorted(in_yaml)}")
print(f"vocabulary lists compared with schema: {len(vocab)}")

print("\nPROBLEMS:" if problems else "\nAll checks passed.")
for p in problems: print(" -", p)
sys.exit(1 if problems else 0)
