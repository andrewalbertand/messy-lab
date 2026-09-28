"""
Messy Lab, step 2: loaders.

Reads the five studies' raw downloads exactly as GEO delivered them and lands everything in one
staging database, db/staging.duckdb, in the same four shapes for every study. Nothing is cleaned,
renamed or interpreted here. That happens in the mapping step. The raw files are opened read-only
and their checksums are compared before and after, so the run proves it never changed them.

There is one loader per file layout:
  * series matrix with values in it ............ S1, S2, S3
  * labels + a tar bundle of per-sample files ... S4
  * labels + one wide counts table .............. S5 (columns matched by sample ID, never position)

Staging tables:
  stg_source_file   every file used, with sha256 and size
  stg_series        every !Series_ line:  study_id, field, n, value
  stg_sample_field  every !Sample_ line, one row per sample per field:  study_id, sample_id, field, value
  stg_feature       each study's row IDs with whatever gene annotation the study supplied
  stg_value         every raw number:  study_id, sample_id, feature_id, raw_value

Run from the Messy_Lab folder:   python3 pipeline/loaders.py
"""
import csv
import gzip
import os
import hashlib
import io
import re
import tarfile
from pathlib import Path

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "studies" / "data" / "raw"
STAGING = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db")) / "staging.duckdb"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo"

# source_file_id: (path relative to Messy_Lab, study_id, file_kind, origin_url)
FILES = {
    "S1_matrix":   ("studies/data/raw/paper_1_GSE4025_series_matrix.txt.gz", "S1", "series_matrix", f"{GEO}/series/GSE4nnn/GSE4025/matrix/GSE4025_series_matrix.txt.gz"),
    "S1_platform": ("studies/data/raw/paper_1_GPL96.annot.gz", "S1", "platform_annotation", f"{GEO}/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"),
    "S2_matrix":   ("studies/data/raw/paper_2_GSE21618_series_matrix.txt.gz", "S2", "series_matrix", f"{GEO}/series/GSE21nnn/GSE21618/matrix/GSE21618_series_matrix.txt.gz"),
    "S2_platform": ("studies/data/raw/paper_2_GPL10371.txt", "S2", "platform_annotation", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10371&targ=self&form=text&view=data"),
    "S3_matrix":   ("studies/data/raw/paper_3_GSE26298_series_matrix.txt.gz", "S3", "series_matrix", f"{GEO}/series/GSE26nnn/GSE26298/matrix/GSE26298_series_matrix.txt.gz"),
    "S3_platform": ("studies/data/raw/paper_3_GPL570.annot.gz", "S3", "platform_annotation", f"{GEO}/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"),
    "S4_matrix":   ("studies/data/raw/paper_4_GSE111151_series_matrix.txt.gz", "S4", "series_matrix", f"{GEO}/series/GSE111nnn/GSE111151/matrix/GSE111151_series_matrix.txt.gz"),
    "S4_raw_tar":  ("studies/data/raw/paper_4_GSE111151_RAW.tar", "S4", "per_sample_bundle", f"{GEO}/series/GSE111nnn/GSE111151/suppl/GSE111151_RAW.tar"),
    "S5_matrix":   ("studies/data/raw/paper_5_GSE117942_series_matrix.txt.gz", "S5", "series_matrix", f"{GEO}/series/GSE117nnn/GSE117942/matrix/GSE117942_series_matrix.txt.gz"),
    "S5_counts":   ("studies/data/raw/paper_5_GSE117942_RNA-seq_counts.tsv.gz", "S5", "counts_table", f"{GEO}/series/GSE117nnn/GSE117942/suppl/GSE117942_RNA-seq_counts.tsv.gz"),
    "HGNC_2026-09-26": ("studies/data/reference/hgnc_complete_set_2026-09-26.txt", None, "gene_reference", "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"),
    "S1_paper": ("studies/papers/paper_1.pdf", "S1", "paper", "https://pubmed.ncbi.nlm.nih.gov/16849584/"),
    "S2_paper": ("studies/papers/paper_2.pdf", "S2", "paper", "https://pubmed.ncbi.nlm.nih.gov/21044952/"),
    "S3_paper": ("studies/papers/paper_3.pdf", "S3", "paper", "https://pubmed.ncbi.nlm.nih.gov/21299862/"),
    "S4_paper": ("studies/papers/paper_4.pdf", "S4", "paper", "https://pubmed.ncbi.nlm.nih.gov/30143015/"),
    "S5_paper": ("studies/papers/paper_5.pdf", "S5", "paper", "https://pubmed.ncbi.nlm.nih.gov/31353221/"),
}
EXPECTED_SAMPLES = {"S1": 17, "S2": 143, "S3": 12, "S4": 11, "S5": 82}

# The five papers are copyrighted and are not in the public repository (studies/papers/README.md links them).
# Nothing in the build reads them; they are recorded as sources. If a PDF is present its checksum must match;
# if it is missing, the recorded checksum and size are used, so the build and its fingerprint are identical.
PAPERS = {
    "S1_paper": ("794bd26eb43d824207ac6bf9bba72e83495f4f3a8723b79709d02de58d590916", 1285837),
    "S2_paper": ("9f4d73b323f1ac3ddb829913449f9da7f69635d36f9fb473a36bf6729290387f", 3598718),
    "S3_paper": ("33bf4f3807a528629439705233d5cba78dfb420f8c81f406a27a6555367b1bc0", 1799545),
    "S4_paper": ("e790de657d271904a5a9f9217fa02a89962195dc0e2c4eb75ec44af6fedbeef4", 5942584),
    "S5_paper": ("cf0b5a2724cfd715badad5c7810c0bf6261005d480553a57b6d3480076884e86", 11259850),
}


def checksum(fid):
    """sha256 of a source file; for a paper that isn't downloaded, its recorded checksum."""
    p = path_of(fid)
    if fid in PAPERS and not p.exists():
        return PAPERS[fid][0]
    h = sha256(p)
    assert fid not in PAPERS or h == PAPERS[fid][0], f"{p} is not the expected version of the paper"
    return h


def size_of(fid):
    p = path_of(fid)
    return PAPERS[fid][1] if fid in PAPERS and not p.exists() else p.stat().st_size


# ----------------------------------------------------------------------------------------------
# shared readers
# ----------------------------------------------------------------------------------------------
def path_of(file_id):
    return ROOT / FILES[file_id][0]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_matrix(file_id):
    """Split a GEO series matrix into series lines, sample lines and the value table (if any)."""
    series, sample_rows, table = [], [], []
    in_table = False
    with gzip.open(path_of(file_id), "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table.append(line)
            elif line.startswith("!Series_"):
                key, _, val = line.partition("\t")
                series.append((key[1:], val.strip('"')))
            elif line.startswith("!Sample_"):
                sample_rows.append(next(csv.reader([line], delimiter="\t")))
    return series, sample_rows, table


def sample_fields(study_id, sample_rows):
    """Long table of every !Sample_ field. Repeated keys (e.g. characteristics) get #1, #2 ..."""
    ids = next(r for r in sample_rows if r[0] == "!Sample_geo_accession")[1:]
    seen, out = {}, []
    for r in sample_rows:
        key = r[0][8:]
        seen[key] = seen.get(key, -1) + 1
        field = key if seen[key] == 0 else f"{key}#{seen[key]}"
        vals = r[1:]
        assert len(vals) == len(ids), f"{study_id} {field}: {len(vals)} values for {len(ids)} samples"
        out += [(study_id, sid, field, v) for sid, v in zip(ids, vals)]
    return pd.DataFrame(out, columns=["study_id", "sample_id", "field", "value"]), ids


def series_frame(study_id, series):
    n = {}
    rows = []
    for field, val in series:
        n[field] = n.get(field, 0) + 1
        rows.append((study_id, field, n[field], val))
    return pd.DataFrame(rows, columns=["study_id", "field", "n", "value"])


def read_gpl_table(file_id):
    lines, on = [], False
    opener = gzip.open if str(path_of(file_id)).endswith(".gz") else open
    with opener(path_of(file_id), "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                on = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if on:
                lines.append(line)
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str, keep_default_na=False)


FEATURE_COLS = ["study_id", "feature_id", "feature_type", "symbol_as_given", "entrez_as_given",
                "ensembl_as_given", "refseq_as_given", "description_as_given", "extra", "source_file_id"]


# ----------------------------------------------------------------------------------------------
# layout 1: series matrix with values (S1, S2, S3)
# ----------------------------------------------------------------------------------------------
def load_matrix_study(study_id):
    series, sample_rows, table = read_matrix(f"{study_id}_matrix")
    fields, ids = sample_fields(study_id, sample_rows)
    values = pd.read_csv(io.StringIO("\n".join(table)), sep="\t", dtype={"ID_REF": str})
    values = values.rename(columns={"ID_REF": "feature_id"})
    assert list(values.columns[1:]) == ids, f"{study_id}: value columns don't match sample IDs"

    gpl = read_gpl_table(f"{study_id}_platform")
    if study_id in ("S1", "S3"):          # standard Affymetrix annotation
        feats = pd.DataFrame({
            "feature_id": gpl["ID"], "feature_type": "affy_probe_set",
            "symbol_as_given": gpl["Gene symbol"], "entrez_as_given": gpl["Gene ID"],
            "ensembl_as_given": "", "refseq_as_given": "", "description_as_given": gpl["Gene title"], "extra": ""})
    else:                                  # S2: custom RefSeq layout; symbol only inside free text
        feats = pd.DataFrame({
            "feature_id": gpl["ID"], "feature_type": "refseq_probe_set",
            "symbol_as_given": "", "entrez_as_given": "", "ensembl_as_given": "",
            "refseq_as_given": gpl["GB_ACC"], "description_as_given": gpl["Description"], "extra": ""})
    feats.insert(0, "study_id", study_id)
    feats["source_file_id"] = f"{study_id}_platform"
    return series_frame(study_id, series), fields, feats[FEATURE_COLS], values, f"{study_id}_matrix"


# ----------------------------------------------------------------------------------------------
# layout 2: labels + tar bundle of per-sample files (S4)
# ----------------------------------------------------------------------------------------------
def load_s4():
    series, sample_rows, _ = read_matrix("S4_matrix")
    fields, ids = sample_fields("S4", sample_rows)
    per_sample, feats, skipped = {}, None, []
    with tarfile.open(path_of("S4_raw_tar"), "r") as tar:
        for m in tar.getmembers():
            gsm = m.name.split("_")[0]
            if gsm not in ids:                    # anything not a sample of this series stays out
                skipped.append(m.name)
                continue
            d = pd.read_csv(io.BytesIO(gzip.decompress(tar.extractfile(m).read())), sep="\t", dtype={"EnsEMBL_GenID": str, "gene_name": str})
            assert list(d.columns) == ["EnsEMBL_GenID", "gene_name", "counts", "CPM (log2)", "CPM_batch (log2)"], m.name
            if feats is None:
                feats = d[["EnsEMBL_GenID", "gene_name"]]
            else:
                assert d["EnsEMBL_GenID"].equals(feats["EnsEMBL_GenID"]), f"{m.name}: gene rows differ"
            # raw counts only: the CPM columns were normalized together with 8 patient samples (M59)
            per_sample[gsm] = d["counts"].values
    assert set(per_sample) == set(ids), f"S4: samples without a file: {set(ids) - set(per_sample)}"
    values = pd.DataFrame({"feature_id": feats["EnsEMBL_GenID"], **{g: per_sample[g] for g in ids}})
    f = pd.DataFrame({"study_id": "S4", "feature_id": feats["EnsEMBL_GenID"], "feature_type": "ensembl_gene",
                      "symbol_as_given": feats["gene_name"], "entrez_as_given": "", "ensembl_as_given": feats["EnsEMBL_GenID"],
                      "refseq_as_given": "", "description_as_given": "", "extra": "", "source_file_id": "S4_raw_tar"})
    # The 8 patient tumor files (GSE111151_GSM1417177..184) were never downloaded, and any other
    # member that isn't one of this series' 11 samples is skipped above.
    return series_frame("S4", series), fields, f[FEATURE_COLS], values, "S4_raw_tar", skipped


# ----------------------------------------------------------------------------------------------
# layout 3: labels + one wide counts table with extra header rows (S5)
# ----------------------------------------------------------------------------------------------
def load_s5():
    series, sample_rows, _ = read_matrix("S5_matrix")
    fields, ids = sample_fields("S5", sample_rows)
    title_to_gsm = dict(zip(fields[fields.field == "title"].value, fields[fields.field == "title"].sample_id))

    raw = pd.read_csv(path_of("S5_counts"), sep="\t", header=None, dtype=str, keep_default_na=False)
    head_cell, head_treat, header = raw.iloc[0], raw.iloc[1], raw.iloc[2]
    assert head_cell[10] == "CELL LINE" and head_treat[10] == "TREATMENT", "S5: header rows not where expected"
    annot_cols = list(header[:11])            # aliases chr start end str coords size type symbol desc Entrez Gene ID
    sam_codes = list(header[11:])
    # match by sample ID, never by position
    missing = [s for s in sam_codes if s not in title_to_gsm]
    assert not missing and len(sam_codes) == len(ids), f"S5: counts columns without a matching sample: {missing}"
    gsm_of_col = [title_to_gsm[s] for s in sam_codes]
    body = raw.iloc[3:].reset_index(drop=True)
    body.columns = annot_cols + gsm_of_col

    values = body[["Entrez Gene ID"] + ids].rename(columns={"Entrez Gene ID": "feature_id"})   # reordered to matrix order by ID
    values[ids] = values[ids].astype("int64")

    extra = body[["aliases", "chr", "start", "end", "str", "size", "type"]].apply(
        lambda r: "|".join(f"{k}={v}" for k, v in r.items()), axis=1)
    f = pd.DataFrame({"study_id": "S5", "feature_id": body["Entrez Gene ID"], "feature_type": "entrez_gene",
                      "symbol_as_given": body["symbol"], "entrez_as_given": body["Entrez Gene ID"], "ensembl_as_given": "",
                      "refseq_as_given": "", "description_as_given": body["desc"], "extra": extra, "source_file_id": "S5_counts"})

    # the counts file's own two label rows become sample fields too, so both label sources are kept
    extra_fields = pd.DataFrame(
        [("S5", g, "counts_file_cell_line", c) for g, c in zip(gsm_of_col, head_cell[11:])] +
        [("S5", g, "counts_file_treatment", t) for g, t in zip(gsm_of_col, head_treat[11:])] +
        [("S5", g, "counts_file_column_position", str(i + 12)) for i, g in enumerate(gsm_of_col)],
        columns=["study_id", "sample_id", "field", "value"])
    fields = pd.concat([fields, extra_fields], ignore_index=True)
    return series_frame("S5", series), fields, f[FEATURE_COLS], values, "S5_counts"


# ----------------------------------------------------------------------------------------------
# run
# ----------------------------------------------------------------------------------------------
def main():
    before = {fid: checksum(fid) for fid in FILES}

    STAGING.parent.mkdir(exist_ok=True)
    for old in (STAGING, STAGING.with_name(STAGING.name + ".wal")):   # always rebuild from scratch
        if old.exists():
            old.unlink()
    con = duckdb.connect(str(STAGING))

    src = pd.DataFrame([(fid, p, s, k, before[fid], size_of(fid), url)
                        for fid, (p, s, k, url) in FILES.items()],
                       columns=["source_file_id", "path", "study_id", "file_kind", "sha256", "size_bytes", "origin_url"])
    con.execute("CREATE TABLE stg_source_file AS SELECT * FROM src")

    results = [load_matrix_study("S1"), load_matrix_study("S2"), load_matrix_study("S3")]
    s4 = load_s4()
    results.append(s4[:5])
    results.append(load_s5())

    first = True
    for series, fields, feats, values, value_src in results:
        study = fields.study_id.iloc[0]
        wide = values
        con.register("wide", wide)
        cols = ", ".join(f'"{c}"' for c in wide.columns[1:])
        long_sql = (f"SELECT '{study}' AS study_id, sample_id, feature_id, CAST(raw_value AS DOUBLE) AS raw_value, "
                    f"'{value_src}' AS source_file_id "
                    f"FROM (UNPIVOT wide ON {cols} INTO NAME sample_id VALUE raw_value)")
        for name, df in (("stg_series", series), ("stg_sample_field", fields), ("stg_feature", feats)):
            con.register("df", df)
            con.execute(f"CREATE TABLE {name} AS SELECT * FROM df" if first else f"INSERT INTO {name} SELECT * FROM df")
            con.unregister("df")
        con.execute(f"CREATE TABLE stg_value AS {long_sql}" if first else f"INSERT INTO stg_value {long_sql}")
        con.unregister("wide")
        first = False

    # ---- checks -------------------------------------------------------------------------------
    after = {fid: checksum(fid) for fid in FILES}
    changed = [f for f in FILES if before[f] != after[f]]
    assert not changed, f"raw files changed during the run: {changed}"

    print(f"Staging written to {STAGING}. Raw files unchanged (sha256 before = after for all "
          f"{len(FILES)} files).\n")
    print(f"{'study':6} {'samples':>8} {'expected':>8} {'features':>9} {'values':>11} {'blank values':>13} {'value IDs not in feature list':>30}")
    for s in ["S1", "S2", "S3", "S4", "S5"]:
        n_samp = con.execute("select count(distinct sample_id) from stg_sample_field where study_id=?", [s]).fetchone()[0]
        n_feat = con.execute("select count(*) from stg_feature where study_id=?", [s]).fetchone()[0]
        n_val, n_null = con.execute("select count(*), count(*) - count(raw_value) from stg_value where study_id=?", [s]).fetchone()
        orphan = con.execute("""select count(distinct v.feature_id) from stg_value v
                                left join stg_feature f on f.study_id=v.study_id and f.feature_id=v.feature_id
                                where v.study_id=? and f.feature_id is null""", [s]).fetchone()[0]
        print(f"{s:6} {n_samp:>8} {EXPECTED_SAMPLES[s]:>8} {n_feat:>9,} {n_val:>11,} {n_null:>13,} {orphan:>30,}")
        assert n_samp == EXPECTED_SAMPLES[s], f"{s}: {n_samp} samples, expected {EXPECTED_SAMPLES[s]}"

    dup = con.execute("select count(*) from (select study_id, sample_id, feature_id from stg_value group by all having count(*) > 1)").fetchone()[0]
    print(f"\nduplicate (sample, feature) values: {dup}")
    print(f"S4 tar members skipped as not part of this series: {s4[5] or 'none'}")
    pos = con.execute("""select count(*) from stg_sample_field a join stg_sample_field b using (study_id, sample_id)
                         where a.field='counts_file_column_position' and b.field='title'""").fetchone()[0]
    print(f"S5 counts columns matched to samples by ID: {pos} of {EXPECTED_SAMPLES['S5']}")
    print("tables:", [r[0] for r in con.execute("show tables").fetchall()])
    con.close()


if __name__ == "__main__":
    main()
