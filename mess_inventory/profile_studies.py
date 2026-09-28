"""
Messy Lab, step 1: profile the five raw studies and print the evidence behind the mess inventory.

Run from the Messy_Lab folder:
    python3 mess_inventory/profile_studies.py > mess_inventory/evidence_output.txt

Reads only the files in studies/data/raw (untouched downloads) and the HGNC reference in
studies/data/reference. Nothing is edited. Every number quoted in mess_inventory.md comes from here.
"""
import csv, gzip, io, re, tarfile
from pathlib import Path
import numpy as np
import pandas as pd

RAW = Path("studies/data/raw")
REF = Path("studies/data/reference/hgnc_complete_set_2026-09-26.txt")
DATE_LIKE = r"^\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$"
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 70); pd.set_option("display.max_rows", 200)


def header(t):
    print("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)


def sample_sheet(path):
    """All !Sample_ rows of a GEO series matrix as a table, one row per sample."""
    rows = {}
    for line in gzip.open(path, "rt", errors="replace"):
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_"):
            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            key, i = parts[0][8:], 0
            k = key
            while k in rows:
                i += 1; k = f"{key}#{i}"
            rows[k] = parts[1:]
    return pd.DataFrame(rows)


def series_fields(path, prefix):
    out = []
    for line in gzip.open(path, "rt", errors="replace"):
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith(prefix):
            out.append(line.rstrip("\n").split("\t", 1)[1].strip('"'))
    return out


def matrix_values(path):
    lines, on = [], False
    for line in gzip.open(path, "rt"):
        if line.startswith("!series_matrix_table_begin"): on = True; continue
        if line.startswith("!series_matrix_table_end"): break
        if on: lines.append(line)
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0)


def gpl_annot(path):
    lines, on = [], False
    opener = gzip.open if str(path).endswith(".gz") else open
    for line in opener(path, "rt", errors="replace"):
        if line.startswith("!platform_table_begin"): on = True; continue
        if line.startswith("!platform_table_end"): break
        if on: lines.append(line)
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", low_memory=False)


def values_summary(df, name):
    v = df.apply(pd.to_numeric, errors="coerce")
    print(f"{name}: {v.shape[0]:,} rows x {v.shape[1]} samples | min {np.nanmin(v.values):.2f} | "
          f"median {np.nanmedian(v.values):.2f} | max {np.nanmax(v.values):,.2f} | zeros {int((v == 0).sum().sum())} | "
          f"blank {int(v.isna().sum().sum())}")


def distinct(sheet, col):
    return sheet[col].value_counts(sort=False).to_dict()


# ------------------------------------------------------------------------------------------------
header("STUDY 1  GSE4025  (Frasor 2006, microarray U133A)")
s1 = sample_sheet(RAW / "paper_1_GSE4025_series_matrix.txt.gz")
print(s1[["geo_accession", "title", "characteristics_ch1", "description"]].to_string())
print("\nsource_name:", distinct(s1, "source_name_ch1"))
print("data_processing:", distinct(s1, "data_processing"))
print("contact_institute:", distinct(s1, "contact_institute"))
print("sample extract protocol (moi):", re.findall(r"multiplicity of infection \(moi\) of [^.]*", s1.extract_protocol_ch1[0]))
print("series overall design (moi):", re.findall(r"multiplicity of infection \(moi\) of [^.]*", series_fields(RAW / "paper_1_GSE4025_series_matrix.txt.gz", "!Series_overall_design")[0]))
print("characteristics moi:", sorted(set(re.findall(r"moi\d+", " ".join(s1.characteristics_ch1)))))
v1 = matrix_values(RAW / "paper_1_GSE4025_series_matrix.txt.gz"); values_summary(v1, "values")

header("STUDY 2  GSE21618  (Oyama 2011, microarray U133 Plus 2 with custom RefSeq CDF)")
s2 = sample_sheet(RAW / "paper_2_GSE21618_series_matrix.txt.gz")
cols = ["geo_accession", "title", "source_name_ch1", "characteristics_ch1#8", "characteristics_ch1#9", "description#1"]
print(s2[cols].to_string())
print("\nphenotype:", distinct(s2, "characteristics_ch1#7"))
print("treatment:", distinct(s2, "characteristics_ch1#8"))
print("time:", distinct(s2, "characteristics_ch1#9"))
print("growth protocol:", s2.growth_protocol_ch1[0])
print("data_processing:", s2.data_processing[0])
# source name vs treatment disagreement
def src_treat(s):
    parts = [p.strip() for p in s.split(",")]
    return parts[2] if len(parts) > 2 else ""
t_norm = {"untreated": "without stimulus", "E2": "E2", "E2, tamoxifen": "E2+Tamoxifen", "HRG": "HRG",
          "HRG, tamoxifen": "HRG+Tamoxifen", "tamoxifen": "Tamoxifen"}
s2["src_treat"] = s2.source_name_ch1.map(src_treat)
s2["field_treat"] = s2["characteristics_ch1#8"].str.replace("treatment: ", "").map(t_norm)
bad = s2[s2.src_treat != s2.field_treat]
print("\nSamples whose source_name treatment disagrees with the treatment field:")
print(bad[["geo_accession", "title", "source_name_ch1", "characteristics_ch1#8"]].to_string())
v2 = matrix_values(RAW / "paper_2_GSE21618_series_matrix.txt.gz"); values_summary(v2, "values")

header("STUDY 3  GSE26298  (Salazar 2011, microarray U133 Plus 2)")
s3 = sample_sheet(RAW / "paper_3_GSE26298_series_matrix.txt.gz")
print(s3[["geo_accession", "title", "characteristics_ch1#1", "characteristics_ch1#2", "characteristics_ch1#3"]].to_string())
core = s3.title.str.replace(r"\s*\(MD\d+\)", "", regex=True)
print("\ntitles that are identical once the (MDxx) code is removed:")
print(s3.assign(core=core)[core.duplicated(keep=False)][["geo_accession", "title"]].to_string())
print("source_name:", distinct(s3, "source_name_ch1"))
print("data_processing:", s3.data_processing[0])
print("contributors:", series_fields(RAW / "paper_3_GSE26298_series_matrix.txt.gz", "!Series_contributor"))
v3 = matrix_values(RAW / "paper_3_GSE26298_series_matrix.txt.gz"); values_summary(v3, "values")
lg = np.log2(v3.clip(lower=1)); lg = lg[lg.mean(axis=1) > 7]
c3 = pd.DataFrame(np.corrcoef(lg.T.values), index=s3.description_1 if "description_1" in s3 else s3["description#1"], columns=s3["description#1"])
print("\nCorrelation of the four vehicle controls with the rest (log2, expressed probes):")
print(c3.loc[["MD31", "MD32", "MD39", "MD40"]].round(3).to_string())
a = [f"MD{i}" for i in (33, 34, 35, 36, 37, 38)]; b = ["MD41", "MD42"]
for ctl in (["MD31", "MD32"], ["MD39", "MD40"]):
    print(f"mean corr {ctl} vs OH-Tam/ER-siRNA arrays {c3.loc[ctl, a].values.mean():.4f} | vs RARa-siRNA arrays {c3.loc[ctl, b].values.mean():.4f}")

header("STUDY 4  GSE111151  (Hultsch 2018, RNA-seq)")
s4 = sample_sheet(RAW / "paper_4_GSE111151_series_matrix.txt.gz")
print(s4[["geo_accession", "title", "characteristics_ch1", "description"]].to_string())
print("treatment_protocol:", distinct(s4, "treatment_protocol_ch1"))
for k in [c for c in s4.columns if c.startswith("data_processing")]:
    print(k, ":", s4[k][0][:300])
print("series supplementary files:", len(series_fields(RAW / "paper_4_GSE111151_series_matrix.txt.gz", "!Series_supplementary_file")))
print("series relations:", series_fields(RAW / "paper_4_GSE111151_series_matrix.txt.gz", "!Series_relation"))
t4 = tarfile.open(RAW / "paper_4_GSE111151_RAW.tar")
f4 = {}
for m in t4.getmembers():
    f4[m.name] = pd.read_csv(io.BytesIO(gzip.decompress(t4.extractfile(m).read())), sep="\t")
for n, d in f4.items():
    print(f"{n}: {d.shape[0]:,} genes | columns {list(d.columns)} | reads {int(d.counts.sum()):,} | duplicate gene names {d.gene_name.duplicated().sum()}")

header("STUDY 5  GSE117942  (Guan 2019, RNA-seq)")
s5 = sample_sheet(RAW / "paper_5_GSE117942_series_matrix.txt.gz")
print("cell line:", distinct(s5, "characteristics_ch1"))
print("ligand:", distinct(s5, "characteristics_ch1#1"))
print("concentration:", distinct(s5, "characteristics_ch1#2"))
print("time:", distinct(s5, "characteristics_ch1#3"))
print("treatment_protocol:", s5.treatment_protocol_ch1[0])
print("growth_protocol:", s5.growth_protocol_ch1[0])
print("data_processing#1:", s5["data_processing#1"][0])
print("series relations:", series_fields(RAW / "paper_5_GSE117942_series_matrix.txt.gz", "!Series_relation"))
print("per line x ligand:\n", pd.crosstab(s5.characteristics_ch1.str[11:], s5["characteristics_ch1#1"].str[8:]).to_string())
raw5 = pd.read_csv(RAW / "paper_5_GSE117942_RNA-seq_counts.tsv.gz", sep="\t", header=None, low_memory=False, nrows=3)
print("\nfirst three rows of the counts file (first 14 columns):\n", raw5.iloc[:, :14].to_string())
c5 = pd.read_csv(RAW / "paper_5_GSE117942_RNA-seq_counts.tsv.gz", sep="\t", header=2, low_memory=False)
lab = pd.DataFrame({"sample": raw5.iloc[2, 11:].values, "counts_file_cell": raw5.iloc[0, 11:].values, "counts_file_treatment": raw5.iloc[1, 11:].values})
m5 = s5.set_index("title")
lab["matrix_ligand"] = lab["sample"].map(m5["characteristics_ch1#1"]); lab["matrix_conc"] = lab["sample"].map(m5["characteristics_ch1#2"])
print("\ncounts-file treatment label vs series-matrix label:\n", lab.groupby(["counts_file_treatment", "matrix_ligand", "matrix_conc"]).size().to_string())
print("counts file column order equals matrix order:", list(lab["sample"]) == list(s5.title))
print(f"counts table: {c5.shape[0]:,} genes | blank symbols {c5.symbol.isna().sum()} | duplicated symbols {c5.symbol.duplicated().sum()}")
dates = c5[c5.symbol.astype(str).str.match(DATE_LIKE)][["symbol", "Entrez Gene ID", "aliases"]]
print(f"symbols turned into dates: {len(dates)}\n", dates.to_string())
lib = c5.iloc[:, 11:].sum() / 1e6
libdf = pd.DataFrame({"reads_M": lib.round(1), "cell": m5.loc[lib.index, "characteristics_ch1"].str[11:], "ligand": m5.loc[lib.index, "characteristics_ch1#1"].str[8:]})
print("\nreads per sample (millions), sorted:\n", libdf.sort_values("reads_M").to_string())

# ------------------------------------------------------------------------------------------------
header("GENE IDENTIFIERS ACROSS STUDIES")
g1 = gpl_annot(RAW / "paper_1_GPL96.annot.gz"); g3 = gpl_annot(RAW / "paper_3_GPL570.annot.gz")
for name, g in (("GPL96 (study 1)", g1), ("GPL570 (study 3)", g3)):
    sym = g["Gene symbol"].astype(str)
    print(f"{name}: {len(g):,} probes | no gene symbol {g['Gene symbol'].isna().sum():,} | map to several genes (///) {sym.str.contains('///').sum():,} | ESR1 probes {int((sym == 'ESR1').sum())}")
g2 = gpl_annot(RAW / "paper_2_GPL10371.txt")
g2 = g2[~g2.ID.str.startswith("AFFX")]
g2["sym"] = g2.Description.map(lambda d: (re.findall(r"\(([A-Za-z0-9\-\.@_/]+)\)", str(d)) or [None])[-1])
print(f"GPL10371 (study 2): {len(g2):,} probe sets | columns {list(g2.columns[:4])} | prefixes {g2.GB_ACC.str[:3].value_counts().to_dict()}")
print("  genes with more than one probe set:", int((g2.sym.value_counts() > 1).sum()), "| ESR1 rows:", g2[g2.sym == "ESR1"].ID.tolist())
print("  rows whose description is not text:", g2[g2.sym.isna()][["ID", "Description"]].to_string())

S = {
    "GSE4025": set(g1["Gene symbol"].dropna()[lambda s: ~s.str.contains("///")]),
    "GSE21618": set(g2.sym.dropna()),
    "GSE26298": set(g3["Gene symbol"].dropna()[lambda s: ~s.str.contains("///")]),
    "GSE111151": set(next(iter(f4.values())).gene_name),
    "GSE117942": set(c5.symbol.dropna()),
}
h = pd.read_csv(REF, sep="\t", low_memory=False, dtype=str)
current = set(h.symbol); prev = {}
for _, r in h[h.prev_symbol.notna()].iterrows():
    for p in r.prev_symbol.split("|"):
        prev.setdefault(p, set()).add(r.symbol)
for k, v in S.items():
    outdated = [g for g in v if g not in current and g in prev]
    print(f"{k}: {len(v):,} gene symbols | current HGNC {sum(g in current for g in v):,} | outdated (renamed since) {len(outdated):,}")
print("symbols present in all five studies, as written:", len(set.intersection(*S.values())))
for g in ["KMT2D", "MLL2", "NSD2", "WHSC1", "SEPTIN7", "SEPT7", "7-Sep", "MARCHF1", "MARCH1", "1-Mar", "ESR1", "GREB1", "PGR", "TFF1"]:
    print(f"  {g:8s}", {k: (g in v) for k, v in S.items()}, "| HGNC:", "current" if g in current else ("renamed to " + ",".join(sorted(prev[g])) if g in prev else "unknown"))
