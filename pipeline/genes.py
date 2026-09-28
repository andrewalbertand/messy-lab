"""
Messy Lab, step 2: map every study's row IDs to HGNC gene IDs.

Why HGNC IDs: gene names change (WHSC1 is now NSD2) and Excel damages them (MARCHF1 became
'1-Mar'). An HGNC ID like HGNC:3467 never changes, so it is the only safe key for joining studies.

How each study is mapped, strongest evidence first:
  S1, S3  Affymetrix probe set -> the Entrez ID in the platform annotation -> HGNC.
          Fallback: the annotation's gene symbol, if it is a current or unambiguous previous symbol.
  S2      RefSeq accession -> HGNC's RefSeq column.
          Fallback: the symbol inside the free-text description '(A2M)', current or previous symbol.
  S4      Ensembl gene ID -> HGNC's Ensembl column.
          Fallback: the file's gene name, current or unambiguous previous symbol.
  S5      Entrez ID -> HGNC. Never the file's gene names (28 were turned into dates by Excel).

Rows that match nothing stay in feature_map as 'unmapped'; rows that match several genes are
'multiple_genes'; Affymetrix controls are 'control_probe'. Nothing is dropped silently.

When several rows map to the same gene in one study, one is chosen to carry that gene's values:
the row with the highest average expression across the study's samples (log2 scale). The rule is
written into every row's choice_rule.

Output (in db/staging.duckdb): std_gene, std_feature_map.
Run from the Messy_Lab folder, after loaders.py:   python3 pipeline/genes.py
"""
import os
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DB_DIR = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db"))
STAGING = DB_DIR / "staging.duckdb"
HGNC_FILE = ROOT / "studies/data/reference/hgnc_complete_set_2026-09-26.txt"
CHOICE_RULE = "highest mean log2 expression across this study's samples among rows mapped to this gene"


def split(v):
    return [x for x in str(v).split("|") if x and x != "nan"]


def load_hgnc():
    h = pd.read_csv(HGNC_FILE, sep="\t", dtype=str, low_memory=False)
    h = h[h.status == "Approved"]
    by_entrez = dict(zip(h.entrez_id.dropna(), h.hgnc_id[h.entrez_id.notna()]))
    by_ensembl = dict(zip(h.ensembl_gene_id.dropna(), h.hgnc_id[h.ensembl_gene_id.notna()]))
    by_refseq, by_symbol, by_prev = {}, dict(zip(h.symbol, h.hgnc_id)), {}
    for hid, rs, prev in zip(h.hgnc_id, h.refseq_accession, h.prev_symbol):
        for r in split(rs):
            by_refseq[r.split(".")[0]] = hid
        for p in split(prev):
            by_prev.setdefault(p, set()).add(hid)
    gene = pd.DataFrame({
        "hgnc_id": h.hgnc_id, "symbol": h.symbol, "name": h.name, "locus_type": h.locus_type,
        "entrez_id": h.entrez_id, "ensembl_gene_id": h.ensembl_gene_id,
        "previous_symbols": h.prev_symbol.map(split), "alias_symbols": h.alias_symbol.map(split),
        "source_file_id": "HGNC_2026-09-26"})
    return gene, by_entrez, by_ensembl, by_refseq, by_symbol, by_prev


def by_name(sym, by_symbol, by_prev):
    """Current symbol first; then a previous symbol only if it points to exactly one gene."""
    if not sym:
        return None, None
    if sym in by_symbol:
        return by_symbol[sym], "current_symbol"
    hits = by_prev.get(sym, set())
    if len(hits) == 1:
        return next(iter(hits)), "previous_symbol"
    return None, None


def map_study(feat, lookups):
    _, by_entrez, by_ensembl, by_refseq, by_symbol, by_prev = lookups
    out = []
    for r in feat.itertuples(index=False):
        fid, study = r.feature_id, r.study_id
        hid = method = None
        status = "mapped"
        if study in ("S1", "S3"):
            sym_given = r.symbol_as_given or ""
            if fid.startswith("AFFX"):
                status = "control_probe"
            elif "///" in (r.entrez_as_given or "") or "///" in sym_given:
                status = "multiple_genes"
            else:
                hid = by_entrez.get(r.entrez_as_given or "")
                method = "entrez_id" if hid else None
                if not hid:
                    hid, method = by_name(sym_given, by_symbol, by_prev)
        elif study == "S2":
            found = re.findall(r"\(([A-Za-z0-9\-\.@_/]+)\)", r.description_as_given or "")
            sym_given = found[-1] if found else ""
            if fid.startswith("AFFX"):
                status = "control_probe"
            elif set((r.description_as_given or "").strip()) == {"#"}:
                status = "corrupt_annotation"
            else:
                hid = by_refseq.get(r.refseq_as_given or "")
                method = "refseq_accession" if hid else None
                if not hid:
                    hid, method = by_name(sym_given, by_symbol, by_prev)
        elif study == "S4":
            sym_given = r.symbol_as_given or ""
            hid = by_ensembl.get(fid)
            method = "ensembl_gene_id" if hid else None
            if not hid:
                hid, method = by_name(sym_given, by_symbol, by_prev)
        else:  # S5: Entrez only, names never used
            sym_given = r.symbol_as_given or ""
            hid = by_entrez.get(fid)
            method = "entrez_id" if hid else None
        if status == "mapped" and not hid:
            status = "unmapped"
        out.append((study, fid, r.feature_type, sym_given or None, hid, status, method))
    return pd.DataFrame(out, columns=["study_id", "feature_id", "feature_type", "symbol_as_given", "hgnc_id",
                                      "map_status", "map_method"])


def main():
    con = duckdb.connect(str(STAGING))
    lookups = load_hgnc()
    gene = lookups[0]

    maps = []
    for study in ["S1", "S2", "S3", "S4", "S5"]:
        feat = con.execute("select * from stg_feature where study_id=?", [study]).df().fillna("")
        fm = map_study(feat, lookups)
        # mean expression per row, on a log2 scale, to pick one row per gene
        scale = con.execute("""select case when ? in ('S2') then 'log2' else 'linear' end""", [study]).fetchone()[0]
        expr = "avg(raw_value)" if scale == "log2" else "avg(log2(raw_value + 1))"
        means = con.execute(f"select feature_id, {expr} as m from stg_value where study_id=? group by 1", [study]).df()
        fm = fm.merge(means, on="feature_id", how="left")
        mapped = fm.map_status == "mapped"
        best = fm[mapped].sort_values(["hgnc_id", "m", "feature_id"], ascending=[True, False, True]).drop_duplicates("hgnc_id")
        fm["chosen_for_gene"] = fm.index.isin(best.index)
        fm["choice_rule"] = np.where(mapped, CHOICE_RULE, None)
        maps.append(fm.drop(columns="m"))
    FM = pd.concat(maps, ignore_index=True)

    for name, df in (("std_gene", gene), ("std_feature_map", FM)):
        con.register("df", df)
        con.execute(f"CREATE OR REPLACE TABLE {name} AS SELECT * FROM df")
        con.unregister("df")
    con.close()

    # ---- prove it fits the schema, then report --------------------------------------------------
    chk = duckdb.connect(":memory:")
    chk.execute((ROOT / "pipeline/schema.sql").read_text())
    chk.execute(f"ATTACH '{STAGING}' AS stg (READ_ONLY)")
    chk.execute("INSERT INTO source_file SELECT source_file_id, path, study_id, file_kind, sha256, size_bytes, origin_url, NULL, NULL FROM stg.stg_source_file")
    chk.execute("INSERT INTO study (study_id, geo_accession, geo_title, technology, platform, raw_value_scale, samples_in_geo) "
                "SELECT DISTINCT study_id, study_id, 'placeholder', 'microarray', 'x', 'log2', 0 FROM stg.std_feature_map")
    chk.execute("INSERT INTO gene SELECT * FROM stg.std_gene")
    chk.execute("INSERT INTO feature_map SELECT * FROM stg.std_feature_map")
    print(f"{len(gene):,} HGNC genes loaded. {len(FM):,} study rows mapped; both pass the schema's rules.\n")

    print(f"{'study':6} {'rows':>7} {'mapped':>7} {'unmapped':>9} {'multi-gene':>11} {'control':>8} {'corrupt':>8} {'genes':>7}   how they were matched")
    for s in ["S1", "S2", "S3", "S4", "S5"]:
        d = FM[FM.study_id == s]
        c = d.map_status.value_counts()
        how = d[d.map_status == "mapped"].map_method.value_counts().to_dict()
        print(f"{s:6} {len(d):>7,} {c.get('mapped',0):>7,} {c.get('unmapped',0):>9,} {c.get('multiple_genes',0):>11,} "
              f"{c.get('control_probe',0):>8,} {c.get('corrupt_annotation',0):>8,} {d[d.chosen_for_gene].shape[0]:>7,}   {how}")

    genes_by_study = [set(FM[(FM.study_id == s) & FM.chosen_for_gene].hgnc_id) for s in ["S1", "S2", "S3", "S4", "S5"]]
    print(f"\nGenes measurable in all five studies: {len(set.intersection(*genes_by_study)):,} "
          f"(by name as written it was 10,687)")

    sym = dict(zip(gene.hgnc_id, gene.symbol))
    dates = FM[(FM.study_id == "S5") & FM.symbol_as_given.fillna("").str.match(r"^\d{1,2}-(Mar|Sep|Dec)$")]
    print(f"\nS5 Excel-date names recovered through Entrez: {dates.hgnc_id.notna().sum()} of {len(dates)}")
    for r in dates.sort_values("symbol_as_given").head(6).itertuples():
        print(f"   '{r.symbol_as_given}' (Entrez {r.feature_id}) -> {r.hgnc_id} {sym.get(r.hgnc_id)}")
    print("\nRenamed genes now joined correctly:")
    for s, g in (("S2", "MLL2"), ("S1", "WHSC1"), ("S4", "SEPT7")):
        r = FM[(FM.study_id == s) & (FM.symbol_as_given == g) & FM.chosen_for_gene]
        for x in r.head(1).itertuples():
            print(f"   {s} '{g}' -> {x.hgnc_id} {sym.get(x.hgnc_id)} (via {x.map_method})")
    esr1 = FM[(FM.hgnc_id == "HGNC:3467")]
    print("\nESR1 rows per study (chosen row marked *):")
    for s, d in esr1.groupby("study_id"):
        print(f"   {s}: " + ", ".join(("*" if c else "") + f for f, c in zip(d.feature_id, d.chosen_for_gene)))


if __name__ == "__main__":
    main()
