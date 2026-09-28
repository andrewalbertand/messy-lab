"""
Messy Lab, step 2: standardize sample labels.

Reads each study's raw sample fields from db/staging.duckdb (written by loaders.py) and applies that
study's mapping file, pipeline/mappings/S1.yaml .. S5.yaml, to produce clean `sample` and
`sample_treatment` rows plus a `provenance` row for every field it fills in.

How a mapping file works
  defaults:  values every sample in the study gets, with their source.
  rules:     applied in order to every sample. `when` holds regular expressions tested against the
             sample's raw GEO fields; all must match. Named groups, e.g. (?P<h>\\d+), can be used in
             `set` values as {h}. A later rule overrides an earlier one for the same field.
             `treatments` adds one row per compound. `cell_line_alias`, `derived_line_alias` and
             `compound_alias` are looked up in pipeline/reference/, so a mapping never invents an ID.
  samples:   fixes for named samples (GSM IDs) that the general rules can't express.
  Every rule names its source (file + location) and, where it fixes a known problem, the
  mess-inventory row (M01..M75).

Nothing is guessed: a field no rule fills stays empty, and required fields left empty stop the run.

Run from the Messy_Lab folder (after loaders.py):   python3 pipeline/standardize.py
"""
import csv
import os
import re
from pathlib import Path

import duckdb
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
DB_DIR = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db"))
STAGING = DB_DIR / "staging.duckdb"
MAPS = ROOT / "pipeline" / "mappings"
REF = ROOT / "pipeline" / "reference"
STUDIES = ["S1", "S2", "S3", "S4", "S5"]

SAMPLE_COLS = ["sample_id", "study_id", "geo_title", "lab_sample_code", "cell_line_id", "derived_line_id",
               "genetic_change", "genetic_change_target", "genetic_change_method", "role", "resistance_phenotype",
               "batch", "batch_basis", "replicate", "replicate_basis", "hormone_deprived", "in_drug_at_harvest",
               "culture_medium", "harvest_note", "include_in_analysis", "exclusion_reason", "source_file_id", "notes"]
REQUIRED = ["cell_line_id", "genetic_change", "role", "resistance_phenotype", "hormone_deprived", "in_drug_at_harvest"]
TREAT_COLS = ["sample_id", "compound_id", "dose_nM", "dose_status", "exposure_type", "exposure_h", "exposure_note",
              "label_as_given", "dose_as_given", "source_file_id"]
NUMERIC = {"replicate": int, "dose_nM": float, "exposure_h": float}


def aliases(name):
    table = {}
    for r in csv.DictReader(open(REF / name)):
        table[(r["alias"], r["seen_in_study"])] = r["cell_line_id" if "cell_line" in name else "compound_id"]
    return table


CELL_ALIAS, COMPOUND_ALIAS = aliases("cell_line_alias.csv"), aliases("compound_alias.csv")


def fill(template, ctx):
    """Substitute {group} placeholders; return numbers as numbers where the field needs one."""
    if not isinstance(template, str):
        return template
    return template.format(**ctx) if "{" in template else template


def match(when, fields):
    ctx = {}
    for field, pattern in (when or {}).items():
        m = re.search(pattern, fields.get(field, ""))
        if not m:
            return None
        ctx.update({k: v for k, v in m.groupdict().items() if v is not None})
    return ctx


def standardize(study, fields_by_sample, series_titles):
    spec = yaml.safe_load(open(MAPS / f"{study}.yaml"))
    samples, treatments, prov = [], [], []

    def record(entity_type, entity_id, field, value, source, raw, rule):
        prov.append((entity_type, entity_id, field, None if value is None else str(value),
                     source["file"], source["location"], raw, rule))

    for gsm, f in fields_by_sample.items():
        row = {c: None for c in SAMPLE_COLS}
        row.update(sample_id=gsm, study_id=study, geo_title=f["title"], include_in_analysis=True,
                   source_file_id=spec["source_file_id"])
        notes, treats = [], {}

        def apply(block, ctx, source, raw, rule_name):
            for key, val in (block.get("set") or {}).items():
                val = fill(val, ctx)
                if key in ("cell_line_alias", "derived_line_alias"):
                    target = "cell_line_id" if key == "cell_line_alias" else "derived_line_id"
                    found = CELL_ALIAS.get((val, study))
                    assert found, f"{study} {gsm}: no cell line alias '{val}' for {study} in cell_line_alias.csv"
                    row[target] = found
                    record("sample", gsm, target, found, source, raw, f"{rule_name}: '{val}' -> {found}")
                    continue
                if key == "note":
                    notes.append(val)
                    continue
                if key in NUMERIC and val not in (None, ""):
                    val = NUMERIC[key](val)
                row[key] = val
                record("sample", gsm, key, val, source, raw, rule_name)
            if block.get("clear_treatments"):
                treats.clear()
            for t in block.get("treatments") or []:
                t = {k: fill(v, ctx) for k, v in t.items()}
                alias = t.pop("compound_alias")
                cid = COMPOUND_ALIAS.get((alias, study))
                assert cid, f"{study} {gsm}: no compound alias '{alias}' for {study} in compound_alias.csv"
                tr = {c: None for c in TREAT_COLS}
                tr.update(sample_id=gsm, compound_id=cid, source_file_id=spec["source_file_id"])
                for k, v in t.items():
                    tr[k] = NUMERIC[k](v) if k in NUMERIC and v not in (None, "") else v
                if tr["label_as_given"] is None:
                    tr["label_as_given"] = alias
                treats[cid] = tr
                for k in ("dose_nM", "dose_status", "exposure_type", "exposure_h"):
                    record("sample_treatment", f"{gsm}|{cid}", k, tr[k], source, raw,
                           f"{rule_name}: '{alias}' -> {cid}")

        # defaults
        for key, d in (spec.get("defaults") or {}).items():
            apply({"set": {key: d["value"]}}, {}, d["source"], d.get("raw"), "study default")
        # rules
        for rule in spec.get("rules") or []:
            ctx = match(rule.get("when"), f)
            if ctx is None:
                continue
            raw = "; ".join(f"{k}='{f.get(k, '')}'" for k in (rule.get("when") or {}))
            name = rule["name"] + (f" [{rule['inventory']}]" if rule.get("inventory") else "")
            apply(rule, ctx, rule["source"], raw, name)
        # named-sample fixes
        fix = (spec.get("samples") or {}).get(gsm)
        if fix:
            name = fix["name"] + (f" [{fix['inventory']}]" if fix.get("inventory") else "")
            apply(fix, {}, fix["source"], f"title='{f['title']}'", name)

        row["notes"] = " ".join(notes) or None
        missing = [c for c in REQUIRED if row[c] in (None, "")]
        assert not missing, f"{study} {gsm} ('{f['title']}'): no rule set {missing}"
        samples.append(row)
        treatments.extend(treats.values())

    # replicate numbers the lab didn't give: 1, 2, 3 within identical conditions, in GSM order
    sdf = pd.DataFrame(samples, columns=SAMPLE_COLS)
    tdf = pd.DataFrame(treatments, columns=TREAT_COLS)
    sig = (tdf.assign(k=tdf.compound_id + "@" + tdf.dose_nM.astype(str) + "@" + tdf.exposure_h.astype(str))
              .groupby("sample_id").k.apply(lambda s: "+".join(sorted(s))))
    sdf["_cond"] = (sdf[["cell_line_id", "derived_line_id", "genetic_change", "genetic_change_target", "role", "batch"]]
                    .astype(str).agg("|".join, axis=1) + "|" + sdf.sample_id.map(sig).fillna("none"))
    need = sdf.replicate.isna()
    sdf.loc[need, "replicate"] = (sdf[need].sort_values("sample_id").groupby("_cond").cumcount() + 1)
    sdf.loc[need, "replicate_basis"] = "assigned"
    for gsm, rep in sdf.loc[need, ["sample_id", "replicate"]].itertuples(index=False):
        prov.append(("sample", gsm, "replicate", str(int(rep)), spec["source_file_id"], "sample order (GSM)", None,
                     "assigned: numbered within identical conditions because the lab gave no replicate number"))
    sdf = sdf.drop(columns="_cond")
    sdf["replicate"] = sdf.replicate.astype("Int64")
    return sdf, tdf, prov


def main():
    con = duckdb.connect(str(STAGING))
    all_s, all_t, all_p = [], [], []
    for study in STUDIES:
        long = con.execute("select sample_id, field, value from stg_sample_field where study_id=?", [study]).df()
        by_sample = {g: dict(zip(d.field, d.value)) for g, d in long.groupby("sample_id", sort=True)}
        s, t, p = standardize(study, by_sample, None)
        all_s.append(s); all_t.append(t); all_p += p
        print(f"{study}: {len(s):>3} samples, {len(t):>3} treatment rows | roles {s.role.value_counts().to_dict()}")
    S = pd.concat(all_s, ignore_index=True)
    T = pd.concat([t for t in all_t if len(t)], ignore_index=True)
    P = pd.DataFrame(all_p, columns=["entity_type", "entity_id", "field", "clean_value", "source_file_id",
                                     "source_location", "raw_text", "rule"])
    for name, df in (("std_sample", S), ("std_sample_treatment", T), ("std_provenance", P)):
        con.register("df", df)
        con.execute(f"CREATE OR REPLACE TABLE {name} AS SELECT * FROM df")
        con.unregister("df")
    con.close()

    # prove the rows satisfy the real schema: load them into a fresh copy with the reference lists
    chk = duckdb.connect(":memory:")
    chk.execute((ROOT / "pipeline" / "schema.sql").read_text())
    chk.execute(f"ATTACH '{STAGING}' AS stg (READ_ONLY)")
    chk.execute("INSERT INTO source_file SELECT source_file_id, path, study_id, file_kind, sha256, size_bytes, origin_url, NULL, NULL FROM stg.stg_source_file")
    for t, where in (("compound", ""), ("compound_alias", ""), ("cell_line", " WHERE NOT is_derived"),
                     ("cell_line", " WHERE is_derived"), ("cell_line_alias", "")):
        chk.execute(f"INSERT INTO {t} SELECT * FROM read_csv('{REF / (t + '.csv')}', header=true, nullstr=''){where}")
    targets = [x for x in S.genetic_change_target.dropna().unique()]
    hg = pd.read_csv(ROOT / "studies/data/reference/hgnc_complete_set_2026-09-26.txt", sep="\t", dtype=str, usecols=["hgnc_id", "symbol"])
    for hid in targets:
        chk.execute("INSERT INTO gene (hgnc_id, symbol, source_file_id) VALUES (?, ?, 'HGNC_2026-09-26')",
                    [hid, hg.loc[hg.hgnc_id == hid, "symbol"].iloc[0]])
    chk.execute("INSERT INTO study (study_id, geo_accession, geo_title, technology, platform, raw_value_scale, samples_in_geo) "
                "SELECT DISTINCT study_id, study_id, 'placeholder', 'microarray', 'x', 'log2', 0 FROM stg.std_sample")
    chk.execute("INSERT INTO sample SELECT * FROM stg.std_sample")
    chk.execute("INSERT INTO sample_treatment SELECT * FROM stg.std_sample_treatment")
    chk.execute("INSERT INTO provenance (entity_type, entity_id, field, clean_value, source_file_id, source_location, raw_text, rule) "
                "SELECT * FROM stg.std_provenance")
    print(f"\nAll {len(S)} samples and {len(T)} treatment rows pass the schema's rules. "
          f"{len(P):,} provenance rows written.")


if __name__ == "__main__":
    main()
