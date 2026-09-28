"""
Messy Lab, step 3: validation checks.

Each check looks for one kind of problem in the loaded data and writes what it finds to the issue
table, tagged with the mess-inventory row it covers (M01..M75). Checks look at the data itself, so
they would find the same kinds of problem in a new study, not just the ones already listed.

Every issue says what was observed, what was expected, and how it was resolved:
  fixed_by_rule  a mapping rule already corrects it in the clean tables
  flagged        kept as is and marked, because the sources can't settle it
  excluded       the record is left out of analysis

Output: std_issue in db/staging.duckdb.
Run from the Messy_Lab folder, after the step 2 scripts:   python3 pipeline/checks.py
"""
import csv
import os
import re
from pathlib import Path

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DB_DIR = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db"))
STAGING = DB_DIR / "staging.duckdb"
REF = ROOT / "pipeline" / "reference"
IMPACT = {r["id"]: r["impact"].lower() for r in csv.DictReader(open(ROOT / "mess_inventory" / "mess_inventory.csv"))}
UNIT_TO_NM = {"M": 1e9, "mM": 1e6, "uM": 1e3, "µM": 1e3, "μM": 1e3, "nM": 1.0, "pM": 1e-3}
DATE_LIKE = r"^\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$"
ISSUES = []


def issue(check, inv, study, etype, eid, field, observed, expected, resolution, message):
    ISSUES.append(dict(check_name=check, inventory_id=inv, impact=IMPACT.get(inv, "minor") if inv else "minor",
                       study_id=study, entity_type=etype, entity_id=eid, field=field,
                       observed=None if observed is None else str(observed)[:500],
                       expected=None if expected is None else str(expected)[:500],
                       resolution=resolution, message=message))


def fields_of(con):
    long = con.execute("select study_id, sample_id, field, value from stg_sample_field").df()
    return {(s, g): dict(zip(d.field, d.value)) for (s, g), d in long.groupby(["study_id", "sample_id"])}


def compound_patterns():
    """Per study: compound_id -> regex of the spellings that study uses in its GEO fields."""
    pats = {}
    for r in csv.DictReader(open(REF / "compound_alias.csv")):
        if "paper only" in r["seen_in_field"] or "lab codes" in r["seen_in_field"]:
            continue
        flags = "" if len(r["alias"]) < 3 else "(?i)"
        pats.setdefault(r["seen_in_study"], {}).setdefault(r["compound_id"], []).append(
            flags + r"(?<![\w-])" + re.escape(r["alias"]) + r"(?![\w-])")
    return pats


def compounds_in(text, study_pats):
    return {cid for cid, ps in study_pats.items() if any(re.search(p, text or "") for p in ps)}


# ------------------------------------------------------------------------------------------------
def check_dose_above_100uM(con, F, pats):
    """1. Any dose of a study compound, written anywhere in the GEO text, that works out above 100 µM.
    A number counts as a dose only when one of the study's compound names follows within a few words,
    so '2 mM L-glutamine' or '30M reads' are not mistaken for doses."""
    texts = con.execute("select study_id, 'series' as sample_id, field, value from stg_series "
                        "union all select study_id, sample_id, field, value from stg_sample_field").df()
    seen = set()
    for r in texts.itertuples():
        text = r.value or ""
        for m in re.finditer(r"(?<![\d.])(\d+(?:\.\d+)?)\s?(M|mM|uM|µM|μM|nM|pM)\b", text):
            num, unit = m.groups()
            after = text[m.end(): m.end() + 40]
            if not compounds_in(after, pats.get(r.study_id, {})):
                continue
            nM = float(num) * UNIT_TO_NM[unit]
            key = (r.study_id, r.field, num + unit)
            if nM > 1e5 and key not in seen:
                seen.add(key)
                issue("dose_above_100uM", "M13", r.study_id, "study", r.study_id, r.field, f"{num} {unit}{after[:25]}",
                      "at most 100 uM (these ligands are used at nM to low uM)", "fixed_by_rule",
                      f"'{num} {unit}' in {r.field} equals {nM:,.0f} nM, a million-fold error from a lost micro sign; "
                      f"the sample fields say 1uM and the clean table stores 1,000 nM.")
    bad = con.execute("select sample_id, compound_id, dose_nM from std_sample_treatment where dose_nM > 100000").fetchall()
    for sid, cid, d in bad:
        issue("dose_above_100uM", "M13", None, "sample_treatment", f"{sid}|{cid}", "dose_nM", d, "<= 100000",
              "flagged", "A stored dose is above 100 uM.")


def check_fields_contradict(con, F, pats):
    """2. The same sample described differently by its own fields (title vs source vs treatment)."""
    candidates = {"title": "title", "source_name_ch1": "source name", "counts_file_treatment": "counts file label"}
    by_study = {}
    for (s, g), f in F.items():
        by_study.setdefault(s, []).append((g, f))
    for s, rows in by_study.items():
        sp = pats.get(s, {})
        treat_fields = [k for k in rows[0][1] if k.startswith("characteristics") or k in candidates]
        # only compare fields that this study uses to describe treatments at all
        used = [k for k in treat_fields
                if len({frozenset(compounds_in(f.get(k, ""), sp)) for g, f in rows}) > 1]   # varies between samples
        if len(used) < 2:
            continue
        for g, f in rows:
            sets = {k: compounds_in(f.get(k, ""), sp) for k in used}
            chars = set().union(*[v for k, v in sets.items() if k.startswith("characteristics")])
            views = {k: v for k, v in sets.items() if not k.startswith("characteristics")}
            views["characteristics"] = chars
            distinct = {frozenset(v) for v in views.values()}
            if len(distinct) > 1:
                odd = [k for k, v in views.items() if list(views.values()).count(v) == 1 and len(views) > 2]
                issue("fields_contradict", "M44", s, "sample", g, ",".join(odd) or ",".join(views),
                      "; ".join(f"{k}: {f.get(k, '') if k != 'characteristics' else 'treatment field'} -> {sorted(v) or 'nothing'}" for k, v in views.items()),
                      "every field describes the same treatment", "fixed_by_rule",
                      f"Fields disagree about what {g} received. The clean table follows the fields that agree "
                      f"(title, treatment field, lab code); the odd one out is {', '.join(odd) or 'unclear'}.")


def check_duplicate_titles(con, F):
    """3. Titles that become identical once the lab's own code is removed."""
    rows = []
    for (s, g), f in F.items():
        core = re.sub(r"\s*\((MD\d+)\)\s*$", "", f["title"]).strip()
        rows.append((s, g, f["title"], core))
    d = pd.DataFrame(rows, columns=["study", "gsm", "title", "core"])
    batch = dict(con.execute("select sample_id, batch from std_sample").fetchall())
    for (s, core), grp in d.groupby(["study", "core"]):
        if len(grp) < 2:
            continue
        batches = {batch.get(g) for g in grp.gsm}
        same_batch = len(batches) == 1
        for r in grp.itertuples():
            others = ", ".join(x for x in grp.gsm if x != r.gsm)
            if same_batch:
                issue("duplicate_titles", "M46", s, "sample", r.gsm, "title", r.title, "a unique title per sample",
                      "fixed_by_rule", f"Same title as {others} in the same experiment; the replicate number is repeated. "
                      "The clean table numbers them 1 and 2.")
            else:
                issue("duplicate_titles", "M47", s, "sample", r.gsm, "title", r.title, "a unique title per sample",
                      "flagged", f"Same title as {others}, but they belong to different experiments "
                      f"({', '.join(sorted(b for b in batches if b))}). Compare each only with controls from its own batch.")


def check_control_in_batch(con):
    """4. Treated samples with no matching control in the same batch."""
    rows = con.execute("""
        select t.sample_id, t.study_id, t.batch, t.cell_line_id, t.derived_line_id, t.genetic_change, t.genetic_change_target
        from std_sample t
        where t.role = 'treated' and not exists (
            select 1 from std_sample c
            where c.study_id = t.study_id and c.role in ('vehicle_control','untreated_control','parental')
              and coalesce(c.batch,'') = coalesce(t.batch,'')
              and c.cell_line_id = t.cell_line_id
              and (coalesce(c.derived_line_id,'') = coalesce(t.derived_line_id,'') or c.role = 'parental')
              and c.genetic_change is not distinct from t.genetic_change
              and c.genetic_change_target is not distinct from t.genetic_change_target)""").fetchall()
    for sid, s, b, cl, dl, gc, gt in rows:
        issue("no_control_in_batch", "M50" if s == "S5" else None, s, "sample", sid, "batch", b,
              "a vehicle, untreated or parental control with the same cells in the same batch", "flagged",
              f"No matching control in batch {b or '(none)'} for these cells ({dl or cl}). "
              "Any comparison has to borrow a control from another run, which mixes run-to-run differences into the result.")


def check_value_scale(con):
    """5. Values whose range doesn't match what the processing note says they are."""
    for s in ["S1", "S2", "S3", "S4", "S5"]:
        note = " ".join(v for (v,) in con.execute(
            "select distinct value from stg_sample_field where study_id=? and field like 'data_processing%'", [s]).fetchall())
        lo, hi, nonint = con.execute("select min(raw_value), max(raw_value), count(*) filter (where raw_value <> round(raw_value)) "
                                     "from stg_value where study_id=?", [s]).fetchone()
        looks_log = hi <= 25 and lo >= -5
        says_log = bool(re.search(r"\bG?C?RMA\b", note))   # methods whose output is log2
        says_normalized = bool(re.search(r"RPKM|FPKM|TPM|CPM", note))
        if says_log and not looks_log:
            issue("value_scale_mismatch", "M55", s, "study", s, "values", f"range {lo:,.2f} to {hi:,.2f}",
                  "about 0 to 20 for log2 values", "fixed_by_rule",
                  f"Processing note says '{re.search(r'G?C?RMA', note).group(0)}' (a log2 method) but values reach {hi:,.0f}, "
                  "so they were un-logged. The clean table applies log2(x+1).")
        if says_normalized and nonint == 0 and s in ("S5",):
            issue("value_scale_mismatch", "M59", s, "study", s, "values", "all whole numbers (raw read counts)",
                  "fractional normalized values, as the processing note describes", "fixed_by_rule",
                  "Processing note describes normalized values (nRPKM), but the file holds raw counts. "
                  "The clean table recomputes log2 CPM from the counts.")
        if s == "S3":
            issue("value_scale_mismatch", "M57", s, "study", s, "values", "GCOS/MAS5-style linear values, 9 zeros",
                  "RMA + quantile normalized values, as the paper describes", "flagged",
                  "The deposited values are not the values the paper analyzed (paper: RMA, quantile, GFP correction). "
                  "The clean table uses the deposited values, log2(x+1).")


def check_date_gene_names(con):
    """6. Gene names that look like dates (Excel damage)."""
    rows = con.execute("select study_id, feature_id, symbol_as_given, entrez_as_given from stg_feature").df()
    hit = rows[rows.symbol_as_given.fillna("").str.match(DATE_LIKE)]
    mapped = dict(con.execute("select study_id || '|' || feature_id, hgnc_id from std_feature_map").fetchall())
    sym = dict(con.execute("select hgnc_id, symbol from std_gene").fetchall())
    for r in hit.itertuples():
        hid = mapped.get(f"{r.study_id}|{r.feature_id}")
        issue("date_like_gene_name", "M65", r.study_id, "feature", r.feature_id, "symbol", r.symbol_as_given,
              "a gene symbol", "fixed_by_rule",
              f"'{r.symbol_as_given}' is a gene name Excel turned into a date. Matched through Entrez {r.entrez_as_given} "
              f"to {sym.get(hid, 'nothing')} ({hid}).")


def check_duplicate_gsm(con):
    """7. The same sample arriving twice, or a series that also sits inside another series."""
    dup = con.execute("""select sample_id, count(distinct study_id) n, string_agg(distinct study_id, ',') studies
                         from stg_sample_field where field='geo_accession' group by 1 having count(*) > 1""").fetchall()
    for sid, n, studies in dup:
        issue("duplicate_sample_id", "M53", None, "sample", sid, "geo_accession", studies, "each GSM once", "flagged",
              f"{sid} arrives {n} times.")
    for s, v in con.execute("select study_id, value from stg_series where field='Series_relation' and value like 'SubSeries of%'").fetchall():
        issue("duplicate_sample_id", "M53", s, "study", s, "Series_relation", v, "samples reachable through one accession",
              "fixed_by_rule", f"This series is also part of {v.split(':')[-1].strip()}. Loading both would double count "
              "its samples, so only the SubSeries is loaded and GSM IDs are kept unique.")
    for s, n in con.execute("""select study_id, count(*) from stg_series where field='Series_relation'
                               and value like 'Reanalysis of%' group by 1""").fetchall():
        issue("out_of_scope_files", "M52", s, "study", s, "Series_relation", f"{n} reanalysed samples from another study",
              "only this study's cell-line samples", "excluded",
              f"The series also carries {n} patient tumor samples from GSE58708. They are not loaded.")


def check_required_and_unknown(con):
    """8. Required fields still blank, and fields the sources leave unknown."""
    req = ["cell_line_id", "genetic_change", "role", "resistance_phenotype", "hormone_deprived", "in_drug_at_harvest", "source_file_id"]
    for c in req:
        for (sid,) in con.execute(f"select sample_id from std_sample where {c} is null or {c} = ''").fetchall():
            issue("required_field_blank", None, None, "sample", sid, c, None, "a value", "flagged", f"{c} is blank.")
    for sid, s in con.execute("""select s.sample_id, s.study_id from std_sample s join stg_sample_field f using(sample_id)
                                 where f.field='characteristics_ch1#9' and f.value=''""").fetchall():
        issue("unknown_value", "M24", s, "sample", sid, "exposure time", "(blank)", "a time in hours", "flagged",
              "GEO gives no time for this control.")
    for sid, cid, s in con.execute("""select t.sample_id, t.compound_id, s.study_id from std_sample_treatment t join std_sample s using(sample_id)
                                      where t.exposure_type='acute' and t.exposure_h is null""").fetchall():
        issue("unknown_value", "M25", s, "sample_treatment", f"{sid}|{cid}", "exposure_h", None, "hours of exposure",
              "flagged", "Exposure length isn't stated anywhere; the paper only says RNA was taken 72 h after transfection.")
    for sid, cid, s, st, given in con.execute("""select t.sample_id, t.compound_id, s.study_id, t.dose_status, t.dose_as_given
                                                 from std_sample_treatment t join std_sample s using(sample_id)
                                                 where t.dose_status in ('unknown','conflicting')""").fetchall():
        inv = "M11" if st == "conflicting" else "M03"
        issue("dose_" + st, inv, s, "sample_treatment", f"{sid}|{cid}", "dose_nM", given, "one stated dose", "flagged",
              "Dose sources disagree; 10 nM (gene-expression methods) is stored." if st == "conflicting"
              else "Dose not stated in GEO or the paper.")
    for sid, s in con.execute("""select distinct t.sample_id, s.study_id from std_sample_treatment t join std_sample s using(sample_id)
                                 where t.compound_id='TAM_UNSPECIFIED'""").fetchall():
        issue("compound_form_unstated", "M03", s, "sample", sid, "compound", "tamoxifen", "4-hydroxytamoxifen or tamoxifen citrate, stated",
              "flagged", "The study never says which form of tamoxifen; kept separate from 4-hydroxytamoxifen.")
    for sid, s in con.execute("select sample_id, study_id from std_sample where in_drug_at_harvest='unknown'").fetchall():
        issue("unknown_value", "M36", s, "sample", sid, "in_drug_at_harvest", "unknown", "yes or no", "flagged",
              "Resistant cells were kept in tamoxifen; the paper doesn't say whether it was removed before this sample.")


def check_inferred_and_assigned(con):
    """Batches worked out from the data, and replicate numbers the lab never gave."""
    for sid, s, b in con.execute("select sample_id, study_id, batch from std_sample where batch_basis='inferred'").fetchall():
        issue("batch_inferred", "M47" if s == "S3" else "M50", s, "sample", sid, "batch", b, "a batch stated by the lab",
              "flagged", f"Batch {b} is inferred from the data (sample IDs, read depth or expression), not stated by the lab.")
    for s, n in con.execute("select study_id, count(*) from std_sample where replicate_basis='assigned' group by 1").fetchall():
        issue("replicate_assigned", "M48" if s == "S5" else ("M51" if s == "S4" else "M45"), s, "study", s, "replicate",
              f"{n} samples without a replicate number", "replicate numbers from the lab", "fixed_by_rule",
              f"{n} samples had no replicate number; numbered 1, 2, 3 within identical conditions in GSM order.")
    for s, n in con.execute("""select study_id, count(*) from std_sample s join stg_sample_field f using(sample_id, study_id)
                               where f.field='title' and regexp_matches(f.value, 'replicate\\d b$') group by 1""").fetchall():
        issue("unexplained_label", "M42", s, "study", s, "title", f"{n} titles end in ' b'", "an explained label", "flagged",
              "The 'b' suffix is never explained.")


def check_genes(con):
    """Gene rows that can't be matched cleanly, and old names."""
    for s, status, n in con.execute("select study_id, map_status, count(*) from std_feature_map "
                                    "where map_status <> 'mapped' group by all order by 1").fetchall():
        inv = {"multiple_genes": "M61", "corrupt_annotation": "M63", "control_probe": "M61"}.get(status,
              {"S1": "M61", "S3": "M61", "S2": "M62", "S4": "M64", "S5": "M66"}[s])
        if status == "corrupt_annotation":
            for (fid,) in con.execute("select feature_id from std_feature_map where study_id=? and map_status='corrupt_annotation'", [s]).fetchall():
                issue("gene_annotation_corrupt", "M63", s, "feature", fid, "description", "only '#' characters",
                      "a gene description", "excluded", "Annotation is unreadable; this row carries no gene.")
            continue
        issue("gene_rows_" + status, inv, s, "study", s, "feature_id", f"{n:,} rows", "one gene per row",
              "excluded", f"{n:,} {s} rows are {status.replace('_', ' ')}; they carry no gene into the clean tables.")
    for s, n in con.execute("""select study_id, count(*) from std_feature_map f join std_gene g using(hgnc_id)
                               where f.symbol_as_given is not null and f.symbol_as_given <> g.symbol
                               and not regexp_matches(f.symbol_as_given, '^\\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$')
                               group by 1 order by 1""").fetchall():
        issue("gene_name_outdated", "M67", s, "study", s, "symbol", f"{n:,} rows use a name that differs from today's official symbol",
              "current HGNC symbols", "fixed_by_rule", f"{n:,} {s} gene names differ from today's HGNC symbol "
              "(renamed genes and aliases). Joined on HGNC ID, so they still match across studies.")
    for s, n in con.execute("""select study_id, count(*) - count(distinct hgnc_id) from std_feature_map
                               where map_status='mapped' group by 1 having count(*) > count(distinct hgnc_id) order by 1""").fetchall():
        issue("several_rows_per_gene", "M61" if s in ("S1", "S3") else ("M62" if s == "S2" else "M64"), s, "study", s,
              "feature_id", f"{n:,} extra rows for genes already measured", "one row per gene", "fixed_by_rule",
              "Genes with several probes or rows keep the one with the highest average expression.")


def check_labels_and_layout(con, F, pats):
    """Same compound spelled differently, S5 column order, per-study naming."""
    for s, sp in pats.items():
        for cid, ps in sp.items():
            if len(ps) > 1:
                inv = {"OHT": {"S1": "M01", "S3": "M04", "S5": "M06"}.get(s), "GNE274": "M07", "E2": "M09",
                       "ETHANOL": "M16", "DMSO": "M21", "TAM_UNSPECIFIED": "M03"}.get(cid)
                inv = inv if isinstance(inv, str) or inv is None else inv
                spell = [re.sub(r"^\(\?i\)|\(\?<!\[\\w-\]\)|\(\?!\[\\w-\]\)|\\", "", p) for p in ps]
                issue("compound_spellings", inv, s, "study", s, "compound", " / ".join(spell), "one name per compound",
                      "fixed_by_rule", f"{cid} is written {len(ps)} ways in {s}; all map to one compound.")
    moved = con.execute("""select count(*) from stg_sample_field p join (
                             select sample_id, row_number() over (order by sample_id) + 11 as pos from
                             (select distinct sample_id from stg_sample_field where study_id='S5')) m using(sample_id)
                           where p.field='counts_file_column_position' and cast(p.value as int) <> m.pos""").fetchone()[0]
    if moved:
        issue("column_order", "M70", "S5", "study", "S5", "counts file columns", f"{moved} of 82 columns out of GSM order",
              "columns in the same order as the sample sheet", "fixed_by_rule",
              "The counts file lists samples in a different order from the series matrix; columns are matched by sample ID.")
    for s, v in con.execute("select study_id, value from stg_sample_field where field='treatment_protocol_ch1' "
                            "and value='see above' group by all").fetchall():
        issue("empty_reference", "M71", s, "study", s, "treatment_protocol", v, "a protocol", "flagged",
              "The treatment protocol just says 'see above'; details were taken from the growth protocol and paper.")
    s1 = [f for (st, g), f in F.items() if st == "S1"]
    mois = {m for f in s1 for m in re.findall(r"moi (?:of )?(\d+(?: or \d+)?)|moi(\d+)", " ".join(f.values()))}
    if len({x for pair in mois for x in pair if x}) > 1:
        issue("fields_contradict", "M41", "S1", "study", "S1", "virus dose (moi)",
              " vs ".join(sorted({x for pair in mois for x in pair if x})), "one moi", "flagged",
              "The sample protocol says moi 5 or 50; characteristics and paper say 10. 10 is recorded.")


def check_design_balance(con):
    """Groups with fewer replicates than the study's usual number, and compounds tested in only some lines."""
    groups = con.execute("""
        with sig as (select sample_id, string_agg(compound_id || '@' || coalesce(cast(dose_nM as varchar),'') || '@' ||
                            coalesce(cast(exposure_h as varchar), exposure_type), '+' order by compound_id) k
                     from sample_treatment_src group by 1)
        select s.study_id, coalesce(s.batch,'') batch, s.cell_line_id, coalesce(s.derived_line_id,'') dl, s.genetic_change,
               coalesce(s.genetic_change_target,'') tgt, s.role, coalesce(sig.k,'none') k, count(*) n, list(s.sample_id order by s.sample_id) ids
        from sample_src s left join sig using(sample_id) group by all""".replace("sample_treatment_src", "std_sample_treatment")
                          .replace("sample_src", "std_sample")).df()
    for study, g in groups.groupby("study_id"):
        usual = int(g.n.mode().max())
        for r in g[g.n < usual].itertuples():
            for sid in r.ids:
                issue("fewer_replicates", "M43" if study == "S1" else None, study, "sample", sid, "replicate",
                      f"{r.n} sample(s) in this condition", f"{usual}, the usual number in {study}", "flagged",
                      f"This condition has {r.n} replicate(s) where {study} usually has {usual}; its average is less certain.")
    gaps = con.execute("""
        with used as (select distinct s.study_id, s.cell_line_id, t.compound_id from std_sample s join std_sample_treatment t using(sample_id)),
             lines as (select distinct study_id, cell_line_id from std_sample where derived_line_id is null),
             cmpds as (select distinct study_id, compound_id from used)
        select l.study_id, c.compound_id, list(l.cell_line_id order by l.cell_line_id) missing
        from lines l join cmpds c using(study_id)
        where not exists (select 1 from used u where u.study_id=l.study_id and u.cell_line_id=l.cell_line_id and u.compound_id=c.compound_id)
        group by all""").fetchall()
    for study, cid, missing in gaps:
        issue("design_gap", "M49" if study == "S5" else None, study, "study", study, "compound x cell line",
              f"{cid} not tested in {', '.join(missing)}", f"{cid} in every line the study used", "flagged",
              f"{cid} was only tested in some of {study}'s cell lines; comparisons across lines can't include it for "
              f"{', '.join(missing)}.")


def main():
    con = duckdb.connect(str(STAGING))
    F = fields_of(con)
    pats = compound_patterns()
    check_dose_above_100uM(con, F, pats)
    check_fields_contradict(con, F, pats)
    check_duplicate_titles(con, F)
    check_control_in_batch(con)
    check_value_scale(con)
    check_date_gene_names(con)
    check_duplicate_gsm(con)
    check_required_and_unknown(con)
    check_inferred_and_assigned(con)
    check_genes(con)
    check_labels_and_layout(con, F, pats)
    check_design_balance(con)

    I = pd.DataFrame(ISSUES)
    con.register("df", I)
    con.execute("CREATE OR REPLACE TABLE std_issue AS SELECT * FROM df")
    con.close()

    chk = duckdb.connect(":memory:")
    chk.execute((ROOT / "pipeline/schema.sql").read_text())
    chk.execute(f"ATTACH '{STAGING}' AS stg (READ_ONLY)")
    chk.execute("INSERT INTO study (study_id, geo_accession, geo_title, technology, platform, raw_value_scale, samples_in_geo) "
                "VALUES ('S1','S1','p','microarray','x','log2',0),('S2','S2','p','microarray','x','log2',0),"
                "('S3','S3','p','microarray','x','log2',0),('S4','S4','p','microarray','x','log2',0),('S5','S5','p','microarray','x','log2',0)")
    chk.execute("INSERT INTO issue (check_name, inventory_id, impact, study_id, entity_type, entity_id, field, observed, expected, "
                "resolution, message) SELECT check_name, inventory_id, impact, study_id, entity_type, entity_id, field, observed, "
                "expected, resolution, message FROM stg.std_issue")
    print(f"{len(I):,} issues written; all pass the schema's rules.\n")
    summary = I.groupby(["check_name", "resolution"]).agg(n=("entity_id", "size"), studies=("study_id", lambda x: ",".join(sorted(set(x.dropna())))),
                                                          inventory=("inventory_id", lambda x: ",".join(sorted(set(x.dropna())))))
    print(summary.to_string())
    covered = sorted(set(I.inventory_id.dropna()))
    print(f"\nInventory rows covered by a check: {len(covered)} of 75")


if __name__ == "__main__":
    main()
