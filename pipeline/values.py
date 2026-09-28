"""
Messy Lab, step 2: put every value on one scale.

For each study, only the row chosen to represent each gene (std_feature_map.chosen_for_gene) is used.
  S1, S3  arrays deposited un-logged (linear)   -> value = log2(raw + 1)     value_type log2_array_signal
  S2      array already log2 (RMA)              -> value = raw               value_type log2_array_signal
  S4, S5  RNA-seq raw read counts               -> value = log2(CPM + 1)     value_type log2_cpm
          CPM = count / (all reads counted for that sample) x 1,000,000, recomputed here from raw counts.
          S4's own CPM columns are ignored: they were normalized together with 8 patient samples.
The +1 keeps zeros defined (log2 of 0 is undefined). The raw number is kept next to every value.

Values are comparable within a study (sample vs sample), not across studies (array vs RNA-seq).

Output: std_measurement in db/staging.duckdb.
Run from the Messy_Lab folder, after loaders.py, standardize.py and genes.py:   python3 pipeline/values.py
"""
import os
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parent.parent
DB_DIR = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db"))
STAGING = DB_DIR / "staging.duckdb"

SCALE = {  # study: (how to get value from raw_value, value_type, raw_value_type)
    "S1": ("log2(v.raw_value + 1)", "log2_array_signal", "linear_array_signal"),
    "S2": ("v.raw_value", "log2_array_signal", "log2_array_signal"),
    "S3": ("log2(v.raw_value + 1)", "log2_array_signal", "linear_array_signal"),
    "S4": ("log2(v.raw_value / lib.total * 1e6 + 1)", "log2_cpm", "raw_count"),
    "S5": ("log2(v.raw_value / lib.total * 1e6 + 1)", "log2_cpm", "raw_count"),
}


def main():
    con = duckdb.connect(str(STAGING))
    con.execute("DROP TABLE IF EXISTS std_measurement")
    # library size = every read counted for the sample, including rows that don't map to a gene
    con.execute("CREATE TEMP TABLE lib AS SELECT sample_id, sum(raw_value) AS total FROM stg_value "
                "WHERE study_id IN ('S4','S5') GROUP BY 1")
    first = True
    for study, (expr, vtype, rtype) in SCALE.items():
        sql = f"""
            SELECT v.sample_id, f.hgnc_id, {expr} AS value, '{vtype}' AS value_type,
                   v.raw_value, '{rtype}' AS raw_value_type, v.feature_id, v.source_file_id
            FROM stg_value v
            JOIN std_feature_map f ON f.study_id = v.study_id AND f.feature_id = v.feature_id AND f.chosen_for_gene
            LEFT JOIN lib ON lib.sample_id = v.sample_id
            WHERE v.study_id = '{study}'"""
        con.execute(("CREATE TABLE std_measurement AS " if first else "INSERT INTO std_measurement ") + sql)
        first = False

    # ---- checks ---------------------------------------------------------------------------------
    print(f"{'study':6} {'values':>10} {'samples x genes':>16} {'min':>7} {'median':>7} {'max':>7}  scale")
    for study, (_, vtype, _) in SCALE.items():
        n, lo, med, hi = con.execute("""select count(*), min(value), median(value), max(value)
                                        from std_measurement m join std_sample s using(sample_id) where s.study_id=?""", [study]).fetchone()
        ns = con.execute("select count(*) from std_sample where study_id=?", [study]).fetchone()[0]
        ng = con.execute("select count(*) from std_feature_map where study_id=? and chosen_for_gene", [study]).fetchone()[0]
        assert n == ns * ng, f"{study}: {n} values, expected {ns} x {ng}"
        print(f"{study:6} {n:>10,} {f'{ns} x {ng:,}':>16} {lo:>7.2f} {med:>7.2f} {hi:>7.2f}  {vtype}")
    bad = con.execute("select count(*) from std_measurement where value is null or isnan(value) or isinf(value)").fetchone()[0]
    assert bad == 0, f"{bad} values are missing, NaN or infinite"
    cpm = con.execute("""select min(t), max(t) from (select sample_id, sum(pow(2, value) - 1) t from std_measurement
                         where value_type='log2_cpm' group by 1)""").fetchone()
    print(f"\nNo missing or infinite values. RNA-seq CPM per sample sums to {cpm[0]:,.0f}-{cpm[1]:,.0f} "
          f"(under 1,000,000 because unmapped rows are left out).")

    # biology sanity check: estradiol should switch on GREB1, a classic estrogen target, in MCF-7
    print("\nSanity check, GREB1 in MCF-7 (mean log2 value, estradiol vs its control):")
    for study, ctrl_role in (("S1", "vehicle_control"), ("S5", "vehicle_control")):
        r = con.execute(f"""
            with g as (select m.sample_id, m.value from std_measurement m join std_gene g using(hgnc_id) where g.symbol='GREB1'),
                 s as (select s.sample_id, s.role, s.genetic_change,
                              bool_or(t.compound_id='E2') as e2 from std_sample s left join std_sample_treatment t using(sample_id)
                       where s.study_id='{study}' and s.cell_line_id='CVCL_0031' and s.genetic_change in ('none','empty_vector')
                       group by all)
            select avg(value) filter (where e2), avg(value) filter (where role='{ctrl_role}') from g join s using(sample_id)""").fetchone()
        print(f"   {study}: E2 {r[0]:.2f} vs control {r[1]:.2f}  ->  {r[0]-r[1]:+.2f} log2 ({2**(r[0]-r[1]):.1f}-fold)")
    con.close()


if __name__ == "__main__":
    main()
