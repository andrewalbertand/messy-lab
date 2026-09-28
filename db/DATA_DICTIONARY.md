# Messy Lab · Data dictionary

This database (`db/messy_lab.duckdb`, DuckDB) holds five published gene-expression studies of
**tamoxifen response in breast cancer cells**, cleaned into one system. It has 265 samples and about
5.5 million expression values.

It was built by `python3 pipeline/build.py` from the raw GEO downloads in `studies/data/raw`. Nothing
was edited by hand, and every cleaned value can be traced back to the file and text it came from.

Read the **ground rules** first. They are what separate a right answer from a plausible wrong one.

---

## Ground rules

1. **Compare samples within a study, never raw values across studies.**
   - Arrays and RNA-seq are on different scales, and even two arrays aren't comparable in absolute terms.
   - To compare studies, compute a change within each study (treated vs its control), then compare the changes.
2. **Compare a treated sample only with controls from the same `batch`.**
   - A batch is one experimental run, and runs differ from each other.
   - Where a batch has no matching control, the `issue` table says so (check `no_control_in_batch`).
3. **"Unknown" means unknown.** Nothing was guessed.
   - When no source gives a value, it is NULL or `'unknown'`, and the reason is in `notes`, `dose_status` or the `issue` table.
   - Don't fill it in.
4. **4-hydroxytamoxifen and "tamoxifen, form not stated" are different compounds.**
   - `OHT` is 4-hydroxytamoxifen, the active form of the drug, used by studies S1, S3, S4 and S5.
   - `TAM_UNSPECIFIED` is study S2, which never says which form it used.
   - Pool them only if you say so explicitly.
5. **Genes are keyed on `hgnc_id`, not on names.**
   - Names change (WHSC1 is now NSD2), and one source had names damaged by Excel ("1-Mar").
   - Look genes up in `gene` (current, previous and alias symbols), then join on `hgnc_id`.
6. **At 0 hours nothing has been added yet.**
   - S2 labels its 0 h samples with the drug that will be added, for example "WT E2 0h".
   - Here they are `role = 'untreated_control'` with no treatment rows.
7. **Check the `issue` table before trusting a surprising result.** It lists every known problem with the source data, what was done about it, and which records it touches.

---

## The five studies

| study_id | GEO | Lab, year | What was done | Samples | Values are |
|---|---|---|---|---|---|
| S1 | GSE4025 | Frasor, Univ. of Illinois, 2006 | MCF-7 given extra estrogen receptor beta (ERβ) or an empty virus, then vehicle, E2 or 4-OHT for 24 h | 17 | array, log2 |
| S2 | GSE21618 | Oyama, Univ. of Tokyo and RIKEN, 2011 | Wild-type and tamoxifen-resistant MCF-7 given E2, heregulin and/or tamoxifen, from 0 to 48 h | 143 | array, log2 |
| S3 | GSE26298 | Salazar, Univ. of Toledo, 2011 | Hormone-starved MCF-7 with ESR1 or RARA knocked down, vehicle or 4-OHT | 12 | array, log2 |
| S4 | GSE111151 | Hultsch, Univ. of Helsinki, 2018 | Four parental lines and seven tamoxifen-resistant derivatives, one sample each | 11 | RNA-seq, log2 CPM |
| S5 | GSE117942 | Guan, Genentech, 2019 | Seven lines given DMSO, E2, 4-OHT or four receptor degraders for 24 h, in duplicate | 82 | RNA-seq, log2 CPM |

---

## How the tables fit together

```
study ─┬─< sample ─┬─< sample_treatment >── compound ──< compound_alias
       │           ├─< measurement >── gene ──< feature_map >── study
       │           └── cell_line (base line, and derived line) ──< cell_line_alias
       └─< issue            provenance  (every cleaned field → source file + raw text)
source_file  (every raw file used, with its checksum; every row above points to one)
```

`A ─< B` means one A has many B. Two views, `v_sample` and `v_expression`, join the common pieces for you.

Common abbreviations:

- **E2:** estradiol, the natural estrogen.
- **4-OHT:** 4-hydroxytamoxifen.
- **HRG:** heregulin, a growth factor.
- **ER:** estrogen receptor.
- **CPM:** counts per million sequencing reads.
- **GSM:** GEO's sample ID.

---

## sample

One row per sample. A sample is one dish of cells whose RNA was measured once. What was added to it
is in `sample_treatment`.

| Column | Meaning | Type / units | Allowed values | Comes from |
|---|---|---|---|---|
| `sample_id` | GEO sample ID, the key | text, e.g. `GSM645710` | unique | GEO |
| `study_id` | Which study | text | `S1`–`S5` | GEO series |
| `geo_title` | Sample name exactly as the lab wrote it | text | | GEO title |
| `lab_sample_code` | The lab's own ID | text, e.g. `MD31`, `TamR.Ctrl.00h`, `SAM24314537` | may be NULL | GEO description or title |
| `cell_line_id` | The base cell line, always the original line even for resistant derivatives | Cellosaurus ID, e.g. `CVCL_0031` = MCF-7 | see `cell_line` | GEO fields, mapped through `cell_line_alias` |
| `derived_line_id` | The lab-made resistant version, if any | e.g. `S2_TAMR_SELECTED`, `S4_MCF7_TAM1` | NULL for ordinary cells | titles, mapped through `cell_line_alias` |
| `genetic_change` | Whether the cells' genes were altered | text | `none`, `overexpression`, `empty_vector`, `knockdown`, `non_targeting_control` | titles and characteristics |
| `genetic_change_target` | The altered gene | HGNC ID, e.g. `HGNC:3468` = ESR2 | NULL when none | titles and characteristics |
| `genetic_change_method` | How it was altered | text, e.g. `adenovirus carrying ESR2, moi 10` | | protocols and paper |
| `role` | The sample's part in the experiment | text | `treated`, `vehicle_control` (solvent only), `untreated_control` (nothing added, including 0 h baselines), `parental` (original line in a resistant-vs-parental comparison) | characteristics, via the study's mapping file |
| `resistance_phenotype` | Whether the cells resist tamoxifen | text | `sensitive`, `resistant`, `unknown` (S5 doesn't say) | characteristics and papers |
| `batch` | The experimental run. Compare within a batch. | text, e.g. `S3_A`, `S2_wt3`, `S5_SAM24322` | NULL when the study is one run | lab codes, or inferred (see `batch_basis`) |
| `batch_basis` | How the batch is known | text | `stated` (the lab says so), `inferred` (worked out from IDs, read depth or expression), NULL | |
| `replicate` | Which copy among identically treated samples in the same batch | integer, 1, 2, 3 … | | titles or lab codes, else assigned |
| `replicate_basis` | How the replicate number is known | text | `stated`, `assigned` (numbered by us in GSM order), NULL | |
| `hormone_deprived` | Whether estrogen was removed from the medium before the experiment. This matters for estrogen-response genes. | text | `yes`, `no`, `unknown` | protocols |
| `in_drug_at_harvest` | Whether a drug or ligand was still present when the RNA was taken | text | `yes`, `no`, `unknown` | protocols and papers |
| `culture_medium` | What the cells were grown in | free text | | protocols and papers |
| `harvest_note` | When and how RNA was collected, where known | free text | | protocols and papers |
| `include_in_analysis` | Whether to use the sample | boolean | TRUE for all 265 | |
| `exclusion_reason` | Why it's excluded. Required when `include_in_analysis` is FALSE. | text | | |
| `source_file_id` | The raw file the row was built from | text | see `source_file` | |
| `notes` | Anything a careful reader needs, e.g. a contradiction in the source | free text | | mapping rules |

## sample_treatment

One row per compound per sample. E2 plus tamoxifen is two rows.

- Vehicle controls have one row naming the solvent (`ETHANOL` or `DMSO`).
- Untreated controls and parental lines have no rows.

| Column | Meaning | Type / units | Allowed values | Comes from |
|---|---|---|---|---|
| `sample_id` | The sample | GSM ID | | |
| `compound_id` | What was added | text | see `compound` | GEO labels, mapped through `compound_alias` |
| `dose_nM` | Concentration | **nanomolar (nM)**. 1 µM = 1,000 nM. | 0–100,000. NULL when unknown or not applicable. | GEO or paper, converted to nM |
| `dose_status` | How trustworthy the dose is | text | `stated` (in GEO), `from_paper` (GEO silent, paper gives it), `conflicting` (sources disagree; the stored value is the most credible), `unknown` (no source gives it), `not_applicable` (solvent) | |
| `exposure_type` | Short treatment or months-long selection | text | `acute`, `chronic` (used to make resistant lines; `exposure_h` is NULL) | |
| `exposure_h` | How long the compound was on the cells before RNA was taken | **hours** | 0 or more. NULL when unknown or chronic. | GEO time fields |
| `exposure_note` | Extra timing detail | free text, e.g. `8-12 months` | | |
| `label_as_given` | The lab's own words for the treatment | text, e.g. `TOT`, `1uM G-03046274` | | GEO |
| `dose_as_given` | The lab's own words for the dose | text, e.g. `1uM (treatment protocol says '1 M')` | | GEO or paper |
| `source_file_id` | The raw file | text | | |

## compound

One row per real compound. Every spelling the labs used is in `compound_alias`.

| Column | Meaning | Allowed values |
|---|---|---|
| `compound_id` | Key | `OHT` (4-hydroxytamoxifen), `TAM_UNSPECIFIED` (tamoxifen, form not stated), `E2` (estradiol), `HRG_B1` (heregulin-β1), `FULVESTRANT`, `GDC0810`, `GDC0927`, `GNE274`, `ETHANOL`, `DMSO` |
| `preferred_name` | Readable name | |
| `compound_class` | Kind of compound | `SERM` (tamoxifen family: blocks the estrogen receptor), `SERD` (degrades the receptor), `hormone`, `growth_factor`, `vehicle` (solvent only), `other` |
| `is_vehicle` | Solvent used for controls | boolean |
| `pubchem_cid` | PubChem compound ID | NULL when there's no single structure |
| `notes` | Caveats | |

## compound_alias

Every way a compound was written, and where. Use it to see how a lab labeled a compound. Join on
`compound_id`, not on the alias.

| Column | Meaning |
|---|---|
| `alias` | Exact spelling, e.g. `TOT`, `OH-Tam`, `4OH-Tamoxifen`, `G-03046274` |
| `compound_id` | The compound it means |
| `seen_in_study` | Study where it appears |
| `seen_in_field` | Where in that study's files, or "paper only" |

## cell_line

Base lines, with their Cellosaurus IDs, plus the resistant lines the labs made from them.

| Column | Meaning | Allowed values / examples |
|---|---|---|
| `cell_line_id` | Key | base lines `CVCL_0031` MCF-7, `CVCL_0553` T-47D, `CVCL_0179` BT-474, `CVCL_0588` ZR-75-1, `CVCL_1254` HCC1500, `CVCL_0619` MDA-MB-330, `CVCL_1115` CAMA-1, `CVCL_0253` EFM-19; derived lines `S2_TAMR_SELECTED`, `S2_TAMR_C1`–`C6`, `S4_<LINE>_TAM1/2` |
| `preferred_name` | Readable name | |
| `is_derived` | Lab-made from another line | boolean |
| `parent_cell_line_id` | The line it was made from | NULL for base lines |
| `derived_by_study` | Which study's lab made it | |
| `derivation` | How it was made | e.g. `1 uM 4-OHT for 8-12 months` |
| `disease` | Tumor type of the base line (Cellosaurus) | |
| `source_bank` | Where the base line came from | `ATCC`, `DSMZ`, `lab-made` |
| `notes` | Caveats | e.g. S2's time-course clone is one of `S2_TAMR_C1`–`C6`, but the paper never says which |

Important: resistant lines from different labs are different cell lines. S2's `TamR` and S4's `Tam1`
are unrelated.

## cell_line_alias

| Column | Meaning |
|---|---|
| `alias` | Exact spelling, e.g. `MCF7`, `MCF-7cells`, `TamR#3`, `WT` |
| `cell_line_id` | The line it means |
| `seen_in_study` | Study where it appears |

## gene

Every approved human gene from HGNC (downloaded 26 Sep 2026), keyed on the permanent HGNC ID.

| Column | Meaning | Type / examples |
|---|---|---|
| `hgnc_id` | Key | `HGNC:3467` |
| `symbol` | Today's official symbol | `ESR1` |
| `name` | Full name | `estrogen receptor 1` |
| `locus_type` | Kind of gene | e.g. `gene with protein product`, `RNA, long non-coding` |
| `entrez_id`, `ensembl_gene_id` | Other databases' IDs | text |
| `previous_symbols` | Names it used to have | list, e.g. `['WHSC1']` for NSD2 |
| `alias_symbols` | Other common names | list |
| `source_file_id` | `HGNC_2026-09-26` | |

## feature_map

How each study's own row IDs map to genes. Arrays measure probe sets, and RNA-seq measures Ensembl or
Entrez genes.

| Column | Meaning | Allowed values |
|---|---|---|
| `study_id`, `feature_id` | The study and its row ID, e.g. `205225_at`, `NM_000125_at`, `ENSG00000091831`, `2099` | |
| `feature_type` | Kind of row ID | `affy_probe_set` (S1, S3), `refseq_probe_set` (S2), `ensembl_gene` (S4), `entrez_gene` (S5) |
| `symbol_as_given` | The study's own gene name, kept for audit only. **Don't join on it.** | |
| `hgnc_id` | The gene, when matched | |
| `map_status` | Match result | `mapped`, `unmapped`, `multiple_genes` (probe matches several genes; not used), `control_probe`, `corrupt_annotation` |
| `map_method` | Evidence used | `entrez_id`, `refseq_accession`, `ensembl_gene_id`, `current_symbol`, `previous_symbol` |
| `chosen_for_gene` | TRUE for the one row per gene whose values are in `measurement` | boolean |
| `choice_rule` | How that row was chosen | highest mean log2 expression among rows for the gene |

## measurement

Expression values, one per sample per gene.

| Column | Meaning | Type / units | Allowed values |
|---|---|---|---|
| `sample_id` | The sample | GSM ID | |
| `hgnc_id` | The gene | HGNC ID | |
| `value` | Expression, on the study's clean scale. Higher means more. A difference of 1 is a 2-fold change. | **log2** units | |
| `value_type` | Which scale | text | `log2_array_signal` (S1–S3: log2 of array signal; S1 and S3 were converted with log2(x+1)), `log2_cpm` (S4, S5: log2(counts per million + 1), recomputed from raw counts) |
| `raw_value` | The number as deposited | raw units | |
| `raw_value_type` | What the raw number is | text | `linear_array_signal`, `log2_array_signal`, `raw_count` |
| `feature_id` | The study row the value came from | text | see `feature_map` |
| `source_file_id` | The raw file | text | |

`value` is comparable **within a study only**. S4 has one sample per line, so there are no replicates
and no within-line statistics.

## study

| Column | Meaning |
|---|---|
| `study_id` | `S1`–`S5` |
| `geo_accession`, `geo_title` | GEO series ID and title as deposited |
| `paper_title`, `first_author`, `year`, `journal`, `pmid` | The publication. Titles can differ from GEO's. |
| `lab_institution` | Where the work was done |
| `technology` | `microarray` or `rna_seq` |
| `platform`, `platform_description` | GEO platform, e.g. `GPL96`, and a readable description |
| `processing_as_stated` | What GEO says was done to the values. Not always true; see `issue`. |
| `raw_value_scale` | What the deposited values actually are: `log2`, `linear` or `raw_count` |
| `samples_in_geo` | Number of samples |
| `parent_series` | Another GEO series that also contains these samples (S5 is inside GSE117943) |
| `notes` | Key caveats |

## source_file

| Column | Meaning |
|---|---|
| `source_file_id` | e.g. `S3_matrix`, `S4_raw_tar`, `S5_counts`, `HGNC_2026-09-26`, `S2_paper` |
| `path` | Location under the project folder |
| `study_id` | NULL for the shared gene reference |
| `file_kind` | `series_matrix`, `platform_annotation`, `per_sample_bundle`, `counts_table`, `gene_reference`, `paper` |
| `sha256`, `size_bytes` | Checksum and size, proving which exact file was used |
| `origin_url`, `downloaded_on` | Where it came from, and when |

## provenance

Why the database says what it says. There is one row per cleaned field per record.

| Column | Meaning |
|---|---|
| `provenance_id` | Key, numbered in build order |
| `entity_type` | `sample` or `sample_treatment` |
| `entity_id` | e.g. `GSM539725`, or `GSM539727|E2` for a treatment |
| `field` | The column that was filled, e.g. `role`, `dose_nM` |
| `clean_value` | What the database stores |
| `source_file_id`, `source_location` | The file, and where in it, e.g. `!Sample_characteristics_ch1 (treatment)` |
| `raw_text` | Exactly what the source said |
| `rule` | The mapping rule applied, with its mess-inventory number, e.g. `[M18]` |

## issue

Every problem the checks found in the source data.

| Column | Meaning | Allowed values |
|---|---|---|
| `issue_id` | Key | |
| `check_name` | Which check found it | e.g. `fields_contradict`, `no_control_in_batch`, `dose_conflicting`, `dose_unknown`, `compound_form_unstated`, `duplicate_titles`, `value_scale_mismatch`, `date_like_gene_name`, `unknown_value`, `batch_inferred` |
| `inventory_id` | Row in the project's mess inventory, `M01`–`M75` | may be NULL |
| `impact` | How much it matters | `breaks` (wrong answers if ignored), `misleads` (plausible answers that mean something else), `minor` |
| `study_id` | Study | |
| `entity_type`, `entity_id` | What it touches | `study`, `sample`, `sample_treatment`, `feature`, `gene`, `source_file`, plus its ID |
| `field`, `observed`, `expected` | What the source says, and what it should say or conflicts with | |
| `resolution` | What was done | `fixed_by_rule` (clean tables already correct it), `flagged` (kept as is and marked; take care), `excluded` (left out) |
| `message` | One plain-English sentence | |
| `created_at` | When the check ran | |

## Views

- **`v_sample`:** one row per sample with readable `cell_line`, `derived_line` and `genetic_change_target` names, all its treatments spelled out in one `treatments` string (e.g. `4-hydroxytamoxifen 10.0 nM 24.0 h`), and `open_issues`, the number of issues touching it.
- **`v_expression`:** `measurement` joined to gene `symbol` and the sample's `study_id`, `role`, `batch`, `cell_line_id` and `derived_line_id`. Samples not included in analysis are left out.

---

## Five example queries

**1. Which samples got 4-hydroxytamoxifen, at what dose, and for how long?**

```sql
SELECT s.study_id, s.sample_id, cl.preferred_name AS cell_line, d.preferred_name AS derived_line,
       t.dose_nM, t.exposure_type, t.exposure_h, s.batch
FROM sample_treatment t
JOIN sample s USING (sample_id)
JOIN cell_line cl ON cl.cell_line_id = s.cell_line_id
LEFT JOIN cell_line d ON d.cell_line_id = s.derived_line_id
WHERE t.compound_id = 'OHT'
ORDER BY s.study_id, s.sample_id;
```

**2. Find a gene by any name it has ever had.**

```sql
SELECT hgnc_id, symbol, name, previous_symbols
FROM gene
WHERE symbol = 'WHSC1'
   OR list_contains(previous_symbols, 'WHSC1')
   OR list_contains(alias_symbols, 'WHSC1');
```

**3. How much does 24 h of estradiol change GREB1 in ordinary MCF-7, study by study, comparing only within each batch?**

```sql
WITH x AS (
  SELECT s.sample_id, s.study_id, s.batch, s.role, m.value,
         EXISTS (SELECT 1 FROM sample_treatment t WHERE t.sample_id = s.sample_id
                 AND t.compound_id = 'E2' AND t.exposure_h = 24) AS e2_24h,
         (SELECT count(*) FROM sample_treatment t WHERE t.sample_id = s.sample_id
                 AND t.compound_id NOT IN ('E2', 'ETHANOL', 'DMSO')) AS other_compounds
  FROM sample s
  JOIN measurement m USING (sample_id)
  JOIN gene g USING (hgnc_id)
  WHERE g.symbol = 'GREB1'
    AND s.cell_line_id = 'CVCL_0031' AND s.derived_line_id IS NULL
    AND s.genetic_change IN ('none', 'empty_vector', 'non_targeting_control')
    AND s.include_in_analysis
)
SELECT study_id, coalesce(batch, '(one run)') AS batch,
       avg(value) FILTER (WHERE e2_24h AND other_compounds = 0)                AS e2_mean,
       avg(value) FILTER (WHERE role IN ('vehicle_control', 'untreated_control')) AS control_mean,
       round(e2_mean - control_mean, 2)                                         AS log2_change
FROM x
GROUP BY study_id, batch
HAVING e2_mean IS NOT NULL AND control_mean IS NOT NULL
ORDER BY study_id, batch;
```

A `log2_change` of 1 means GREB1 doubled. Change values like this can be compared across studies;
raw `value`s can't.

**4. What's wrong with a sample, and was it handled?**

```sql
SELECT check_name, impact, resolution, field, observed, message
FROM issue
WHERE entity_type = 'sample' AND entity_id = 'GSM539725'
UNION ALL
SELECT check_name, impact, resolution, field, observed, message
FROM issue
WHERE entity_type = 'study' AND study_id = (SELECT study_id FROM sample WHERE sample_id = 'GSM539725')
ORDER BY impact, check_name;
```

**5. Why does the database say what it says about a sample?**

```sql
SELECT field, clean_value, source_file_id, source_location, raw_text, rule
FROM provenance
WHERE entity_id = 'GSM645721' OR entity_id LIKE 'GSM645721|%'
ORDER BY field;
```

---

## Short glossary

- **Sample:** one dish of cells, measured once.
- **Vehicle control:** cells given only the solvent the drug is dissolved in.
- **Untreated control:** nothing added.
- **Parental:** the original line a resistant line was made from.
- **Resistant:** cells that keep growing in tamoxifen. The labs made these by growing cells in the drug for months.
- **Batch:** one experimental run. Samples from different runs differ for technical reasons.
- **log2:** the scale for expression. +1 is double, -1 is half.
- **CPM:** RNA-seq reads for a gene per million reads in the sample.
- **HGNC ID:** the permanent ID of a human gene, which survives renaming.
- **Knockdown:** a gene silenced with siRNA.
- **Overexpression:** extra copies of a gene added, here by a virus.
