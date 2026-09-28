-- =====================================================================================
-- Messy Lab · clean schema (DuckDB)
--
-- One set of tables for five studies of tamoxifen response in MCF-7 and related cells.
-- Rules:
--   * Nothing is guessed. A value the sources don't give is NULL plus a status of 'unknown',
--     and the reason goes in the notes or the issue table.
--   * Every row points back to the raw file it came from (source_file_id), and every
--     standardized field can be traced to the exact field and raw text in `provenance`.
--   * Controlled vocabularies are enforced with CHECK constraints. The allowed values are
--     listed next to each column and in db/DATA_DICTIONARY.md.
--
-- Build order: reference tables first, then study, sample, treatments, measurements, flags.
-- Run:  duckdb db/messy_lab.duckdb < pipeline/schema.sql   (pipeline/build.py does this)
-- =====================================================================================


-- -------------------------------------------------------------------------------------
-- SOURCE FILES: every file the database was built from, untouched, with a checksum.
-- Includes the GEO downloads, the gene reference, and the papers (for facts GEO lacks).
-- -------------------------------------------------------------------------------------
CREATE TABLE source_file (
    source_file_id   VARCHAR PRIMARY KEY,          -- e.g. 'S3_matrix', 'S4_raw_tar', 'HGNC_2026-09-26', 'S3_paper'
    path             VARCHAR NOT NULL UNIQUE,      -- relative to the Messy_Lab folder
    study_id         VARCHAR,                      -- NULL for shared reference files
    file_kind        VARCHAR NOT NULL CHECK (file_kind IN (
                         'series_matrix', 'platform_annotation', 'per_sample_bundle',
                         'counts_table', 'gene_reference', 'paper')),
    sha256           VARCHAR NOT NULL,
    size_bytes       BIGINT  NOT NULL,
    origin_url       VARCHAR,                      -- where it was downloaded from
    downloaded_on    DATE,
    notes            VARCHAR
);


-- -------------------------------------------------------------------------------------
-- STUDY: one row per GEO series.
-- -------------------------------------------------------------------------------------
CREATE TABLE study (
    study_id            VARCHAR PRIMARY KEY,       -- 'S1' .. 'S5'
    geo_accession       VARCHAR NOT NULL UNIQUE,   -- 'GSE4025'
    geo_title           VARCHAR NOT NULL,          -- title as written in GEO
    paper_title         VARCHAR,                   -- title as published (can differ from GEO)
    first_author        VARCHAR,
    year                INTEGER,
    journal             VARCHAR,
    pmid                VARCHAR,
    lab_institution     VARCHAR,                   -- where the work was done (from the paper)
    technology          VARCHAR NOT NULL CHECK (technology IN ('microarray', 'rna_seq')),
    platform            VARCHAR NOT NULL,          -- GEO platform, e.g. 'GPL96'
    platform_description VARCHAR,                  -- e.g. 'Affymetrix U133A', 'HiSeq 2000'
    processing_as_stated VARCHAR,                  -- what GEO says was done to the values
    raw_value_scale     VARCHAR NOT NULL CHECK (raw_value_scale IN (
                            'log2', 'linear', 'raw_count')),  -- what the values actually are
    samples_in_geo      INTEGER NOT NULL,          -- samples GEO lists for this series
    parent_series       VARCHAR,                   -- e.g. 'GSE117943' when this is a SubSeries
    notes               VARCHAR
);


-- -------------------------------------------------------------------------------------
-- REFERENCE: compounds. One row per real compound; every spelling is in compound_alias.
-- Vehicles (ethanol, DMSO) are compounds too, so a vehicle control can name its solvent.
-- -------------------------------------------------------------------------------------
CREATE TABLE compound (
    compound_id      VARCHAR PRIMARY KEY,          -- e.g. 'OHT', 'TAM_UNSPECIFIED', 'E2', 'HRG_B1', 'DMSO'
    preferred_name   VARCHAR NOT NULL,             -- e.g. '4-hydroxytamoxifen'
    compound_class   VARCHAR NOT NULL CHECK (compound_class IN (
                         'SERM', 'SERD', 'hormone', 'growth_factor', 'vehicle', 'other')),
    is_vehicle       BOOLEAN NOT NULL DEFAULT FALSE,
    pubchem_cid      VARCHAR,                      -- NULL when the form is not known
    notes            VARCHAR                       -- e.g. 'tamoxifen, form not stated; never pooled with OHT'
);

CREATE TABLE compound_alias (
    alias            VARCHAR NOT NULL,             -- exactly as written, e.g. 'TOT', 'G-03046274', '4OH-Tamoxifen'
    compound_id      VARCHAR NOT NULL REFERENCES compound (compound_id),
    seen_in_study    VARCHAR NOT NULL,             -- 'S1' .. 'S5', or 'paper' when only the paper uses it
    seen_in_field    VARCHAR,                      -- e.g. 'characteristics_ch1', 'counts header row 2'
    PRIMARY KEY (alias, seen_in_study)
);


-- -------------------------------------------------------------------------------------
-- REFERENCE: cell lines. Base lines carry their Cellosaurus ID. Lines a lab derived
-- (tamoxifen-resistant clones) get their own row that points to the parent.
-- -------------------------------------------------------------------------------------
CREATE TABLE cell_line (
    cell_line_id        VARCHAR PRIMARY KEY,       -- base: Cellosaurus ID 'CVCL_0031'; derived: 'S2_TAMR', 'S4_MCF7_TAM1'
    preferred_name      VARCHAR NOT NULL,          -- 'MCF-7', 'MCF-7 TamR (RIKEN, selected clone)'
    is_derived          BOOLEAN NOT NULL DEFAULT FALSE,
    parent_cell_line_id VARCHAR REFERENCES cell_line (cell_line_id),
    derived_by_study    VARCHAR,                   -- the study whose lab made it
    derivation          VARCHAR,                   -- e.g. '1 uM tamoxifen for 3 months, one of six clones'
    disease             VARCHAR,                   -- e.g. 'invasive ductal carcinoma'
    source_bank         VARCHAR,                   -- 'ATCC', 'DSMZ', 'unknown'
    notes               VARCHAR
);

CREATE TABLE cell_line_alias (
    alias            VARCHAR NOT NULL,             -- 'MCF7', 'MCF-7cells', 'T47 D', 'Cama-1', 'TamR#3'
    cell_line_id     VARCHAR NOT NULL REFERENCES cell_line (cell_line_id),
    seen_in_study    VARCHAR NOT NULL,
    PRIMARY KEY (alias, seen_in_study)
);


-- -------------------------------------------------------------------------------------
-- REFERENCE: genes, keyed on the HGNC ID so renamed genes and Excel-damaged names
-- (e.g. '1-Mar', '7-Sep') can't split or merge genes.
-- -------------------------------------------------------------------------------------
CREATE TABLE gene (
    hgnc_id          VARCHAR PRIMARY KEY,          -- 'HGNC:3467'
    symbol           VARCHAR NOT NULL,             -- current official symbol, e.g. 'ESR1'
    name             VARCHAR,
    locus_type       VARCHAR,                      -- e.g. 'gene with protein product'
    entrez_id        VARCHAR,
    ensembl_gene_id  VARCHAR,
    previous_symbols VARCHAR[],                    -- e.g. ['WHSC1'] for NSD2
    alias_symbols    VARCHAR[],
    source_file_id   VARCHAR NOT NULL REFERENCES source_file (source_file_id)
);

-- How each study's own row IDs (probe sets, RefSeq, Ensembl, Entrez) map to genes.
-- Unmapped and ambiguous rows are kept here with a status, never silently dropped.
CREATE TABLE feature_map (
    study_id         VARCHAR NOT NULL REFERENCES study (study_id),
    feature_id       VARCHAR NOT NULL,             -- '201291_s_at', 'NM_000014_at', 'ENSG00000091831', '2099'
    feature_type     VARCHAR NOT NULL CHECK (feature_type IN (
                         'affy_probe_set', 'refseq_probe_set', 'ensembl_gene', 'entrez_gene')),
    symbol_as_given  VARCHAR,                      -- the study's own name for it, kept for audit
    hgnc_id          VARCHAR REFERENCES gene (hgnc_id),
    map_status       VARCHAR NOT NULL CHECK (map_status IN (
                         'mapped', 'unmapped', 'multiple_genes', 'control_probe', 'corrupt_annotation')),
    map_method       VARCHAR,                      -- e.g. 'entrez_id', 'refseq_accession', 'previous_symbol'
    chosen_for_gene  BOOLEAN NOT NULL DEFAULT FALSE,  -- the one row per gene that feeds `measurement`
    choice_rule      VARCHAR,                      -- e.g. 'highest mean signal among probes for this gene'
    PRIMARY KEY (study_id, feature_id)
);


-- -------------------------------------------------------------------------------------
-- SAMPLE: one row per GEO sample (GSM ID). Everything about what the sample IS.
-- What was added to it lives in sample_treatment.
-- -------------------------------------------------------------------------------------
CREATE TABLE sample (
    sample_id            VARCHAR PRIMARY KEY,      -- GSM ID, e.g. 'GSM645710'
    study_id             VARCHAR NOT NULL REFERENCES study (study_id),
    geo_title            VARCHAR NOT NULL,         -- as written, e.g. 'MCF-7 +scrambled siRNA + vehicle rep1 (MD31)'
    lab_sample_code      VARCHAR,                  -- the lab's own ID: 'MD31', 'TamR.Ctrl.00h', 'SAM24314537'
    cell_line_id         VARCHAR NOT NULL REFERENCES cell_line (cell_line_id),  -- always the base line
    derived_line_id      VARCHAR REFERENCES cell_line (cell_line_id),           -- resistant derivative, if any
    genetic_change       VARCHAR NOT NULL CHECK (genetic_change IN (
                             'none', 'overexpression', 'empty_vector', 'knockdown', 'non_targeting_control')),
    genetic_change_target VARCHAR REFERENCES gene (hgnc_id),  -- e.g. ESR2 for AdERb, ESR1 or RARA for siRNA
    genetic_change_method VARCHAR,                 -- e.g. 'adenovirus, moi 10', 'siRNA, Dharmafect 1'
    role                 VARCHAR NOT NULL CHECK (role IN (
                             'treated', 'vehicle_control', 'untreated_control', 'parental')),
    resistance_phenotype VARCHAR NOT NULL CHECK (resistance_phenotype IN (
                             'sensitive', 'resistant', 'unknown')),
    batch                VARCHAR,                  -- e.g. 'S3_A', 'S2_wt3', 'S5_run2_inferred'
    batch_basis          VARCHAR CHECK (batch_basis IS NULL OR batch_basis IN ('stated', 'inferred')),
    replicate            INTEGER,                  -- 1, 2, 3 within its group, assigned when the source is silent
    replicate_basis      VARCHAR CHECK (replicate_basis IS NULL OR replicate_basis IN ('stated', 'assigned')),
    hormone_deprived     VARCHAR NOT NULL CHECK (hormone_deprived IN ('yes', 'no', 'unknown')),
    in_drug_at_harvest   VARCHAR NOT NULL CHECK (in_drug_at_harvest IN ('yes', 'no', 'unknown')),
    culture_medium       VARCHAR,                  -- e.g. 'DMEM, phenol red free, 5% charcoal-stripped FBS'
    harvest_note         VARCHAR,                  -- e.g. 'RNA 72 h after siRNA transfection (paper)'
    include_in_analysis  BOOLEAN NOT NULL DEFAULT TRUE,
    exclusion_reason     VARCHAR,                  -- required when include_in_analysis is FALSE
    source_file_id       VARCHAR NOT NULL REFERENCES source_file (source_file_id),
    notes                VARCHAR,
    CHECK (include_in_analysis OR exclusion_reason IS NOT NULL)
);


-- -------------------------------------------------------------------------------------
-- SAMPLE_TREATMENT: one row per compound per sample. 'E2 + tamoxifen' is two rows.
-- Vehicle controls get one row naming the solvent. Untreated controls and parental
-- lines with nothing added get no rows.
-- -------------------------------------------------------------------------------------
CREATE TABLE sample_treatment (
    sample_id        VARCHAR NOT NULL REFERENCES sample (sample_id),
    compound_id      VARCHAR NOT NULL REFERENCES compound (compound_id),
    dose_nM          DOUBLE,                       -- NULL when unknown or not applicable
    dose_status      VARCHAR NOT NULL CHECK (dose_status IN (
                         'stated', 'from_paper', 'conflicting', 'unknown', 'not_applicable')),
    exposure_type    VARCHAR NOT NULL CHECK (exposure_type IN ('acute', 'chronic')),
    exposure_h       DOUBLE,                       -- hours for acute; NULL for chronic or unknown
    exposure_note    VARCHAR,                      -- e.g. '8-12 months', '0 h: baseline, drug never added'
    label_as_given   VARCHAR NOT NULL,             -- raw text, e.g. 'TOT', '1uM G-03046274', 'E2, tamoxifen'
    dose_as_given    VARCHAR,                      -- raw text, e.g. '1 M', '10nM', 'n/a'
    source_file_id   VARCHAR NOT NULL REFERENCES source_file (source_file_id),
    PRIMARY KEY (sample_id, compound_id),
    CHECK (dose_nM IS NULL OR (dose_nM >= 0 AND dose_nM <= 100000)),   -- above 100 uM is a data error
    CHECK (exposure_h IS NULL OR exposure_h >= 0)
);


-- -------------------------------------------------------------------------------------
-- MEASUREMENT: expression, long format. One value per sample per gene, on one scale.
--   arrays:  log2 signal (linear studies converted with log2(x + 1))
--   RNA-seq: log2 counts per million, recomputed from raw counts within each study
-- The raw number is kept beside it. Values are comparable WITHIN a study, not across.
-- -------------------------------------------------------------------------------------
CREATE TABLE measurement (
    sample_id        VARCHAR NOT NULL REFERENCES sample (sample_id),
    hgnc_id          VARCHAR NOT NULL REFERENCES gene (hgnc_id),
    value            DOUBLE NOT NULL,
    value_type       VARCHAR NOT NULL CHECK (value_type IN ('log2_array_signal', 'log2_cpm')),
    raw_value        DOUBLE NOT NULL,
    raw_value_type   VARCHAR NOT NULL CHECK (raw_value_type IN (
                         'linear_array_signal', 'log2_array_signal', 'raw_count')),
    feature_id       VARCHAR NOT NULL,             -- the probe set or gene row the value came from
    source_file_id   VARCHAR NOT NULL REFERENCES source_file (source_file_id),
    PRIMARY KEY (sample_id, hgnc_id)
);


-- -------------------------------------------------------------------------------------
-- PROVENANCE: for every standardized field, where it came from and what it said.
-- Answers "why does the database say this?" for any cell.
-- -------------------------------------------------------------------------------------
CREATE SEQUENCE provenance_seq START 1;
CREATE TABLE provenance (
    provenance_id    INTEGER PRIMARY KEY DEFAULT nextval('provenance_seq'),
    entity_type      VARCHAR NOT NULL CHECK (entity_type IN (
                         'study', 'sample', 'sample_treatment', 'cell_line', 'compound')),
    entity_id        VARCHAR NOT NULL,             -- e.g. 'GSM539725' or 'GSM539725|E2'
    field            VARCHAR NOT NULL,             -- column in the clean table, e.g. 'role'
    clean_value      VARCHAR,                      -- what the database stores
    source_file_id   VARCHAR NOT NULL REFERENCES source_file (source_file_id),
    source_location  VARCHAR NOT NULL,             -- e.g. '!Sample_characteristics_ch1', 'Methods, Cell culture'
    raw_text         VARCHAR,                      -- exactly what the source said
    rule             VARCHAR                       -- mapping rule applied, e.g. 'S2: 0h -> untreated baseline'
);


-- -------------------------------------------------------------------------------------
-- ISSUE: every problem the checks find (step 3). Each links to its mess-inventory row.
-- -------------------------------------------------------------------------------------
CREATE SEQUENCE issue_seq START 1;
CREATE TABLE issue (
    issue_id         INTEGER PRIMARY KEY DEFAULT nextval('issue_seq'),
    check_name       VARCHAR NOT NULL,             -- e.g. 'dose_above_100uM', 'fields_contradict'
    inventory_id     VARCHAR,                      -- 'M01' .. 'M75'; NULL for a problem the inventory missed
    impact           VARCHAR NOT NULL CHECK (impact IN ('breaks', 'misleads', 'minor')),
    study_id         VARCHAR REFERENCES study (study_id),
    entity_type      VARCHAR NOT NULL CHECK (entity_type IN (
                         'study', 'sample', 'sample_treatment', 'feature', 'gene', 'source_file')),
    entity_id        VARCHAR NOT NULL,
    field            VARCHAR,
    observed         VARCHAR,                      -- what the source says
    expected         VARCHAR,                      -- what it should say, or what it conflicts with
    resolution       VARCHAR NOT NULL CHECK (resolution IN (
                         'fixed_by_rule', 'flagged', 'excluded')),
    message          VARCHAR NOT NULL,             -- one plain-English sentence
    created_at       TIMESTAMP NOT NULL DEFAULT current_timestamp
);


-- -------------------------------------------------------------------------------------
-- VIEWS: the two shapes most questions need, already joined.
-- -------------------------------------------------------------------------------------

-- One row per sample, with its line and every treatment spelled out.
CREATE VIEW v_sample AS
SELECT
    s.sample_id,
    s.study_id,
    st.geo_accession,
    s.geo_title,
    base.preferred_name                      AS cell_line,
    der.preferred_name                       AS derived_line,
    s.genetic_change,
    g.symbol                                 AS genetic_change_target,
    s.role,
    s.resistance_phenotype,
    s.batch,
    s.replicate,
    s.hormone_deprived,
    s.in_drug_at_harvest,
    s.include_in_analysis,
    string_agg(c.preferred_name
               || coalesce(' ' || CAST(t.dose_nM AS VARCHAR) || ' nM', '')
               || CASE WHEN t.exposure_type = 'chronic' THEN ' (chronic)'
                       ELSE coalesce(' ' || CAST(t.exposure_h AS VARCHAR) || ' h', '') END,
               ' + ' ORDER BY c.preferred_name) AS treatments,
    (SELECT count(*) FROM issue i
      WHERE i.entity_type = 'sample' AND i.entity_id = s.sample_id) AS open_issues
FROM sample s
JOIN study st            ON st.study_id = s.study_id
JOIN cell_line base      ON base.cell_line_id = s.cell_line_id
LEFT JOIN cell_line der  ON der.cell_line_id = s.derived_line_id
LEFT JOIN gene g         ON g.hgnc_id = s.genetic_change_target
LEFT JOIN sample_treatment t ON t.sample_id = s.sample_id
LEFT JOIN compound c     ON c.compound_id = t.compound_id
GROUP BY ALL;

-- Expression with the gene symbol and the sample context needed to compare within a study.
CREATE VIEW v_expression AS
SELECT
    m.sample_id,
    s.study_id,
    g.symbol,
    m.hgnc_id,
    m.value,
    m.value_type,
    s.role,
    s.batch,
    s.cell_line_id,
    s.derived_line_id
FROM measurement m
JOIN sample s ON s.sample_id = m.sample_id
JOIN gene g   ON g.hgnc_id  = m.hgnc_id
WHERE s.include_in_analysis;
