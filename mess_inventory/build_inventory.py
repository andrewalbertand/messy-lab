"""
Messy Lab, step 1: the mess inventory.

One list of every way the five studies record the same things differently, or record them wrongly.
This file is the single source. Running it writes mess_inventory.csv (for the pipeline) and
mess_inventory.md (for reading). Evidence for each row is in evidence_output.txt (from
profile_studies.py) or the named paper.

    python3 mess_inventory/build_inventory.py
"""
import csv
from pathlib import Path

HERE = Path(__file__).parent

STUDIES = [
    ("S1", "GSE4025", "Frasor 2006, Cancer Research", "University of Illinois", "Microarray, Affymetrix U133A", 17),
    ("S2", "GSE21618", "Oyama 2011, J Biol Chem", "RIKEN, Japan", "Microarray, U133 Plus 2 read with a custom RefSeq layout", 143),
    ("S3", "GSE26298", "Salazar 2011, Breast Cancer Research", "University of Toledo", "Microarray, Affymetrix U133 Plus 2", 12),
    ("S4", "GSE111151", "Hultsch 2018, BMC Cancer", "University of Helsinki", "RNA-seq, Illumina HiSeq 2000", 11),
    ("S5", "GSE117942", "Guan 2019, Cell", "Genentech", "RNA-seq, Illumina HiSeq 4000", 82),
]

# impact: Breaks = gives wrong answers if ignored. Misleads = answers look fine but mean something else. Minor = cosmetic or bookkeeping.
F = ["id", "category", "study", "field", "what_they_wrote", "what_it_means", "fix", "impact", "evidence"]
ROWS = [
    # ---------------- DRUG NAMES ----------------
    ("Drug", "S1", "Sample title, characteristics", '"TOT"', "trans-4-hydroxytamoxifen, the active form of tamoxifen. Only the sample protocol spells it out.", "drug = 4-hydroxytamoxifen (4-OHT)", "Breaks", "S1 sample protocol"),
    ("Drug", "S1", "Paper", '"Tam"', "The paper uses \"Tam\" to mean trans-hydroxytamoxifen, not tamoxifen itself.", "drug = 4-OHT", "Misleads", "paper_1.pdf, Methods"),
    ("Drug", "S2", "Treatment field, titles, codes", '"tamoxifen", "Tamoxifen", "Tam" (in codes like E2_Tam)', "The paper says only \"tamoxifen (Sigma-Aldrich)\", never 4-hydroxytamoxifen. It may be the parent drug, which is a different compound from the one the other four studies used.", "drug = tamoxifen, form not stated, flagged. Never pooled with 4-OHT without the flag.", "Misleads", "S2 treatment field; paper_2.pdf, Methods"),
    ("Drug", "S3", "Characteristics, titles, protocol", '"4-hydroxytamoxifen", "OH-Tam"', "4-OHT", "drug = 4-OHT", "Minor", "S3 sample sheet"),
    ("Drug", "S4", "Description, protocol, titles", '"4-OH-tamoxifen", resistant lines named "Tam1", "Tam2"', "4-OHT, given for months to make resistant lines", "drug = 4-OHT, exposure = chronic", "Minor", "S4 sample sheet"),
    ("Drug", "S5", "Series matrix vs counts file", '"4-OH tamoxifen" (matrix), "4OH-Tamoxifen" (counts file)', "Same drug, spelled two ways in two files of one study.", "drug = 4-OHT", "Minor", "S5 counts-file vs matrix label table"),
    ("Drug", "S5", "Series matrix vs counts file", '"GNE-274" (matrix), "G-03046274" (counts file)', "The same compound under two different names. Nothing in the files says so. You only know by matching sample IDs.", "One compound record with both names as aliases", "Breaks", "S5 counts-file vs matrix label table"),
    ("Drug", "S2", "Treatment field", '"HRG"', "Heregulin beta-1, a growth factor, not a drug. 10 nM per the paper.", "ligand = heregulin-beta1 (HRG-b1 176-246)", "Minor", "paper_2.pdf, Methods"),
    ("Drug", "S1, S2, S3, S5", "Treatment fields", '"E2", "17beta-estradiol"', "Estradiol, the natural hormone. Some samples are E2 plus a drug, written \"E2+Tamoxifen\", \"E2, tamoxifen\" or \"E2_Tam\".", "One row per compound per sample, so combinations are two rows, not one string", "Breaks", "S2 sample sheet"),

    # ---------------- DOSE ----------------
    ("Dose", "S1", "Sample fields", "No dose on any sample. \"10nM\" appears only in the free-text protocol.", "10 nM for E2 and for 4-OHT", "dose_nM = 10", "Breaks", "S1 series overall design"),
    ("Dose", "S2", "GEO and paper", "No dose anywhere in GEO. The paper says E2 at 10 nM in one section and 100 nM in another. The tamoxifen dose used in the time course is never stated.", "E2 dose is contradictory. Tamoxifen dose is unknown.", "E2 dose_nM = 10 with a conflict flag. Tamoxifen dose = unknown.", "Breaks", "paper_2.pdf, Methods (Cell culture vs Gene expression analysis)"),
    ("Dose", "S3", "Characteristics", '"100 nM", "500 nM", "n/a"', "The paper mentions only 100 nM. The 500 nM samples exist only in GEO.", "dose_nM as a number. Vehicle = 0. 500 nM flagged as not in the paper.", "Misleads", "S3 sample sheet; paper_3.pdf"),
    ("Dose", "S5", "Treatment protocol", '"1 M 4-OH tamoxifen"', "The micro sign was lost. As written it says one molar, a million times the real dose. The sample fields say 1uM.", "dose_nM = 1000. Reject any dose that parses above 100 uM.", "Breaks", "S5 treatment_protocol"),
    ("Dose", "S5", "Counts file vs matrix", 'Dose inside the label ("1uM GDC-0810") vs its own field ("concentration: 1uM"), and "concentration: 0" for DMSO', "Same information stored two ways", "dose_nM parsed from either, stored once", "Minor", "S5 counts-file vs matrix label table"),
    ("Dose", "S4", "Protocol and description", '"1uM" in the descriptions, "1 µM" in the summary', "1 uM", "dose_nM = 1000", "Minor", "S4 sample sheet"),

    # ---------------- CONTROLS ----------------
    ("Control", "S1", "Characteristics", '"vehicle" and "veh" in the same study', "0.1% ethanol, per the paper", "role = vehicle control, solvent = ethanol", "Minor", "S1 sample sheet; paper_1.pdf"),
    ("Control", "S2", "Title, treatment, source, description, code", '"Control", "untreated", "without stimulus", "unstimulated", "Ctrl"', "Five words for the same thing in one study", "role = untreated control", "Minor", "S2 sample sheet"),
    ("Control", "S2", "Treatment field", '16 of the 17 samples at "0h" carry a drug label, for example "WT E2 0h" = treatment E2', "At 0 hours the drug was never added. These are baselines, not treated samples.", "exposure_h = 0 means untreated, whatever the label says", "Breaks", "S2 time field"),
    ("Control", "S2", "Sample set", "Wild-type cells have 1 untreated sample. Resistant cells have a full untreated time course.", "Wild-type changes over time can only be measured against 0h samples", "Controls matched per group, gaps recorded", "Misleads", "S2 sample sheet"),
    ("Control", "S3", "Characteristics", '"control (vehicle)", "dose: n/a", siRNA "control (scrambled siRNA)"', "Ethanol vehicle plus a non-targeting siRNA. Two separate control dimensions.", "role = vehicle control plus siRNA control", "Minor", "S3 sample sheet"),
    ("Control", "S5", "Matrix vs counts file", '"DMSO", "DMSO control", "concentration: 0"', "DMSO vehicle. Note S1 and S3 used ethanol.", "role = vehicle control, solvent = DMSO", "Minor", "S5 label table"),
    ("Control", "S4", "Description", "Blank for the 4 parental lines", "Untreated parental cells", "role = parental control", "Minor", "S4 sample sheet"),

    # ---------------- TIME ----------------
    ("Time", "S1", "Characteristics", '"24hr"', "24 h. The paper also ran 4-hour arrays that were never deposited.", "exposure_h = 24", "Misleads", "S1 sample sheet; paper_1.pdf, Methods"),
    ("Time", "S2", "Time field, codes", '"0h" to "48h" in fields, "00h" in codes, blank for 7 control samples', "Time course. The 7 blanks are unknown.", "exposure_h as a number, blank = unknown and flagged", "Breaks", "S2 time field"),
    ("Time", "S3", "GEO and paper", "No harvest time anywhere in GEO", "The paper's figure legend says RNA was taken 72 hours after transfection", "harvest_h_after_transfection = 72, source = paper", "Breaks", "paper_3.pdf, Figure legend"),
    ("Time", "S4", "Description", '"eight to twelve months"', "Chronic exposure, no acute time point", "exposure = chronic, 8-12 months", "Misleads", "S4 description"),
    ("Time", "S5", "Characteristics, protocol", '"24h", "24 hours"', "24 h", "exposure_h = 24", "Minor", "S5 sample sheet"),

    # ---------------- CELL LINE AND MODEL ----------------
    ("Cell line", "All", "Various", '"MCF-7 breast cancer cells", "cell line: MCF-7", "MCF-7 breast cancer cell line", "MCF7", "MCF-7cells"', "All the same line", "cell_line = MCF-7, Cellosaurus CVCL_0031", "Breaks", "all sample sheets"),
    ("Cell line", "S4, S5", "Various", '"T-47D", "T47D", "T47 D"; "BT-474", "BT474"; "CAMA-1", "Cama-1"', "Same lines, different spellings", "Cellosaurus IDs for every line", "Breaks", "S4, S5 sample sheets; paper_3.pdf"),
    ("Cell line", "All", "Source name", '"MCF-7 breast cancer cells", "breast epithelial cell, adherent, E2, 24h", "MCF-7cells treated with...", "mammary gland", "cell line"', "Each lab uses this field for something different: the line, the tissue, or the whole condition", "Not used for identity. Kept as raw text only.", "Misleads", "all sample sheets"),
    ("Cell line", "S1", "Title, characteristics", '"Ad", "AdERb"', "Half the samples are MCF-7 with extra estrogen receptor beta added by a virus. \"Ad\" is the empty-virus control. Neither is plain MCF-7.", "genetic_modification = ERb overexpression or empty vector", "Breaks", "S1 sample sheet"),
    ("Cell line", "S3", "Characteristics", '"ER siRNA", "RARalpha siRNA", "scrambled siRNA"', "Knockdowns of the estrogen receptor or of RAR-alpha", "genetic_modification = knockdown target", "Breaks", "S3 sample sheet"),
    ("Cell line", "S2", "Titles, codes, phenotype", '"TamR", "TamR#1" to "TamR#6", "TamR1"; "WT", "wt1" to "wt4"; "resistant", "sensitive"', "The paper made six resistant clones and picked one for the time course, but never says which. TamR#1 to #6 are the six clones, not replicates.", "model = derived resistant line, clone_id kept, parent = MCF-7", "Misleads", "S2 sample sheet; paper_2.pdf, Methods"),
    ("Cell line", "S4", "Titles", '"MCF-7 Tam1" (one MCF-7 resistant line); other lines have Tam1 and Tam2', "Resistant lines from a different lab. S4's \"Tam1\" and S2's \"TamR\" are unrelated cell lines.", "Each derived line gets its own ID with lab and parent", "Breaks", "S4 sample sheet"),
    ("Cell line", "S4", "Paper only", "Resistant lines were kept in 1 uM 4-OHT, parental lines without it", "Resistant vs parental also means drug present vs absent at harvest. GEO does not say this.", "in_drug_at_harvest = yes/no", "Misleads", "paper_4.pdf, Methods"),
    ("Cell line", "S2", "Paper only", "TamR cells routinely kept in 1 uM tamoxifen", "Unknown whether the drug was removed before the time course", "in_drug_at_harvest = unknown, flagged", "Misleads", "paper_2.pdf, Methods"),
    ("Cell line", "S1, S3", "GEO and paper", "GEO gives no source for the MCF-7 cells. S3's paper says ATCC. S1's paper names none.", "Different MCF-7 stocks are known to drift apart", "cell_source recorded where known, else unknown", "Minor", "S1, S3 sample sheets; papers"),

    # ---------------- CULTURE CONDITIONS ----------------
    ("Conditions", "All", "Growth and treatment protocols", "S1 5% charcoal-stripped calf serum. S2 serum-free after charcoal-stripped serum. S3 phenol-red-free, 5% charcoal-stripped serum. S4 full 10% serum plus insulin. S5 hormone deprivation for 3+ days.", "Four studies starved the cells of hormone. S4 did not. That changes the estrogen baseline.", "hormone_deprived = yes/no, medium, serum as fields", "Misleads", "all protocols"),
    ("Conditions", "S2", "GEO vs paper", 'GEO says "serum-free RPMI", starved 16-18 h. Paper says DMEM, starved 16-24 h.', "GEO and the paper disagree on the medium and the timing", "Paper value used, conflict flagged", "Minor", "S2 growth protocol; paper_2.pdf"),
    ("Conditions", "S4", "Growth protocol", '"0,1 % bovine insulin"', "0.1%. European decimal comma.", "Numbers parsed with locale check", "Minor", "S4 growth protocol"),
    ("Conditions", "S1", "Sample protocol vs characteristics", 'Sample protocol: virus at "moi of 5 or 50". Characteristics: "moi10". Series and paper: 10.', "The sample protocol contradicts everything else", "moi = 10, conflict flagged", "Minor", "S1 moi lines"),

    # ---------------- SAMPLES AND REPLICATES ----------------
    ("Samples", "S1", "Titles", '"replicate1 b", "replicate2 b" on 4 of 17 samples', "Meaning of \"b\" is never explained. Possibly re-run arrays.", "Raw title kept, flag = unexplained suffix", "Minor", "S1 sample sheet"),
    ("Samples", "S1", "Sample set", "Empty-virus vehicle has 2 samples. Every other group has 3.", "The paper says three per treatment. One is missing from GEO.", "n per group recorded", "Misleads", "S1 sample sheet; paper_1.pdf"),
    ("Samples", "S2", "Source name vs other fields", 'GSM539725 "WT Control" has source "E2, 24h". GSM539726 "WT E2 24h" has "E2+Tamoxifen". GSM539727 "WT E2+Tamoxifen 24h" has "Tamoxifen". GSM539728 "WT Tamoxifen 0h, 2" has "without stimulus".', "The source field slipped one row for these four samples. Title, treatment field and code agree with each other.", "Majority of fields wins. Conflict flagged on all four.", "Breaks", "S2 disagreement table"),
    ("Samples", "S2", "Titles, codes", '"WT E2 0h rep1", "WT E2 0h", "WT E2 24h, 2"; batches only in codes "wt1" to "wt4"', "Replicates and batches are written three ways, and batch is only in a free-text code", "batch = wt1 to wt4 from codes, replicate index assigned", "Misleads", "S2 sample sheet"),
    ("Samples", "S3", "Titles", 'MD41 and MD42 are both "RARalpha siRNA + vehicle rep1"', "The second should be rep2", "Replicate index from order", "Minor", "S3 identical-titles list"),
    ("Samples", "S3", "Titles, source name", 'MD31, MD32 and MD39, MD40 have identical titles, "scrambled siRNA + vehicle rep1/rep2". The source field lumps all four as one group.', "Two separate experiments. MD39 and MD40 correlate with the RAR-alpha arrays (0.986), MD31 and MD32 with the drug and ER arrays (0.980). The paper calls them \"Control sample\" and \"Control sample 2\".", "batch A = MD31-38, batch B = MD39-42. Each treatment compared only to its own batch control.", "Breaks", "S3 correlation table; paper_3.pdf, Additional files"),
    ("Samples", "S5", "Title", 'Titles are codes, "SAM24314537". No replicate number anywhere.', "Duplicate pairs can only be found by matching cell line and ligand", "Replicate index assigned per cell line and ligand", "Misleads", "S5 sample sheet"),
    ("Samples", "S5", "Sample set", "82 samples. GDC-0927 and GNE-274 were run in 3 of 7 lines. The paper lists 8 lines. MDA-MB-134-VI is not in GEO.", "Unbalanced design, one line missing", "Design recorded per line, gaps listed", "Misleads", "S5 line x ligand table; paper_5.pdf, Methods"),
    ("Samples", "S5", "Read counts", "The 4 MCF-7 GDC-0927 and GNE-274 samples have about 40-43 million reads and IDs from a different range (SAM243220xx). The other 10 MCF-7 samples have 16-25 million (SAM243145xx).", "Very likely a separate sequencing run. Undocumented. Two library kits are named without saying which sample got which.", "batch = inferred run, flagged as inferred", "Misleads", "S5 reads per sample"),
    ("Samples", "S4", "Sample set", "One sample per line, no replicates", "The paper confirms no biological replicates. No statistics within a line.", "n = 1 recorded, blocks any replicate-based test", "Misleads", "S4 sample sheet; paper_4.pdf"),
    ("Samples", "S4", "Series files", "The series download also holds 8 patient tumor files (GSE58708) from another study", "Not cell lines, not this study", "Excluded, recorded as out of scope", "Breaks", "S4 series relations"),
    ("Samples", "S5", "Series relation", '"SubSeries of: GSE117943"', "The same 82 samples also sit under a second accession. Pulling both double counts.", "Deduplicate on GSM ID", "Breaks", "S5 series relations"),
    ("Samples", "All", "IDs", '"GSM92147", "MD31", "TamR.Ctrl.00h", "SAM24314537", SRX and SAMN numbers', "Each study has its own internal ID system on top of GEO's", "GSM ID is the key. Every other ID kept as an alias.", "Minor", "all sample sheets"),

    # ---------------- VALUES ----------------
    ("Values", "S1", "Data values vs processing note", 'Says "GCRMA" (normally log2, about 2 to 15). Values run 1.2 to 19,450.', "The numbers were un-logged. Anyone who assumes log2 is off by a huge margin.", "value_scale = linear, converted to log2", "Breaks", "S1 values summary"),
    ("Values", "S2", "Data values", "RMA, log2, 3.4 to 15.0", "Standard log2 array values", "value_scale = log2", "Minor", "S2 values summary"),
    ("Values", "S3", "Data values vs paper", "GEO: GCOS scaled to housekeeping genes, values 0 to 25,381, 9 zeros. Paper: RMA plus quantile normalization plus a GFP correction.", "The deposited numbers are not the numbers the paper analyzed", "value_scale = linear MAS5-type, log2(x+1), flagged as not the paper's values", "Breaks", "S3 values summary; paper_3.pdf, mRNA profiling"),
    ("Values", "S4", "Per-sample files", "Raw counts, log2 CPM, and batch-adjusted log2 CPM. Normalized together with 8 patient tumors not in this study. GEO says CPM means \"Counts per Megabase\".", "The CPM columns depend on samples we don't have. CPM means counts per million.", "Raw counts only, normalized again inside the clean system", "Breaks", "S4 data_processing"),
    ("Values", "S5", "Counts file vs processing note", "The file is raw counts. The processing note says values are nRPKM.", "The description does not match the file", "Raw counts, normalized again inside the clean system", "Misleads", "S5 data_processing"),
    ("Values", "All", "Platforms", "Three chip designs (22,283, 54,675 and 26,871 probe sets) and two RNA-seq pipelines", "Raw values can't be compared across studies", "Only within-study changes vs a matched control are compared across studies", "Breaks", "values summaries"),

    # ---------------- GENE IDS ----------------
    ("Genes", "S1, S3", "Probe IDs", '"201291_s_at". 1,223 (S1) and 2,214 (S3) probes match several genes. 1,127 and 9,557 match none. ESR1 has 9 probes.', "One gene, many rows, and some rows are ambiguous", "Probes mapped to HGNC IDs with one fixed rule for collapsing", "Breaks", "Gene identifiers section"),
    ("Genes", "S2", "Probe IDs and gene key", '"NM_000014_at". The gene key has no gene-symbol column. The symbol is inside free text: "...(A2M), mRNA".', "Gene names must be pulled out of sentences. 4,963 genes span several rows (ESR1 has 4). 3,433 rows are predicted models (XM_, XR_).", "RefSeq accession mapped to HGNC ID", "Breaks", "Gene identifiers section"),
    ("Genes", "S2", "Gene key", 'Probe XR_015118_at has a description of only "#" characters', "A corrupted entry", "Flagged, left unmapped", "Minor", "Gene identifiers section"),
    ("Genes", "S4", "Per-sample files", 'Ensembl IDs plus names. GEO says "EnsEMBL v80 or v82", the paper says v80. 1,893 gene names repeat.', "Annotation version unclear. Names are not unique.", "Ensembl ID as key, version recorded as v80 (paper)", "Misleads", "S4 per-sample files; paper_4.pdf"),
    ("Genes", "S5", "Counts file", '28 gene names turned into dates by Excel: "1-Mar", "7-Sep", "1-Dec". Two different genes are both "1-Mar", two are both "2-Mar".', "MARCHF1 and MTARC1 became the same \"1-Mar\". SEPTIN7 became \"7-Sep\".", "Entrez ID as key, names taken from HGNC, never from the file", "Breaks", "S5 dates table"),
    ("Genes", "S5", "Counts file", "135 blank gene names, 137 repeated names", "Names are not a usable key", "Entrez ID as key", "Misleads", "S5 counts table"),
    ("Genes", "All", "Gene names", "Renamed genes still using old names: 493 (S1), 2,471 (S2), 1,115 (S3), 2,080 (S4), 1,759 (S5). KMT2D is \"MLL2\" in S2 only. NSD2 is \"WHSC1\" in all five. SEPTIN7 is \"SEPT7\" or \"7-Sep\".", "Joining studies on gene name silently drops or mismatches genes. No study uses today's official names.", "Every gene mapped to an HGNC ID (reference downloaded 26 Sep 2026)", "Breaks", "Gene identifiers section"),
    ("Genes", "All", "Coverage", "Only 10,687 gene names appear in all five studies as written", "Many genes can't be compared across all studies until names are fixed", "Coverage counted after HGNC mapping", "Misleads", "Gene identifiers section"),

    # ---------------- FILE LAYOUT ----------------
    ("Files", "All", "Layout", "S1-S3: one file with the sample sheet on top and values below. S4: labels in one file, values in 11 separate files inside a tar bundle. S5: labels in one file, values in one wide table with 2 label rows above the header and 11 gene-info columns before the samples.", "Three layouts for the same kind of data", "One loader per layout, all writing the same tables", "Breaks", "all raw files"),
    ("Files", "S5", "Counts file", "Sample columns are in a different order from the series matrix", "Matching by position mislabels samples", "Match by sample ID only, never position", "Breaks", "S5 column-order check"),
    ("Files", "S4", "Treatment protocol", '"see above"', "Points to text that isn't there", "Taken from growth protocol and paper", "Minor", "S4 sample sheet"),
    ("Files", "S5", "Paper", '"50 base pair reads", then "reads were trimmed to 75 bp"', "The paper contradicts itself on read length", "Recorded as a conflict, not used", "Minor", "paper_5.pdf, Methods"),
    ("Files", "S1", "Contact institute", '"Kite Pharma"', "The submitter's later employer. The work was done at the University of Illinois.", "Institution taken from the paper", "Minor", "S1 sample sheet; paper_1.pdf"),
    ("Files", "S3", "Contributors", '"Ratnam,,Maya", "Patki,,Mugdha", "d\'Alincourt Salazara"', "First and last names swapped for two authors, one surname misspelled", "Authors taken from the paper", "Minor", "S3 contributors"),
    ("Files", "S2, S3", "Titles", "GEO series titles differ from the paper titles", "Searching by paper title won't find the dataset", "Both titles stored", "Minor", "series titles; papers"),
]

rows = [dict(zip(F, (f"M{i+1:02d}",) + r)) for i, r in enumerate(ROWS)]
with open(HERE / "mess_inventory.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=F); w.writeheader(); w.writerows(rows)

cats = []
for r in rows:
    if r["category"] not in cats: cats.append(r["category"])
count = {k: sum(r["impact"] == k for r in rows) for k in ("Breaks", "Misleads", "Minor")}
esc = lambda s: s.replace("|", "/").replace("\n", " ")

md = []
md.append("# Messy Lab · Mess inventory\n")
md.append("Step 1 of the build plan. Every way the five studies record the same things differently, "
          "or record them wrongly, with the fix the clean system will apply.\n")
md.append("Checked 26 September 2026 against the raw downloads in `studies/data/raw` and the five papers. "
          "Nothing in the raw folder was edited.\n")
md.append("## The five studies\n")
md.append("| # | GEO | Paper | Lab | Technology | Samples |\n|---|---|---|---|---|---|")
for s in STUDIES:
    md.append(f"| {s[0]} | {s[1]} | {s[2]} | {s[3]} | {s[4]} | {s[5]} |")
md.append(f"\nTotal samples: {sum(s[5] for s in STUDIES)}.\n")
md.append("## In short\n")
md.append(f"- **{len(rows)} problems found.** {count['Breaks']} would give wrong answers if ignored, "
          f"{count['Misleads']} would give answers that look right but mean something else, and {count['Minor']} are bookkeeping.")
md.append("- **The biology may agree. The records don't.** The active form of tamoxifen is written 9 different ways, untreated and vehicle controls 10 different ways, "
          "one dose is written as a million times too high, and one study's numbers are un-logged while its label says log.")
md.append("- **Some labels are wrong, not just different.** Four samples in S2 have a mislabeled source field, "
          "and S3 hides two separate experiments under identical control labels.")
md.append("- **Gene names can't be trusted for joining.** Excel turned 28 gene names into dates in S5, "
          "and every study uses gene names that have since been renamed.")
md.append("- **Every problem has a fix.** Each one becomes a rule in the step 2 pipeline.\n")
md.append("## How to read the table\n")
md.append("- **Impact.** *Breaks* means a wrong answer if ignored. *Misleads* means a plausible answer that means something else. *Minor* means bookkeeping.")
md.append("- **Evidence.** Section names refer to `evidence_output.txt`, which `profile_studies.py` regenerates from the raw files. "
          "\"paper_N.pdf\" means the paper in `studies/papers`.\n")
for c in cats:
    md.append(f"## {c}\n")
    md.append("| ID | Study | Field | What they wrote | What it actually means | How the clean system fixes it | Impact | Evidence |")
    md.append("|---|---|---|---|---|---|---|---|")
    for r in rows:
        if r["category"] == c:
            md.append("| " + " | ".join(esc(r[k]) for k in ["id", "study", "field", "what_they_wrote", "what_it_means", "fix", "impact", "evidence"]) + " |")
    md.append("")
md.append("## Files in this folder\n")
md.append("| File | What it is |\n|---|---|")
md.append("| `mess_inventory.md` | This page |")
md.append("| `mess_inventory.csv` | The same rows, for the pipeline to read in step 2 |")
md.append("| `build_inventory.py` | The single source for both. Edit rows here, then run it. |")
md.append("| `profile_studies.py` | Reads the raw files and prints the evidence |")
md.append("| `evidence_output.txt` | Output of `profile_studies.py` |")
md.append("\nReference used for gene names: `studies/data/reference/hgnc_complete_set_2026-09-26.txt` (HGNC, downloaded 26 Sep 2026).\n")
(HERE / "mess_inventory.md").write_text("\n".join(md))
print(f"{len(rows)} rows written | {count}")
