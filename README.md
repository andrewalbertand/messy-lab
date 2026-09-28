# Messy Lab

**Same AI, same questions, cleaner data. Correct answers went from 83% to 99%.**

Messy Lab takes five public gene-expression studies on tamoxifen response in ER-positive breast cancer cells,
finds every way their files disagree, and rebuilds them as one clean, traceable database. Two leading AI models
then answer the same 25 research questions three times each, once with the raw files and once with the clean
database. Every answer is graded automatically against a hand-checked answer key, and every transcript is saved.

| | Raw files | Clean database |
|---|---|---|
| Claude (`claude-opus-5-5`) | 90% | 100% |
| GPT (`gpt-6-astra`) | 75% | 97% |
| Both models | **83%** | **99%** |
| Wrong answers | 24 of 138 | 2 of 138 |
| Average time / cost per answer | 45 s · $0.89 | 15 s · $0.22 |

Scored on 23 of the 25 questions; Q10 and Q23 were dropped after the run as ambiguous (see `runs/summary.md`).
Case study with every transcript: https://claude.ai/artifact/TYBTheiuCcropkEhAHAgqA

## Build the database (one command)

```bash
pip install -r requirements.txt
python3 pipeline/build.py
```

This reads the raw GEO files in `studies/data/raw/` and writes:

- `db/messy_lab.duckdb`, the clean database (265 samples, 5.5 million expression values, 45,111 genes, 399 flagged issues)
- `db/build_fingerprint.txt`, a checksum of every table, so you can confirm your build matches this one
- `db/coverage.md`, how each of the 75 problems in the mess inventory is handled (63 fixed by a rule, 12 flagged, 0 uncovered)

It takes about 1–2 minutes. The raw files are never edited; the build checks their SHA-256 before and after.
Row counts and every table checksum should match `db/build_fingerprint.txt`, except that the `measurement`
checksum can differ in the last floating-point digits across machines and library versions.

Then read `db/DATA_DICTIONARY.md` before querying. Its seven ground rules (compare within a batch, "unknown"
means unknown, 4-OHT is not unstated tamoxifen, genes are keyed on HGNC ID, and so on) are what make the answers right.

## What's here

| Folder | What it does |
|---|---|
| `studies/` | The raw inputs: five GEO series matrices, platform annotations, the per-sample RNA-seq files, and the HGNC gene list. `studies/papers/README.md` links the five papers. |
| `mess_inventory/` | The 75 problems found in the raw files (M01–M75), each with evidence, impact and fix. `mess_inventory.md` is the readable version. |
| `pipeline/` | The build. `schema.sql` defines the tables; `loaders.py` → `standardize.py` (with `mappings/S1–S5.yaml`) → `genes.py` → `values.py` → `checks.py`; `reference/` holds compounds, cell lines and allowed values. `build.py` runs it all. |
| `db/` | Data dictionary, coverage report and build fingerprint. The database itself is built locally. |
| `questions/` | The 25 questions, the saved queries that answer them (`answers.py`), an independent hand check against the raw files and papers (`hand_check.py`), the answer key and grading rules (`make_answer_key.py` → `answer_key.json`), and the grader (`grade.py`, which tests itself when run). |
| `runs/` | The before/after test harness (`harness.py`, models in `models.yaml`), all 300 graded answers with full transcripts (`results/<condition>/<model>/run<k>/Qxx.json`), the re-grader and the summary (`summary.md`). |
| `dashboard/` | A one-page dashboard built from the database: samples, treatments, open issues and a gene lookup. |
| `assistant/` | An assistant that answers questions through the database, cites the samples, files and issues it used, and says when the data can't answer. Includes a 3-question recording. |
| `case_study/` | Source for the case study page. |

## Check the answers yourself

```bash
python3 questions/answers.py        # answers all 25 questions from the database
python3 questions/grade.py          # grader self-test: the key's answers pass, known wrong answers fail
python3 questions/hand_check.py     # independent check from the raw files (needs the papers and pdftotext)
```

## Rerun the AI test

Put `ANTHROPIC_API_KEY` and `OPENAI_API_KEY` in a `.env` file at the top of the repo (it is git-ignored), then:

```bash
python3 runs/harness.py --models claude,gpt --conditions before,after --runs 3 --workers 6
python3 runs/regrade.py && python3 runs/summarize.py
```

Answers already saved are skipped, and failed calls are retried on the next run. A full run of 300 answers
cost about $169 in API fees on 27 September 2026. "Before" gives each model the `studies/` folder as downloaded;
"after" gives it the database and the data dictionary. Instructions, tools and limits are otherwise identical.

## Try the dashboard and the assistant

```bash
python3 dashboard/build_data.py && python3 dashboard/build_page.py    # then open dashboard/index.html
python3 assistant/app.py                                              # then open http://localhost:8765
```

## Data sources

All data is public:

- NCBI GEO series GSE4025, GSE21618, GSE26298, GSE111151 and GSE117942, with their platform files
- The HGNC complete gene set (downloaded 26 September 2026)

The papers are cited, not redistributed. Their PubMed links are in `studies/papers/README.md`, and nothing in
the build needs them.

## About

Built by Andrew Albert, AA Bio Consulting (https://www.aabioconsulting.com). I get life-science companies' data
ready for AI. Code is released under the MIT License (see `LICENSE`).
