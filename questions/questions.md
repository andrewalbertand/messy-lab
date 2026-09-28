# Messy Lab · The 25 questions

These are the questions a biotech team would ask of this data, in their own words. Every model gets
exactly this text in both runs: first with the raw downloads (`studies/`), then with the clean
database (`db/messy_lab.duckdb` + `db/DATA_DICTIONARY.md`).

Each question states the **answer format**, so answers can be graded by code against the answer key
(`questions/answer_key.json`, step 4). Study names are given by GEO accession, the only label a model
sees in both runs:

| GEO ID | Short description |
|---|---|
| GSE4025 | Frasor 2006, arrays |
| GSE21618 | Oyama 2011, RIKEN, arrays |
| GSE26298 | Salazar 2011, Toledo, arrays |
| GSE111151 | Hultsch 2018, Helsinki, RNA-seq |
| GSE117942 | Guan 2019, Genentech, RNA-seq |

"4-OHT" means 4-hydroxytamoxifen. "Ordinary MCF-7" means MCF-7 that is not a tamoxifen-resistant
derivative and has had no gene added or knocked down (empty-vector and scrambled-siRNA controls count
as ordinary).

---

## A. What we have (6)

**Q01.** *(CEO)* How many MCF-7 samples were exposed to 4-OHT, across all five studies? Count
resistant MCF-7 lines grown in 4-OHT. Don't count samples whose tamoxifen form isn't stated.
**Answer format:** one number, plus the count per GEO ID.

**Q02.** *(Head of data)* How many control samples does each study have? A control is a sample given
no drug, hormone or growth factor (only the solvent, or nothing at all), including 0-hour baselines and
parental lines, whatever was done to its genes.
**Answer format:** a count per GEO ID.

**Q03.** *(Head of data)* Which cell lines are in the data, and how many samples does each have? Count
tamoxifen-resistant derivatives under the line they were made from.
**Answer format:** cell line → number of samples.

**Q04.** *(Scientist)* Which compounds were added to cells, and in how many samples each? Count the
solvents too. Keep tamoxifen of an unstated form separate from 4-OHT.
**Answer format:** compound → number of samples.

**Q05.** *(CEO)* How many distinct tamoxifen-resistant cell lines do we have, and which lab made each?
**Answer format:** one number, plus the lines per GEO ID.

**Q06.** *(Head of data)* How many genes can we compare across all five studies?
**Answer format:** one number.

## B. Experimental design (6)

**Q07.** *(Scientist)* Which samples are the right controls for the RAR-alpha knockdown samples in
GSE26298?
**Answer format:** a list of GSM IDs.

**Q08.** *(Scientist)* How many samples in GSE21618 never had anything added to them before their RNA
was taken?
**Answer format:** one number.

**Q09.** *(Scientist)* What concentration of 4-OHT did each study use?
**Answer format:** GEO ID → dose(s) in nM, or "not 4-OHT / not stated".

**Q10.** *(Scientist)* How long were cells exposed to 4-OHT before RNA was collected, in each study that
used 4-OHT?
**Answer format:** GEO ID → hours, "chronic (months)", or "not stated".

**Q11.** *(Head of data)* Which treated samples have no control (solvent-only or untreated) with the same
cells in the same experimental run?
**Answer format:** a list of GSM IDs.

**Q12.** *(Scientist)* Which studies removed estrogen from the growth medium before treating the cells,
and which didn't?
**Answer format:** two lists of GEO IDs, deprived and not deprived.

## C. Quality and conflicts (5)

**Q13.** *(Head of data)* Which samples have labels that contradict each other?
**Answer format:** a list of GSM IDs.

**Q14.** *(Head of data)* Is any dose written at an impossible level anywhere in the data? If so, where,
and what should it be?
**Answer format:** GEO ID, the dose as written, and the corrected dose in nM (or "none").

**Q15.** *(Head of data)* In which studies are the deposited expression values not what GEO's own
processing description says they are?
**Answer format:** a list of GEO IDs.

**Q16.** *(Head of data)* How many gene names have been corrupted into dates by Excel, and in which study?
**Answer format:** one number and a GEO ID.

**Q17.** *(Scientist)* Which samples carry a replicate label that is repeated or wrong?
**Answer format:** a list of GSM IDs.

## D. Biology (5)

**Q18.** *(Scientist)* Is GREB1 induced by 24 hours of 4-OHT in ordinary MCF-7, compared with its own
control, in every study that tested this?
**Answer format:** yes or no, plus the log2 change per GEO ID.

**Q19.** *(Scientist)* In GSE117942, which compound raised GREB1 in every cell line, compared with DMSO
from the same run?
**Answer format:** a list of compounds.

**Q20.** *(Scientist)* Does 24 hours of estradiol (E2) raise GREB1 in ordinary MCF-7 in every study
that tested it, compared within each experimental run?
**Answer format:** yes or no.

**Q21.** *(Scientist)* In GSE26298, how much lower is ESR1 in the ESR1-knockdown samples than in their
own experiment's control?
**Answer format:** a log2 difference.

**Q22.** *(Scientist)* In GSE111151, which tamoxifen-resistant line shows the largest drop in ESR1
compared with its parental line?
**Answer format:** one line name.

## E. Combining studies (3)

**Q23.** *(CEO)* Which studies can't simply be pooled with the others to study 4-OHT response, and why?
**Answer format:** GEO IDs, each with the reason(s).

**Q24.** *(Head of data)* Which studies measured ordinary MCF-7 given 4-OHT for a stated 24 hours with a
matched control in the same run, so that their responses can be combined?
**Answer format:** a list of GEO IDs.

**Q25.** *(CEO)* Combining only the studies that can be combined, what is the average change in GREB1
after 24 hours of 4-OHT in ordinary MCF-7?
**Answer format:** a log2 change, the average of each study's own treated-vs-control change.


**Not scored (dropped 27 Sep, after the full run):** Q10 (answer key too strict: 48 h is a reasonable reading of paper_3 for GSE26298) and Q23 (ambiguous wording; "all five studies" is a defensible answer). Both were still asked; their answers are saved in `runs/results/`.
