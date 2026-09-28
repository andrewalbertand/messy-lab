# Messy Lab — step 5 results

Models (picked 2026-09-27): claude = `claude-opus-5-5`, gpt = `gpt-6-astra`.
Each model answered 25 questions × 3 runs × 2 conditions; 300 answers in all, graded automatically against `questions/answer_key.json`.
**Scored on 23 questions (276 answers).** Q10 dropped: key too strict: the paper supports 48 h as a reasonable reading for GSE26298, where the key said 'not stated'. Q23 dropped: ambiguous: 'can't simply be pooled' has a defensible 'all five studies' reading; missed in all 12 runs, before and after. Both were dropped after seeing the results; their answers stay in `runs/results/`.

## Results

| Model | Before (raw files) | After (clean database) | Change | Avg cost per answer (before → after) | Avg time per answer (before → after) |
|---|---|---|---|---|---|
| claude | 90% (62/69) | 100% (69/69) | +10 pts | $0.73 → $0.19 | 44 s → 17 s |
| gpt | 75% (52/69) | 97% (67/69) | +22 pts | $1.06 → $0.24 | 46 s → 14 s |

**Headline:** on the raw public files the two models were right 83% of the time; on the clean database, 99%. Wrong answers fell from 24 to 2 out of 138, and each answer cost 4.1× less.

Spend: $168.63 on the 300 graded answers, plus $0.85 on attempts cut off by API errors and re-run (kept in `runs/failed_attempts/`).

## By question (correct out of 3 runs; Q10 and Q23 not scored)

| Q | Claude before | GPT before | Claude after | GPT after |
|---|---|---|---|---|
| Q01 | 3 | 3 | 3 | 3 |
| Q02 | 3 | 3 | 3 | 3 |
| Q03 | 3 | 3 | 3 | 3 |
| Q04 | 0 | 0 | 3 | 3 |
| Q05 | 3 | 3 | 3 | 3 |
| Q06 | 2 | 0 | 3 | 3 |
| Q07 | 3 | 3 | 3 | 3 |
| Q08 | 3 | 0 | 3 | 3 |
| Q09 | 3 | 3 | 3 | 3 |
| Q11 | 2 | 0 | 3 | 1 |
| Q12 | 3 | 3 | 3 | 3 |
| Q13 | 2 | 3 | 3 | 3 |
| Q14 | 3 | 3 | 3 | 3 |
| Q15 | 3 | 3 | 3 | 3 |
| Q16 | 3 | 3 | 3 | 3 |
| Q17 | 2 | 0 | 3 | 3 |
| Q18 | 3 | 3 | 3 | 3 |
| Q19 | 3 | 3 | 3 | 3 |
| Q20 | 3 | 3 | 3 | 3 |
| Q21 | 3 | 1 | 3 | 3 |
| Q22 | 3 | 3 | 3 | 3 |
| Q24 | 3 | 3 | 3 | 3 |
| Q25 | 3 | 3 | 3 | 3 |

**Right on the raw files too (6/6 before):** 16 of 23 — Q01, Q02, Q03, Q05, Q07, Q09, Q12, Q14, Q15, Q16, Q18, Q19, Q20, Q22, Q24, Q25.

## Five real wrong answers from the "before" runs

1. **Q04, both models, 6 of 6 runs.** "Which compounds were added, and to how many samples?" Both said E2 80, tamoxifen 59
   and heregulin 46; the answers are 76, 49 and 42. In GSE21618 the 0 h samples carry a drug name ("WT E2 0h"), but at
   0 hours nothing had been added yet. Both models counted them as dosed. Inventory **M18**.
2. **Q08, GPT, 3 of 3 runs: said 7, answer is 31.** "How many GSE21618 samples never had anything added?" Same trap from
   the other side: GPT counted the untreated samples plus only the wild-type 0 h ones, and left out every drug-labelled
   0 h sample. Claude got this right all 3 times. Inventory **M18**.
3. **Q13, Claude run 1: named GSM645720/21, answer is GSM539725–28.** "Which samples have labels that contradict each other?"
   In GSE21618 the source field slipped one row for four samples (GSM539725 is titled "WT Control" but its source says
   "E2, 24h"). Claude compared titles with each other and missed the source column. Inventory **M44**.
4. **Q21, GPT runs 2 and 3: said −0.72, answer is −2.82.** "How much lower is ESR1 after ESR1 knockdown in GSE26298?"
   GPT averaged all nine array probes labelled ESR1. Most barely measure ESR1, so averaging shrank a 7-fold drop to under
   2-fold. The database keeps one probe per gene by a fixed rule. Inventory **M61**.
5. **Q11, GPT 3 of 3 runs and Claude run 2.** "Which treated samples have no control in the same run?" The answer is 8
   GSE117942 samples (the GDC-0927 and GNE-274 arms) sequenced in a separate, undocumented run with no DMSO sample of their
   own. Nothing in GEO says so; it shows only in the sample-ID range and read depth. GPT never found them; instead it listed
   72 GSE21618 samples (reading the 0 h baselines as dosed again) and the 7 GSE111151 resistant lines. Claude found the 8
   every time but in run 2 added the same 79 extras. Inventory **M50** (and **M18**).

## Notes

- **Grader fix, 27 Sep (before summarizing):** the first grading marked correct answers wrong when a model wrote a zero
  entry ("GSE21618: 0") or a name with a code in brackets ("17beta-estradiol (E2)"). `questions/grade.py` now accepts both;
  its self-test gained those real forms plus a new trap (25/25 right, 15/15 traps caught). `runs/regrade.py` re-graded all
  300 saved answers with no new API calls; 21 grades changed, all wrong to right, in both conditions.
- Every answer's full transcript (each code call and its output) and score is in `runs/results/<condition>/<model>/run<k>/Qxx.json`.
- Gemini was dropped after the pilot (250 requests/day cap); its 2-question pilot files remain under `results/*/gemini/`.

