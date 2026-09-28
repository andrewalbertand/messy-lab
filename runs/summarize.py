"""Messy Lab step 5: build runs/summary.md from the saved, graded answers in runs/results/."""
import glob, json, statistics as st, yaml
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODELS = ["claude", "gpt"]
CONDS = ["before", "after"]
spec = yaml.safe_load(open(ROOT / "runs/models.yaml"))
Qtext = {q["id"]: q["question"] for q in json.load(open(ROOT / "questions/answer_key.json"))["questions"]}

# Dropped from scoring on 27 Sep, after the full run, because the question or key was ambiguous (answers still saved).
EXCLUDED = {
    "Q10": "key too strict: the paper supports 48 h as a reasonable reading for GSE26298, where the key said 'not stated'",
    "Q23": "ambiguous: 'can't simply be pooled' has a defensible 'all five studies' reading; missed in all 12 runs, before and after",
}
recs = defaultdict(list)
for m in MODELS:
    for c in CONDS:
        for f in sorted(glob.glob(str(ROOT / f"runs/results/{c}/{m}/run*/Q??.json"))):
            r = json.load(open(f))
            if r["question"] not in EXCLUDED:
                recs[c, m].append(r)
failed_cost = sum(json.load(open(f)).get("cost_usd", 0) for f in glob.glob(str(ROOT / "runs/failed_attempts/**/*.json"), recursive=True))

WRONG = """## Five real wrong answers from the "before" runs

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
"""


def pct(rs): return 100 * sum(r["correct"] for r in rs) / len(rs)

L = ["# Messy Lab — step 5 results", ""]
mids = {k: v for k, v in (spec.get("models") or spec).items() if k in MODELS} if isinstance(spec, dict) else {}
L += [f"Models (picked {spec.get('picked_on', '27 Sep 2026')}): " + ", ".join(f"{k} = `{v.get('model', '?')}`" for k, v in mids.items()) + ".",
      "Each model answered 25 questions × 3 runs × 2 conditions; 300 answers in all, graded automatically against `questions/answer_key.json`.",
      f"**Scored on {25 - len(EXCLUDED)} questions (276 answers).** " + " ".join(f"{q} dropped: {why}." for q, why in EXCLUDED.items())
      + " Both were dropped after seeing the results; their answers stay in `runs/results/`.", ""]

L += ["## Results", "", "| Model | Before (raw files) | After (clean database) | Change | Avg cost per answer (before → after) | Avg time per answer (before → after) |",
      "|---|---|---|---|---|---|"]
for m in MODELS:
    b, a = recs["before", m], recs["after", m]
    L.append(f"| {m} | {pct(b):.0f}% ({sum(r['correct'] for r in b)}/{len(b)}) | {pct(a):.0f}% ({sum(r['correct'] for r in a)}/{len(a)}) | "
             f"+{pct(a) - pct(b):.0f} pts | ${st.mean(r['cost_usd'] for r in b):.2f} → ${st.mean(r['cost_usd'] for r in a):.2f} | "
             f"{st.mean(r['seconds'] for r in b):.0f} s → {st.mean(r['seconds'] for r in a):.0f} s |")
allb = recs["before", "claude"] + recs["before", "gpt"]; alla = recs["after", "claude"] + recs["after", "gpt"]
spent = sum(json.load(open(f)).get("cost_usd", 0) for f in glob.glob(str(ROOT / "runs/results/*/[cg][lp]*/run*/Q??.json")))
L += ["", f"**Headline:** on the raw public files the two models were right {pct(allb):.0f}% of the time; on the clean database, {pct(alla):.0f}%. "
      f"Wrong answers fell from {sum(not r['correct'] for r in allb)} to {sum(not r['correct'] for r in alla)} out of 138, "
      f"and each answer cost {st.mean(r['cost_usd'] for r in allb) / st.mean(r['cost_usd'] for r in alla):.1f}× less.", "",
      f"Spend: ${spent:.2f} on the 300 graded answers, plus ${failed_cost:.2f} on attempts cut off by API errors and re-run (kept in `runs/failed_attempts/`).", ""]

L += ["## By question (correct out of 3 runs; Q10 and Q23 not scored)", "", "| Q | Claude before | GPT before | Claude after | GPT after |", "|---|---|---|---|---|"]
grid = defaultdict(dict)
for (c, m), rs in recs.items():
    for q in sorted({r["question"] for r in rs}):
        grid[q][c, m] = sum(r["correct"] for r in rs if r["question"] == q)
for q in sorted(grid):
    L.append(f"| {q} | " + " | ".join(str(grid[q].get(k, "-")) for k in [("before", "claude"), ("before", "gpt"), ("after", "claude"), ("after", "gpt")]) + " |")
both_raw = [q for q in sorted(grid) if grid[q][("before", "claude")] == 3 and grid[q][("before", "gpt")] == 3]
L += ["", f"**Right on the raw files too (6/6 before):** {len(both_raw)} of {25 - len(EXCLUDED)} — " + ", ".join(both_raw) + ".", ""]
L += [WRONG]
open(ROOT / "runs/summary.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
