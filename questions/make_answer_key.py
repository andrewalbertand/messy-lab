"""
Messy Lab, step 4: the answer key.

Combines each question, its hand-checked answer, the query that produced it, the one-line hand check,
and a grading rule into questions/answer_key.json. Every rule is one of three kinds:

  exact_number   the number must match exactly
  exact_list     the set of items must match exactly: no missing items, no extra ones, order ignored.
                 Per-study answers are graded as a list of "key: value" pairs, so every pair must be right.
  tolerance      a number within +/- the stated tolerance

A question can have several parts; all must pass for the question to count as correct. Spellings are
normalized before comparing (case, spaces, units, GSM/GEO IDs, and the aliases listed in `aliases`).

Run from the Messy_Lab folder:   python3 questions/make_answer_key.py
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).parent
STUDIES = ["GSE4025", "GSE21618", "GSE26298", "GSE111151", "GSE117942"]

# Accepted spellings. Keys are the canonical form used in the rules.
ALIASES = {
    "E2": ["e2", "estradiol", "17beta-estradiol", "17β-estradiol", "17-beta-estradiol", "oestradiol"],
    "4-OHT": ["4-oht", "4-hydroxytamoxifen", "4oht", "4-oh tamoxifen", "4-oh-tamoxifen", "4oh-tamoxifen", "hydroxytamoxifen",
              "trans-hydroxytamoxifen", "tot", "afimoxifene", "oht"],
    "tamoxifen (form not stated)": ["tamoxifen", "tamoxifen (form not stated)", "tamoxifen (unspecified)", "tam_unspecified", "tamoxifen, form not stated", "tamoxifen form not stated", "tamoxifen (form unspecified)"],
    "heregulin": ["heregulin", "hrg", "heregulin-beta1", "heregulin beta1", "neuregulin", "hrg-b1", "hrg_b1"],
    "fulvestrant": ["fulvestrant", "ici 182,780", "faslodex"],
    "GDC-0810": ["gdc-0810", "gdc0810", "brilanestrant"],
    "GDC-0927": ["gdc-0927", "gdc0927"],
    "GNE-274": ["gne-274", "gne274", "g-03046274"],
    "ethanol": ["ethanol", "etoh", "vehicle (ethanol)"],
    "DMSO": ["dmso", "dimethyl sulfoxide"],
    "not stated": ["not stated", "unknown", "not given", "none", "not 4-oht", "not 4-oht / not stated",
                   "not 4-oht: tamoxifen of unstated form, dose not stated", "no 4-oht", "unspecified"],
    "chronic": ["chronic", "chronic (months)", "months", "8-12 months", "eight to twelve months"],
    "yes": ["yes", "true", "y"],
    "no": ["no", "false", "n"],
    "ZR-75-1 Tam2": ["zr-75-1 tam2", "zr751 tam2", "zr-75-1_tam2", "zr-75-1 tam2 (helsinki)"],
}

ALL_LINES = ["MCF-7", "CAMA-1", "HCC1500", "BT-474", "T-47D", "EFM-19", "MDA-MB-330", "ZR-75-1"]

# qid: (parts, why this rule)
RULES = {
    "Q01": ([{"field": "per_study", "type": "exact_list",
              "expected": ["GSE4025: 6", "GSE26298: 4", "GSE111151: 1", "GSE117942: 2"]},
             {"field": "total", "type": "exact_number", "expected": 13}],
            "Counts are exact; per-study pairs catch answers that reach 13 by the wrong route. GSE21618 must not appear."),
    "Q02": ([{"field": "per_study", "type": "exact_list",
              "expected": ["GSE4025: 5", "GSE21618: 31", "GSE26298: 8", "GSE111151: 4", "GSE117942: 14"]}],
            "Exact counts per study; GSE21618 = 31 requires treating 0 h samples as untreated."),
    "Q03": ([{"field": "per_line", "type": "exact_list",
              "expected": ["MCF-7: 188", "CAMA-1: 14", "HCC1500: 14", "BT-474: 13", "T-47D: 13", "EFM-19: 10",
                           "MDA-MB-330: 10", "ZR-75-1: 3"]}],
            "Exact counts per line, resistant derivatives counted under their parent."),
    "Q04": ([{"field": "per_compound", "type": "exact_list",
              "expected": ["E2: 76", "tamoxifen (form not stated): 49", "heregulin: 42", "4-OHT: 31", "GDC-0810: 14",
                           "DMSO: 14", "fulvestrant: 14", "ethanol: 13", "GDC-0927: 6", "GNE-274: 6"]}],
            "Exact counts per compound; the spellings above are all accepted for each compound."),
    "Q05": ([{"field": "total", "type": "exact_number", "expected": 13}],
            "6 RIKEN clones + 7 Helsinki lines; RIKEN's time-course clone is one of its six, so 14 is wrong."),
    "Q06": ([{"field": "genes", "type": "tolerance", "expected": 11781, "tolerance": 236}],
            "±2% (11,545-12,017): a proper gene-ID mapping by another route gives 11,703; matching raw names gives "
            "10,687 and fails."),
    "Q07": ([{"field": "samples", "type": "exact_list", "expected": ["GSM645718", "GSM645719"]}],
            "Only the same-experiment controls; MD31/MD32 (GSM645710/11) are the other experiment's and fail."),
    "Q08": ([{"field": "count", "type": "exact_number", "expected": 31}],
            "15 labeled untreated + 16 labeled with a drug at 0 h; answering 15 fails."),
    "Q09": ([{"field": "per_study", "type": "exact_list",
              "expected": ["GSE4025: 10", "GSE26298: 100", "GSE26298: 500", "GSE111151: 1000", "GSE117942: 1000",
                           "GSE21618: not stated"]}],
            "Doses in nM, every dose per study; GSE21618 must be marked not stated / not 4-OHT, never given a number."),
    "Q10": ([{"field": "per_study", "type": "exact_list",
              "expected": ["GSE4025: 24", "GSE117942: 24", "GSE111151: chronic", "GSE26298: not stated"]}],
            "Hours as a number; GSE26298 must be 'not stated' (guessing 72 fails); GSE21618 must not appear."),
    "Q11": ([{"field": "samples", "type": "exact_list",
              "expected": [f"GSM33156{n}" for n in range(88, 96)]}],
            "Exactly the 8 GSE117942 samples in the run with no DMSO."),
    "Q12": ([{"field": "per_study", "type": "exact_list",
              "expected": ["GSE4025: yes", "GSE21618: yes", "GSE26298: yes", "GSE117942: yes", "GSE111151: no"]}],
            "Each study marked deprived (yes) or not (no)."),
    "Q13": ([{"field": "samples", "type": "exact_list", "expected": ["GSM539725", "GSM539726", "GSM539727", "GSM539728"]}],
            "Exactly the four RIKEN samples whose source name disagrees with their title and treatment field."),
    "Q14": ([{"field": "study", "type": "exact_list", "expected": ["GSE117942"]},
             {"field": "corrected_nM", "type": "exact_number", "expected": 1000}],
            "The study and the corrected dose (1 uM = 1,000 nM)."),
    "Q15": ([{"field": "studies", "type": "exact_list", "expected": ["GSE4025", "GSE117942"]}],
            "Only mismatches with GEO's own processing text; GSE26298's mismatch is with the paper, so including it fails."),
    "Q16": ([{"field": "count", "type": "exact_number", "expected": 28},
             {"field": "study", "type": "exact_list", "expected": ["GSE117942"]}],
            "Exact count and the one study."),
    "Q17": ([{"field": "samples", "type": "exact_list", "expected": ["GSM645720", "GSM645721"]}],
            "The two RAR-alpha samples both labeled rep1."),
    "Q18": ([{"field": "answer", "type": "exact_list", "expected": ["no"]},
             {"field": "log2_change.GSE4025", "type": "tolerance", "expected": -0.07, "tolerance": 0.3},
             {"field": "log2_change.GSE117942", "type": "tolerance", "expected": -2.52, "tolerance": 0.3}],
            "No, plus each study's log2 change within ±0.3 (room for probe choice or normalization)."),
    "Q19": ([{"field": "compounds", "type": "exact_list", "expected": ["E2"]}],
            "Only E2 rises in all seven lines against same-run DMSO."),
    "Q20": ([{"field": "answer", "type": "exact_list", "expected": ["yes"]}],
            "Yes: positive in every study and run that tested 24 h E2."),
    "Q21": ([{"field": "log2_change", "type": "tolerance", "expected": -2.82, "tolerance": 0.3}],
            "±0.3 log2."),
    "Q22": ([{"field": "line", "type": "exact_list", "expected": ["ZR-75-1 Tam2"]}],
            "One line."),
    "Q23": ([{"field": "studies", "type": "exact_list", "expected": ["GSE21618", "GSE26298", "GSE111151"]}],
            "The three studies that can't be pooled for 4-OHT response; reasons are recorded but not graded."),
    "Q24": ([{"field": "studies", "type": "exact_list", "expected": ["GSE4025", "GSE117942"]}],
            "Exactly the two studies with ordinary MCF-7, 24 h 4-OHT and a same-run control."),
    "Q25": ([{"field": "log2_change", "type": "tolerance", "expected": -1.29, "tolerance": 0.3}],
            "±0.3 log2; pooling other studies or skipping same-run controls moves it outside."),
}


def main():
    ans = {a["id"]: a for a in json.load(open(HERE / "answers.json"))["answers"]}
    checks = {}
    for line in (HERE / "hand_checks.md").read_text().splitlines():
        m = re.match(r"^\| (Q\d\d) \|.*\| (✅|❌) \| (.*) \|$", line)
        if m:
            checks[m.group(1)] = (m.group(2) == "✅", m.group(3))
    key = []
    for qid in sorted(RULES):
        parts, why = RULES[qid]
        ok, how = checks[qid]
        assert ok, f"{qid} failed its hand check"
        key.append({"id": qid, "question": ans[qid]["question"], "answer": ans[qid]["answer"],
                    "grading": {"parts": parts, "rule": why}, "hand_check": how, "query": ans[qid]["query"]})
    out = {"version": "2026-09-27", "rule_types": {
        "exact_number": "the number must match exactly",
        "exact_list": "the set of items must match exactly (no missing, no extra, order ignored)",
        "tolerance": "a number within ± the stated tolerance"},
        "aliases": ALIASES, "questions": key}
    (HERE / "answer_key.json").write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str))
    kinds = {}
    for q in key:
        for p in q["grading"]["parts"]:
            kinds[p["type"]] = kinds.get(p["type"], 0) + 1
    print(f"answer_key.json: {len(key)} questions, {sum(kinds.values())} graded parts {kinds}")


if __name__ == "__main__":
    main()
