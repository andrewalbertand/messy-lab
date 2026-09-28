"""
Messy Lab: grade one model answer against the answer key.

    from grade import grade
    grade("Q07", {"samples": ["GSM645718", "GSM645719"]})   ->  (True, [part results])

A model answer is a JSON object with the fields named in each question's grading parts (e.g.
{"per_study": {"GSE4025": 6, ...}, "total": 13}). Per-study fields may be given as an object or as a
list of "key: value" strings. Running this file tests the grader: the key's own answers must all pass,
and known wrong answers (the traps) must fail.
"""
import json
import re
from pathlib import Path

KEY = json.load(open(Path(__file__).parent / "answer_key.json", encoding="utf-8"))
Q = {q["id"]: q for q in KEY["questions"]}
ALIAS = {a.lower(): canon for canon, al in KEY["aliases"].items() for a in al + [canon]}


def norm(x):
    """Canonical spelling for one item: trims, lowercases, maps aliases, strips units and '.0'."""
    s = str(x).strip()
    s = re.sub(r"\s*(nm|nanomolar|h|hr|hrs|hours)\s*$", "", s, flags=re.I)
    s = re.sub(r"\.0+$", "", s)
    low = s.lower()
    if low in ALIAS:
        return ALIAS[low]
    m = re.fullmatch(r"(.+?)\s*\((.+)\)", low)      # "17beta-estradiol (E2)" -> try each half
    if m:
        for part in (m.group(1), m.group(2)):
            if part.strip() in ALIAS:
                return ALIAS[part.strip()]
    if re.fullmatch(r"(gsm|gse)\d+", low):
        return low.upper()
    try:
        f = float(s.replace(",", ""))
        return str(int(f)) if f == int(f) else str(round(f, 4))
    except ValueError:
        return low


def as_items(v):
    """Turn a dict, list or scalar into a set of normalized items; dict values that are lists expand to pairs."""
    items = set()
    if isinstance(v, dict):
        for k, val in v.items():
            for one in (val if isinstance(val, list) else [val]):
                items.add(f"{norm(k)}: {norm(one)}")
    elif isinstance(v, list):
        for one in v:
            if isinstance(one, str) and ":" in one and not one.lower().startswith("not 4-oht:"):
                k, val = one.split(":", 1)
                items.add(f"{norm(k)}: {norm(val)}")
            else:
                items.add(norm(one))
    else:
        items.add(norm(v))
    return items


def expected_items(exp):
    out = set()
    for e in exp:
        if isinstance(e, str) and ": " in e:
            k, v = e.split(": ", 1)
            out.add(f"{norm(k)}: {norm(v)}")
        else:
            out.add(norm(e))
    return out


def grade(qid, response):
    results = []
    for part in Q[qid]["grading"]["parts"]:
        got = response
        for key in part["field"].split("."):          # "log2_change.GSE4025" reads a nested value
            got = got.get(key) if isinstance(got, dict) else None
        if got is None:
            results.append((part["field"], False, "missing"))
            continue
        try:
            if part["type"] == "exact_number":
                ok = float(str(got).replace(",", "")) == float(part["expected"])
            elif part["type"] == "tolerance":
                ok = abs(float(str(got).replace(",", "")) - float(part["expected"])) <= float(part["tolerance"])
            else:
                exp = expected_items(part["expected"])
                items = {i for i in as_items(got) if not (i.endswith(": 0") and i not in exp)}   # "X: 0" = X not counted
                ok = items == exp
        except (TypeError, ValueError):
            ok = False
        results.append((part["field"], ok, got))
    return all(ok for _, ok, _ in results), results


if __name__ == "__main__":
    right = {
        "Q01": {"per_study": {"GSE4025": 6, "GSE26298": 4, "GSE111151": 1, "GSE117942": 2}, "total": 13},
        "Q02": {"per_study": {"GSE4025": 5, "GSE21618": 31, "GSE26298": 8, "GSE111151": 4, "GSE117942": 14}},
        "Q03": {"per_line": {"MCF-7": 188, "CAMA-1": 14, "HCC1500": 14, "BT-474": 13, "T-47D": 13, "EFM-19": 10, "MDA-MB-330": 10, "ZR-75-1": 3}},
        "Q04": {"per_compound": {"17beta-estradiol": 76, "tamoxifen (form not stated)": 49, "heregulin-beta1": 42, "4-hydroxytamoxifen": 31,
                                 "brilanestrant": 14, "dimethyl sulfoxide": 14, "fulvestrant": 14, "ethanol": 13, "GDC-0927": 6, "GNE-274": 6}},
        "Q05": {"total": 13}, "Q06": {"genes": 11781}, "Q07": {"samples": ["GSM645718", "GSM645719"]}, "Q08": {"count": 31},
        "Q09": {"per_study": {"GSE4025": 10, "GSE26298": [100, 500], "GSE111151": "1000 nM", "GSE117942": 1000, "GSE21618": "not stated"}},
        "Q10": {"per_study": {"GSE4025": "24 h", "GSE117942": 24, "GSE111151": "chronic (months)", "GSE26298": "not stated"}},
        "Q11": {"samples": [f"GSM33156{n}" for n in range(88, 96)]},
        "Q12": {"per_study": {"GSE4025": "yes", "GSE21618": "yes", "GSE26298": "yes", "GSE117942": "yes", "GSE111151": "no"}},
        "Q13": {"samples": ["GSM539725", "GSM539726", "GSM539727", "GSM539728"]},
        "Q14": {"study": ["GSE117942"], "corrected_nM": 1000}, "Q15": {"studies": ["GSE4025", "GSE117942"]},
        "Q16": {"count": 28, "study": "GSE117942"}, "Q17": {"samples": ["GSM645720", "GSM645721"]},
        "Q18": {"answer": "no", "log2_change": {"GSE4025": -0.07, "GSE117942": -2.52}}, "Q19": {"compounds": ["estradiol"]},
        "Q20": {"answer": "yes"}, "Q21": {"log2_change": -2.82}, "Q22": {"line": "ZR-75-1 Tam2"},
        "Q23": {"studies": ["GSE21618", "GSE26298", "GSE111151"]}, "Q24": {"studies": ["GSE117942", "GSE4025"]},
        "Q25": {"log2_change": -1.29},
    }
    # forms real models used in the full run (27 Sep): zero entries, names with codes in brackets
    right_forms = {
        "Q01": {"per_study": {"GSE4025": 6, "GSE21618": 0, "GSE26298": 4, "GSE111151": 1, "GSE117942": 2}, "total": 13},
        "Q04": {"per_compound": {"17beta-estradiol (E2)": 76, "tamoxifen, form not stated (TAM_UNSPECIFIED)": 49, "heregulin-beta1 (HRG_B1)": 42,
                                 "4-hydroxytamoxifen (4-OHT)": 31, "fulvestrant": 14, "GDC-0810 (brilanestrant)": 14, "DMSO": 14,
                                 "ethanol": 13, "GDC-0927": 6, "GNE-274": 6}},
    }
    for q, a in right_forms.items():
        if not grade(q, a)[0]:
            print("real answer form FAILED:", q, grade(q, a)[1])
    traps = {
        "Q01": {"per_study": {"GSE4025": 6, "GSE26298": 4, "GSE111151": 1, "GSE117942": 2, "GSE21618": 49}, "total": 62},
        "Q02": {"per_study": {"GSE4025": 5, "GSE21618": 15, "GSE26298": 8, "GSE111151": 4, "GSE117942": 14}},
        "Q05": {"total": 14}, "Q06": {"genes": 10687}, "Q07": {"samples": ["GSM645710", "GSM645711"]}, "Q08": {"count": 15},
        "Q09": {"per_study": {"GSE4025": 10, "GSE26298": 100, "GSE111151": 1000, "GSE117942": 1000, "GSE21618": 1000}},
        "Q10": {"per_study": {"GSE4025": 24, "GSE117942": 24, "GSE111151": "chronic", "GSE26298": 72}},
        "Q11": {"samples": []}, "Q15": {"studies": ["GSE4025", "GSE117942", "GSE26298"]},
        "Q18": {"answer": "yes", "log2_change": {"GSE4025": -0.07, "GSE117942": -2.52}},
        "Q19": {"compounds": ["E2", "GDC-0927"]}, "Q24": {"studies": ["GSE4025", "GSE117942", "GSE26298"]},
        "Q25": {"log2_change": -2.52},
        "Q04": {"per_compound": {"17beta-estradiol (E2)": 76, "tamoxifen, form not stated (TAM_UNSPECIFIED)": 80, "heregulin-beta1 (HRG_B1)": 42,
                                 "fulvestrant": 14, "GDC-0810 (brilanestrant)": 14, "DMSO": 14, "ethanol": 13, "GDC-0927": 6, "GNE-274": 6}},
    }
    passed = sum(grade(q, a)[0] for q, a in right.items())
    caught = sum(not grade(q, a)[0] for q, a in traps.items())
    for q, a in right.items():
        if not grade(q, a)[0]:
            print("key answer FAILED:", q, grade(q, a)[1])
    for q, a in traps.items():
        if grade(q, a)[0]:
            print("trap answer PASSED (bad):", q)
    print(f"Right answers graded correct: {passed}/{len(right)}. Known wrong answers graded wrong: {caught}/{len(traps)}.")
