"""
Messy Lab case study page: rebuild runs_index.json from runs/results and put it into page.html.
Writes case_study/index.html. Publish it with the transcripts in runs/results/ served as runs/<condition>/<model>/run<k>/Qxx.json,
plus questions/answer_key.json, questions/grade.py, runs/harness.py and mess_inventory/mess_inventory.csv at the top level.

    python3 case_study/build_page.py
"""
import glob, json
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
EXCLUDED = {"Q10", "Q23"}   # dropped after the run as ambiguous; see runs/summary.md
Q = {q["id"]: q["question"] for q in json.load(open(ROOT / "questions/answer_key.json"))["questions"]}
rows = []
for f in sorted(glob.glob(str(ROOT / "runs/results/*/[cg][lp]*/run*/Q??.json"))):
    d = json.load(open(f))
    a = d.get("answer")
    rows.append({"c": d["condition"], "m": d["model_name"], "r": int(d["run"]), "q": d["question"],
                 "ok": d["correct"] in (True, "True"), "sc": d["question"] not in EXCLUDED,
                 "s": float(d["seconds"]), "usd": float(d["cost_usd"]), "n": int(d["tool_calls"]),
                 "a": json.dumps(a, ensure_ascii=False)[:160] if a is not None else "",
                 "p": f.split("runs/results/")[1]})
idx = json.dumps({"runs": rows, "questions": Q}, separators=(",", ":"))
(HERE / "runs_index.json").write_text(idx)
(HERE / "index.html").write_text((HERE / "page.html").read_text().replace("/*__INDEX__*/", idx.replace("</", "<\\/")))
print(f"case_study/index.html: {len(rows)} answers")
