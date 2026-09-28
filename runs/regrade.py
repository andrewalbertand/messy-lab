"""Re-grade every saved answer in runs/results/ with the current questions/grade.py (no API calls)."""
import glob, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "questions"))
from grade import grade

changed = 0
for f in sorted(glob.glob(str(ROOT / "runs/results/*/*/run*/Q??.json"))):
    d = json.load(open(f))
    ok, parts = grade(d["question"], d["answer"] or {})
    ok = bool(ok) and not d.get("outside_access_attempt")   # same rule as harness.py
    parts = [{"field": a, "ok": b, "got": c} for a, b, c in parts]
    if ok != d.get("correct"):
        changed += 1
        print(f"{Path(f).relative_to(ROOT)}: {d.get('correct')} -> {ok}")
    d["correct"], d["grade_parts"] = ok, parts
    json.dump(d, open(f, "w"), indent=2)
print(f"{changed} grades changed")
