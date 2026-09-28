"""
Messy Lab assistant: answers questions about the five studies by querying db/messy_lab.duckdb,
cites the sample, file and issue records it used, and says so when the data can't answer.

    python3 assistant/assistant.py "How many MCF-7 samples got 4-OHT?"      # one question in the terminal
    python3 assistant/app.py                                                  # the browser version

Needs ANTHROPIC_API_KEY in .env (same file as runs/harness.py). The database is opened read-only.
"""
import json
import os
import re
import sys
import time
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parent.parent
DB = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db")) / "messy_lab.duckdb"
MODEL = os.environ.get("MESSY_LAB_ASSISTANT_MODEL", "claude-opus-5-5")
MAX_STEPS, ROW_CAP = 12, 60


def load_dotenv(path=ROOT / ".env"):
    if path.exists():
        for line in path.read_text().splitlines():
            if "=" in line and not line.lstrip().startswith("#"):
                k, v = line.split("=", 1)
                os.environ[k.strip()] = v.strip().strip('"').strip("'")


SYSTEM = f"""You are the data assistant for a research program on tamoxifen response in ER-positive breast cancer cells.
The program's data is one clean DuckDB database built from five public GEO studies. You answer questions from
scientists and managers by querying that database with the run_sql tool, then calling final_answer.

How to work:
- Answer only from what your queries return. Never use outside knowledge for facts about these samples.
- Keep queries small and specific. Check the issue table when a result depends on a dose, time, control or label.
- Cite what you used: every sample ID behind a count or a number (up to 40; if more, cite a representative set and
  say how many in total), the source_file IDs the values came from, and the inventory IDs (M01-M75) of any issue
  that affects the answer.
- If the database records the value as unknown, NULL, 'not stated' or flagged, or the question asks for something
  the data doesn't contain, set can_answer to false. Say plainly what is missing and why, citing the issue rows.
  Do not estimate, infer from other studies or fill gaps. A clear "the data can't answer this" is a correct answer.
- Write the answer for a busy scientist: the direct answer first, in one or two sentences, then at most three short
  sentences of support. Use study GEO accessions (GSE...) and sample GSM IDs, not internal study IDs, in the text.

The data dictionary follows. Its ground rules are binding.

{(ROOT / "db" / "DATA_DICTIONARY.md").read_text()}
"""

TOOLS = [
    {"name": "run_sql",
     "description": f"Run one read-only DuckDB SQL query (SELECT or WITH) against the program database. Returns up to {ROW_CAP} rows.",
     "input_schema": {"type": "object", "properties": {
         "purpose": {"type": "string", "description": "One short line, in plain English, saying what this query checks."},
         "sql": {"type": "string"}}, "required": ["purpose", "sql"]}},
    {"name": "final_answer",
     "description": "Give the final answer with its citations. Call exactly once, at the end.",
     "input_schema": {"type": "object", "properties": {
         "can_answer": {"type": "boolean", "description": "false if the database can't answer the question as asked"},
         "answer": {"type": "string"},
         "sample_ids": {"type": "array", "items": {"type": "string"}, "description": "GSM IDs of the samples used"},
         "source_file_ids": {"type": "array", "items": {"type": "string"}, "description": "source_file IDs the values came from"},
         "inventory_ids": {"type": "array", "items": {"type": "string"}, "description": "M01-M75 rows of issues that affect the answer"},
         "total_samples": {"type": "integer", "description": "How many samples the answer rests on, if more were used than cited"}},
         "required": ["can_answer", "answer", "sample_ids", "source_file_ids", "inventory_ids"]}},
]

WRITE = re.compile(r"\b(insert|update|delete|drop|create|alter|attach|detach|copy|export|import|install|load|pragma|set|call)\b", re.I)


def connect():
    return duckdb.connect(str(DB), read_only=True)


def run_sql(con, sql):
    s = sql.strip().rstrip(";")
    if not re.match(r"(?is)^\s*(select|with)\b", s) or ";" in s or WRITE.search(re.sub(r"'[^']*'", "''", s)):
        return {"error": "Only one read-only SELECT or WITH query is allowed."}
    try:
        cur = con.execute(s)
        cols = [d[0] for d in cur.description]
        rows = cur.fetchmany(ROW_CAP + 1)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {str(e)[:400]}"}
    return {"columns": cols, "rows": [[None if v is None else (round(v, 4) if isinstance(v, float) else str(v) if not isinstance(v, (int, bool)) else v) for v in r] for r in rows[:ROW_CAP]],
            "truncated": len(rows) > ROW_CAP}


def as_text(res):
    if "error" in res:
        return "ERROR " + res["error"]
    lines = [" | ".join(res["columns"])] + [" | ".join("NULL" if v is None else str(v) for v in r) for r in res["rows"]]
    return "\n".join(lines) + (f"\n(first {ROW_CAP} rows only)" if res["truncated"] else "") + f"\n({len(res['rows'])} rows)"


def resolve(con, fa):
    """Check every citation against the database and attach the records behind it."""
    ids = [s for s in dict.fromkeys(x.strip() for x in fa.get("sample_ids", [])) if s]
    samples = con.execute("""
        SELECT s.sample_id, st.geo_accession, s.geo_title, s.role, s.batch, s.source_file_id,
               string_agg(DISTINCT t.compound_id || COALESCE(' ' || CAST(t.dose_nM AS VARCHAR) || ' nM', '') ||
                          COALESCE(' ' || CAST(t.exposure_h AS INTEGER) || ' h', ''), ', ') AS treatment
        FROM sample s JOIN study st USING (study_id) LEFT JOIN sample_treatment t USING (sample_id)
        WHERE s.sample_id IN (SELECT unnest(?)) GROUP BY ALL ORDER BY s.sample_id""", [ids]).fetchall() if ids else []
    found = {r[0] for r in samples}
    prov = con.execute("""
        SELECT entity_id, field, source_file_id, source_location, raw_text, rule FROM provenance
        WHERE split_part(entity_id, '|', 1) IN (SELECT unnest(?))
          AND field IN ('role', 'compound_id', 'dose_nM', 'exposure_h', 'batch')
        ORDER BY entity_id, field""", [ids]).fetchall() if ids else []
    files = set(fa.get("source_file_ids", [])) | {r[5] for r in samples}
    frows = con.execute("SELECT source_file_id, path, sha256, origin_url FROM source_file WHERE source_file_id IN (SELECT unnest(?)) ORDER BY 1",
                        [sorted(files)]).fetchall() if files else []
    inv = [i for i in dict.fromkeys(x.strip().upper() for x in fa.get("inventory_ids", [])) if re.fullmatch(r"M\d\d", i)]
    irows = con.execute("""SELECT inventory_id, any_value(impact), any_value(message), count(*), any_value(resolution)
                           FROM issue WHERE inventory_id IN (SELECT unnest(?)) GROUP BY 1 ORDER BY 1""", [inv]).fetchall() if inv else []
    return {
        "samples": [dict(zip(["sample_id", "study", "title", "role", "batch", "source_file_id", "treatment"], r)) for r in samples],
        "not_found": [s for s in ids if s not in found],
        "provenance": [dict(zip(["entity_id", "field", "source_file_id", "source_location", "raw_text", "rule"], r)) for r in prov],
        "files": [dict(zip(["source_file_id", "path", "sha256", "origin_url"], r)) for r in frows],
        "issues": [dict(zip(["inventory_id", "impact", "message", "records", "resolution"], r)) for r in irows],
    }


def ask(question, on_event=lambda kind, data: None):
    """Run one question to completion. on_event(kind, data) streams progress: 'step', 'result', 'answer', 'error'."""
    import anthropic
    load_dotenv()
    client = anthropic.Anthropic()
    con = connect()
    msgs = [{"role": "user", "content": question}]
    t0, steps, usage = time.time(), [], {"input": 0, "output": 0}
    for _ in range(MAX_STEPS):
        r = client.messages.create(model=MODEL, max_tokens=4000, system=SYSTEM, tools=TOOLS, messages=msgs)
        usage["input"] += r.usage.input_tokens; usage["output"] += r.usage.output_tokens
        msgs.append({"role": "assistant", "content": [b.model_dump() for b in r.content]})
        results = []
        for b in r.content:
            if b.type != "tool_use":
                continue
            if b.name == "final_answer":
                out = {"question": question, **b.input, "citations": resolve(con, b.input), "steps": steps,
                       "seconds": round(time.time() - t0, 1), "usage": usage, "model": MODEL}
                on_event("answer", out)
                return out
            res = run_sql(con, b.input.get("sql", ""))
            step = {"purpose": b.input.get("purpose", ""), "sql": b.input.get("sql", ""),
                    "rows": len(res.get("rows", [])), "error": res.get("error")}
            steps.append(step); on_event("step", step)
            results.append({"type": "tool_result", "tool_use_id": b.id, "content": as_text(res)[:12000]})
        if not results:
            msgs.append({"role": "user", "content": "Call final_answer now with your answer and citations."})
        else:
            msgs.append({"role": "user", "content": results})
    err = {"question": question, "error": f"No answer after {MAX_STEPS} steps."}
    on_event("error", err)
    return err


if __name__ == "__main__":
    q = " ".join(sys.argv[1:]) or "How many MCF-7 samples were given 4-OHT, and in which studies?"
    out = ask(q, lambda k, d: print(f"  · {d['purpose']} ({d['rows']} rows)") if k == "step" else None)
    print()
    print(("ANSWER: " if out.get("can_answer") else "CAN'T ANSWER FROM THE DATA: ") + out.get("answer", out.get("error", "")))
    c = out.get("citations", {})
    if c:
        print("Samples:", ", ".join(s["sample_id"] for s in c["samples"]) or "none")
        print("Files:  ", ", ".join(f"{f['source_file_id']} ({f['path']})" for f in c["files"]) or "none")
        print("Issues: ", ", ".join(f"{i['inventory_id']} ({i['impact']})" for i in c["issues"]) or "none")
        if c["not_found"]:
            print("Cited but not in the database (dropped):", ", ".join(c["not_found"]))
    print(f"({out.get('seconds')} s, {len(out.get('steps', []))} queries)")
