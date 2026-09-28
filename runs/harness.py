"""
Messy Lab, step 5: the before/after harness.

Asks each model each of the 25 questions, once with the raw files ("before") and once with the clean
database ("after"), and grades every answer with questions/grade.py. Everything that could tilt the
result is held the same for every model and both conditions:

  * the same instructions (only the "your data" paragraph differs between before and after)
  * the same two tools: run_python (code runs in a fresh folder holding only that condition's files;
    pandas, numpy and duckdb available; no internet) and submit_answer
  * the same limits: TOOL_CALLS per question, SECONDS per question, per-call timeout and output cap
  * a fresh conversation per question, so one question can't leak into another

  before  the working folder holds a copy of studies/ (raw GEO downloads, papers, gene reference)
  after   the working folder holds a copy of db/messy_lab.duckdb and db/DATA_DICTIONARY.md

Each run is saved to runs/results/<condition>/<model>/run<k>/<Qxx>.json with the full transcript,
every piece of code the model ran and its output, the submitted answer, the grade, tokens and time.
Existing results are skipped, so an interrupted batch picks up where it stopped.

Usage (from the Messy_Lab folder):
    python3 runs/harness.py --models mock --runs 1 --questions Q07,Q08          # offline self-test
    python3 runs/harness.py --models claude --conditions before --runs 1 --questions Q07   # one real call
    python3 runs/harness.py --models claude,gpt,gemini --conditions before,after --runs 3 --workers 4
"""
import argparse
import site
import concurrent.futures as cf
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "questions"))
from grade import grade  # noqa: E402

def load_dotenv(path=ROOT / ".env"):
    """Read KEY=value lines from .env into the environment (existing variables win)."""
    if path.exists():
        for line in path.read_text().splitlines():
            if "=" in line and not line.lstrip().startswith("#"):
                k, v = line.split("=", 1)
                os.environ[k.strip()] = v.strip().strip('"').strip("'")
    # .env wins over the shell; drop stray endpoint overrides not set in .env
    for var in ("ANTHROPIC_BASE_URL", "OPENAI_BASE_URL", "OPENAI_API_BASE"):
        if path.exists() and var + "=" not in path.read_text():
            os.environ.pop(var, None)
    for var in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY"):
        v = os.environ.get(var, "")
        print(f"[keys] {var}: {'…' + v[-4:] if v else 'MISSING'}", flush=True)


load_dotenv()
RESULTS = Path(os.environ.get("MESSY_LAB_RESULTS", ROOT / "runs" / "results"))
DB_DIR = Path(os.environ.get("MESSY_LAB_DB_DIR", ROOT / "db"))

# ---- limits: identical for every model and condition ------------------------------------------
TOOL_CALLS = 30          # run_python calls allowed per question
SECONDS = 900            # wall-clock budget per question
CALL_TIMEOUT = 180       # seconds per run_python call
OUTPUT_CAP = 20_000      # characters of output returned to the model per call
MAX_TOKENS = 16_000      # per model reply (room for reasoning tokens)

SYSTEM = """You are a data analyst at a biotech company. A colleague has asked you a question about the \
company's gene-expression data. Answer it from the files in your working folder.

You have two tools.
- run_python: runs Python 3 code in your working folder and returns what it prints. pandas, numpy and duckdb \
and pypdf are installed. There is no internet. Nothing is remembered between calls, so re-load what you need each time.
- submit_answer: call this exactly once, when you are done, with your answer as a JSON object in the format \
the question asks for.

Work only with the files in your working folder. If the data can't settle part of the question, say so in \
the answer (for example "not stated") rather than guessing."""

CONDITION_TEXT = {
    "before": "Your data: the folder `studies/` holds everything the team has on this program. `studies/data/raw` "
              "has the raw downloads of five GEO studies (series matrix files, platform annotation files, a tar of "
              "per-sample files and a counts table). `studies/papers` has the five published papers as PDFs "
              "(paper_1.pdf ... paper_5.pdf, in the same order as the raw files). `studies/data/reference` has the "
              "HGNC gene list.",
    "after": "Your data: `messy_lab.duckdb` is a DuckDB database holding the same five studies, cleaned into one "
             "system, and `DATA_DICTIONARY.md` explains every table and column. Read the dictionary first.",
}

# The JSON fields each answer must use. These describe the shape only, never the answer.
FIELDS = {
    "Q01": '{"per_study": {GEO ID: number of samples}, "total": number}',
    "Q02": '{"per_study": {GEO ID: number of control samples}}',
    "Q03": '{"per_line": {cell line name: number of samples}}',
    "Q04": '{"per_compound": {compound name: number of samples}}',
    "Q05": '{"total": number, "per_study": {GEO ID: [line names]}}',
    "Q06": '{"genes": number}',
    "Q07": '{"samples": [GSM IDs]}',
    "Q08": '{"count": number}',
    "Q09": '{"per_study": {GEO ID: dose in nM (a list if several), or "not stated"}}',
    "Q10": '{"per_study": {GEO ID: hours as a number, or "chronic", or "not stated"}}',
    "Q11": '{"samples": [GSM IDs]}',
    "Q12": '{"per_study": {GEO ID: "yes" if deprived, "no" if not}}',
    "Q13": '{"samples": [GSM IDs]}',
    "Q14": '{"study": [GEO IDs, empty if none], "dose_as_written": text, "corrected_nM": number}',
    "Q15": '{"studies": [GEO IDs]}',
    "Q16": '{"count": number, "study": [GEO IDs]}',
    "Q17": '{"samples": [GSM IDs]}',
    "Q18": '{"answer": "yes" or "no", "log2_change": {GEO ID: number}}',
    "Q19": '{"compounds": [compound names]}',
    "Q20": '{"answer": "yes" or "no", "log2_change": {GEO ID or GEO ID + run: number}}',
    "Q21": '{"log2_change": number}',
    "Q22": '{"line": line name}',
    "Q23": '{"studies": [GEO IDs], "reasons": {GEO ID: text}}',
    "Q24": '{"studies": [GEO IDs]}',
    "Q25": '{"log2_change": number, "per_study": {GEO ID: number}}',
}

TOOLS = [
    {"name": "run_python",
     "description": "Run Python 3 code in the working folder and return its printed output (stdout and stderr). "
                    f"Each call starts fresh. Time limit {CALL_TIMEOUT} s. Output over {OUTPUT_CAP:,} characters is cut.",
     "parameters": {"type": "object", "properties": {"code": {"type": "string", "description": "Python code to run"}},
                    "required": ["code"]}},
    {"name": "submit_answer",
     "description": "Submit your final answer. Call once, at the end.",
     "parameters": {"type": "object", "properties": {
         "answer_json": {"type": "string", "description": "The answer as a JSON object, in the format the question asks for"},
         "explanation": {"type": "string", "description": "One or two sentences on how you got it"}},
                    "required": ["answer_json", "explanation"]}},
]


# ---- the question text -----------------------------------------------------------------------
def load_questions():
    md = (ROOT / "questions" / "questions.md").read_text()
    intro = md.split("## A.")[0]
    preamble = intro.split("---")[0].split("sees in both runs:")[-1] if "sees in both runs:" in intro else ""
    table = "\n".join(l for l in intro.splitlines() if l.startswith("|"))
    defs = intro[intro.index('"4-OHT" means'):].split("---")[0].strip()
    qs, cur = {}, None
    for line in md.splitlines():
        m = re.match(r"^\*\*(Q\d\d)\.\*\* (?:\*\([^)]*\)\* )?(.*)", line)
        if m:
            cur = m.group(1)
            qs[cur] = m.group(2)
        elif cur and line.startswith("**Answer format:**"):
            qs[cur] += "\n" + line.replace("**", "")
            cur = None
        elif cur and line.strip():
            qs[cur] += " " + line.strip()
    context = "The five studies, by GEO ID:\n" + table + "\n\n" + defs
    return context, qs


# ---- the sandbox -----------------------------------------------------------------------------
def make_workspace(condition):
    ws = Path(tempfile.mkdtemp(prefix=f"messylab_{condition}_"))
    if condition == "before":
        shutil.copytree(ROOT / "studies", ws / "studies", ignore=shutil.ignore_patterns(".DS_Store"))
    else:
        shutil.copy2(DB_DIR / "messy_lab.duckdb", ws / "messy_lab.duckdb")
        shutil.copy2(ROOT / "db" / "DATA_DICTIONARY.md", ws / "DATA_DICTIONARY.md")
    return ws


FORBIDDEN = re.compile(r"answer_key|answers\.json|hand_check|mess_inventory|Messy_Lab|messy_lab/|/pipeline/|\.\./")


def run_python(code, ws):
    """Run code in the workspace with a clean environment. Flags any attempt to reach outside it."""
    flagged = bool(FORBIDDEN.search(code))
    env = {"PATH": os.environ.get("PATH", ""), "HOME": str(ws), "PYTHONDONTWRITEBYTECODE": "1",
           "MPLBACKEND": "Agg", "TMPDIR": str(ws), "PYTHONUSERBASE": site.getuserbase()}
    try:
        p = subprocess.run([sys.executable, "-c", code], cwd=ws, env=env, capture_output=True, text=True,
                           timeout=CALL_TIMEOUT)
        out = (p.stdout + ("\n[stderr]\n" + p.stderr if p.stderr.strip() else "")).strip()
    except subprocess.TimeoutExpired:
        out = f"[timed out after {CALL_TIMEOUT} s]"
    if len(out) > OUTPUT_CAP:
        out = out[:OUTPUT_CAP] + f"\n[output cut at {OUTPUT_CAP:,} characters]"
    return out or "[no output]", flagged


# ---- model adapters: one conversation, same tools, same limits -----------------------------------
class Model:
    """Common shape: start(system, user) -> step() returns (text, [tool calls]); add_results([...])."""
    usage = None


class Anthropic(Model):
    def nudge(self, text):
        self.msgs.append({"role": "user", "content": text})

    def __init__(self, model_id, key):
        import anthropic
        ws = os.environ.get("ANTHROPIC_WORKSPACE_ID")   # needed when the key isn't scoped to a workspace
        self.c = anthropic.Anthropic(api_key=key, default_headers={"anthropic-workspace-id": ws} if ws else None)
        self.m = model_id
        self.tools = [{"name": t["name"], "description": t["description"], "input_schema": t["parameters"]} for t in TOOLS]
        self.usage = {"input": 0, "output": 0}

    def start(self, system, user):
        self.system, self.msgs = system, [{"role": "user", "content": user}]

    def step(self):
        r = self.c.messages.create(model=self.m, max_tokens=MAX_TOKENS, system=self.system, tools=self.tools, messages=self.msgs)
        self.usage["input"] += r.usage.input_tokens; self.usage["output"] += r.usage.output_tokens
        self.msgs.append({"role": "assistant", "content": [b.model_dump() for b in r.content]})
        text = "".join(b.text for b in r.content if b.type == "text")
        calls = [(b.id, b.name, b.input) for b in r.content if b.type == "tool_use"]
        return text, calls

    def add_results(self, results):
        self.msgs.append({"role": "user", "content": [{"type": "tool_result", "tool_use_id": i, "content": out}
                                                      for i, out in results]})


class OpenAI(Model):
    """Uses the Responses API, which OpenAI requires for reasoning models that call tools."""
    def nudge(self, text):
        self.input.append({"role": "user", "content": text})

    def __init__(self, model_id, key):
        import openai
        self.c, self.m = openai.OpenAI(api_key=key), model_id
        self.tools = [{"type": "function", **t} for t in TOOLS]
        self.usage = {"input": 0, "output": 0}

    def start(self, system, user):
        self.system, self.input = system, [{"role": "user", "content": user}]

    def step(self):
        r = self.c.responses.create(model=self.m, instructions=self.system, input=self.input, tools=self.tools,
                                    max_output_tokens=MAX_TOKENS, store=False, include=["reasoning.encrypted_content"])
        self.usage["input"] += r.usage.input_tokens; self.usage["output"] += r.usage.output_tokens
        for item in r.output:                      # reasoning, messages and calls all go back in next turn
            self.input.append(item.model_dump(exclude_none=True))
        calls = [(it.call_id, it.name, json.loads(it.arguments or "{}")) for it in r.output if it.type == "function_call"]
        return r.output_text or "", calls

    def add_results(self, results):
        for i, out in results:
            self.input.append({"type": "function_call_output", "call_id": i, "output": out})


class Google(Model):
    def nudge(self, text):
        self.contents.append(self.t.Content(role="user", parts=[self.t.Part(text=text)]))

    def __init__(self, model_id, key):
        from google import genai
        from google.genai import types
        self.t, self.c, self.m = types, genai.Client(api_key=key), model_id
        decls = [types.FunctionDeclaration(name=t["name"], description=t["description"], parameters=t["parameters"])
                 for t in TOOLS]
        self.cfg = types.GenerateContentConfig(
            tools=[types.Tool(function_declarations=decls)], max_output_tokens=MAX_TOKENS,
            automatic_function_calling=types.AutomaticFunctionCallingConfig(disable=True))
        self.usage = {"input": 0, "output": 0}

    def start(self, system, user):
        self.cfg.system_instruction = system
        self.contents = [self.t.Content(role="user", parts=[self.t.Part(text=user)])]

    def step(self):
        r = self.c.models.generate_content(model=self.m, contents=self.contents, config=self.cfg)
        um = r.usage_metadata
        self.usage["input"] += um.prompt_token_count or 0
        self.usage["output"] += (um.candidates_token_count or 0) + (getattr(um, "thoughts_token_count", 0) or 0)
        content = r.candidates[0].content
        self.contents.append(content)
        parts = content.parts or []
        text = "".join(p.text for p in parts if getattr(p, "text", None))
        calls = [(p.function_call.name + f"#{n}", p.function_call.name, dict(p.function_call.args or {}))
                 for n, p in enumerate(parts) if getattr(p, "function_call", None)]
        return text, calls

    def add_results(self, results):
        self.contents.append(self.t.Content(role="user", parts=[
            self.t.Part.from_function_response(name=i.split("#")[0], response={"output": out}) for i, out in results]))


class Mock(Model):
    def nudge(self, text):
        pass

    """Offline stand-in for testing the harness: looks at the files, then submits a fixed answer."""
    def __init__(self, *_):
        self.usage, self.n = {"input": 0, "output": 0}, 0

    def start(self, system, user):
        self.after = "messy_lab.duckdb" in user

    def step(self):
        self.n += 1
        if self.n == 1:
            code = ("import duckdb; c = duckdb.connect('messy_lab.duckdb', read_only=True); "
                    "print(c.execute(\"select count(*) from sample where study_id='S2' and role in "
                    "('untreated_control')\").fetchone())") if self.after else \
                   "import os; print(sorted(os.listdir('studies/data/raw')))"
            return "Looking at the data.", [("call1", "run_python", {"code": code})]
        ans = {"count": 31} if self.after else {"count": 15}
        return "", [("call2", "submit_answer", {"answer_json": json.dumps(ans), "explanation": "mock"})]

    def add_results(self, results):
        pass


ADAPTERS = {"anthropic": Anthropic, "openai": OpenAI, "google": Google, "mock": Mock}


def load_model(name):
    if name == "mock":
        return {"provider": "mock", "model": "mock"}, lambda: Mock()
    spec = yaml.safe_load(open(ROOT / "runs" / "models.yaml"))[name]
    if not spec.get("model"):
        sys.exit(f"Set the model ID for '{name}' in runs/models.yaml first.")
    key = os.environ.get(spec["key_env"])
    if not key:
        sys.exit(f"Set {spec['key_env']} in your environment first.")
    return spec, lambda: ADAPTERS[spec["provider"]](spec["model"], key)


# ---- rate limits ------------------------------------------------------------------------------
RETRYABLE = re.compile(r"429|RESOURCE_EXHAUSTED|rate.?limit|overloaded|529|503|UNAVAILABLE", re.I)


def step_with_retry(m, tries=8):
    """Call the model; on a rate-limit or overload error wait (as the provider asks) and retry."""
    waited = 0.0
    for attempt in range(tries):
        try:
            text, calls = m.step()
            return text, calls, waited
        except Exception as e:  # noqa: BLE001
            if attempt == tries - 1 or not RETRYABLE.search(str(e)):
                raise
            if re.search(r"per_day|PerDay", str(e)):
                raise RuntimeError("daily request quota used up for this model; it resets in "
                                   + (re.search(r"retry in ([\w.]+)", str(e)) or [None, "several hours"])[1]) from e
            hint = re.search(r"retryDelay': '(\d+)s", str(e))
            delay = float(hint.group(1)) + 2 if hint else min(60, 5 * 2 ** attempt)
            if delay > 300:
                raise RuntimeError(f"provider asks to wait {delay:.0f} s; stopping instead of hanging") from e
            print(f"   rate limited, waiting {delay:.0f} s", flush=True)
            time.sleep(delay)
            waited += delay


# ---- one question, one run ---------------------------------------------------------------------
def ask(model_name, spec, make_model, condition, run, qid, context, qtext):
    out_path = RESULTS / condition / model_name / f"run{run}" / f"{qid}.json"
    if out_path.exists():
        return json.loads(out_path.read_text())
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ws = make_workspace(condition)
    user = (f"{CONDITION_TEXT[condition]}\n\n{context}\n\nQuestion {qid}: {qtext}\n\n"
            f"Submit your answer with submit_answer, as JSON in exactly this shape: {FIELDS[qid]}")
    rec = {"model_name": model_name, "provider": spec["provider"], "model_id": spec["model"], "condition": condition,
           "run": run, "question": qid, "prompt": {"system": SYSTEM, "user": user},
           "limits": {"tool_calls": TOOL_CALLS, "seconds": SECONDS, "call_timeout": CALL_TIMEOUT, "output_cap": OUTPUT_CAP},
           "steps": [], "answer": None, "explanation": None, "stopped": None, "outside_access_attempt": False}
    t0, calls_used = time.time(), 0
    try:
        m = make_model()
        m.start(SYSTEM, user)
        while True:
            if time.time() - t0 > SECONDS:
                rec["stopped"] = "time limit"; break
            text, calls, waited = step_with_retry(m)
            t0 += waited                 # time spent waiting out a provider's rate limit doesn't count
            rec["rate_limit_wait_s"] = rec.get("rate_limit_wait_s", 0) + round(waited, 1)
            step = {"t": round(time.time() - t0, 1), "text": text, "calls": []}
            if not calls:
                rec["steps"].append(step)
                if not rec.get("nudged"):
                    rec["nudged"] = True     # same single reminder for every model
                    m.nudge("Please call submit_answer now with your answer as JSON.")
                    continue
                rec["stopped"] = "ended without submitting"; break
            results, submitted = [], False
            for cid, name, args in calls:
                if name == "submit_answer":
                    raw = args.get("answer_json", "")
                    try:
                        rec["answer"] = json.loads(raw)
                    except (json.JSONDecodeError, TypeError):
                        rec["answer"] = {"_unparseable": raw}
                    rec["explanation"] = args.get("explanation")
                    step["calls"].append({"tool": name, "args": args})
                    submitted = True
                    results.append((cid, "Answer received."))
                elif name == "run_python":
                    calls_used += 1
                    if calls_used > TOOL_CALLS:
                        out, flag = f"[tool-call limit of {TOOL_CALLS} reached; submit your answer now]", False
                    else:
                        out, flag = run_python(args.get("code", ""), ws)
                    rec["outside_access_attempt"] |= flag
                    step["calls"].append({"tool": name, "code": args.get("code", ""), "output": out, "flagged": flag})
                    results.append((cid, out))
                else:
                    results.append((cid, f"Unknown tool {name}."))
            rec["steps"].append(step)
            if submitted:
                rec["stopped"] = "submitted"; break
            if calls_used > TOOL_CALLS + 2:
                rec["stopped"] = "tool-call limit"; break
            m.add_results(results)
        rec["usage"] = m.usage
    except Exception as e:  # an API error is recorded, not hidden
        rec["stopped"] = f"error: {type(e).__name__}: {str(e)[:500]}"
        rec["usage"] = getattr(locals().get("m"), "usage", None)
        rec["traceback"] = traceback.format_exc()[-3000:]
    finally:
        shutil.rmtree(ws, ignore_errors=True)
    rec["seconds"], rec["tool_calls"] = round(time.time() - t0, 1), calls_used
    price, u = spec.get("price_per_mtok"), rec.get("usage") or {}
    rec["cost_usd"] = round(u.get("input", 0) / 1e6 * price["input"] + u.get("output", 0) / 1e6 * price["output"], 4) if price else 0.0
    ok, parts = grade(qid, rec["answer"] or {})
    rec["correct"] = bool(ok) and not rec["outside_access_attempt"]
    rec["grade_parts"] = [{"field": f, "ok": o, "got": g} for f, o, g in parts]
    if rec["stopped"] and rec["stopped"].startswith("error"):
        out_path.with_suffix(".error.json").write_text(json.dumps(rec, indent=2, default=str))  # retried next time
    else:
        out_path.write_text(json.dumps(rec, indent=2, default=str))
        out_path.with_suffix(".error.json").unlink(missing_ok=True)   # an earlier failed try is superseded
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="mock")
    ap.add_argument("--conditions", default="before,after")
    ap.add_argument("--runs", type=int, default=1)
    ap.add_argument("--questions", default="all")
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()

    check = subprocess.run([sys.executable, "-c", "import pandas, numpy, duckdb, pypdf"], capture_output=True, text=True,
                           env={"PATH": os.environ.get("PATH", ""), "PYTHONUSERBASE": site.getuserbase()})
    if check.returncode:
        sys.exit("The sandbox is missing a library the instructions promise. Run: pip install -r requirements.txt\n" + check.stderr)
    context, qs = load_questions()
    qids = sorted(qs) if a.questions == "all" else a.questions.split(",")
    assert set(qids) <= set(FIELDS), "unknown question id"
    jobs = []
    for name in a.models.split(","):
        spec, make = load_model(name)
        for cond in a.conditions.split(","):
            for run in range(1, a.runs + 1):
                for q in qids:
                    jobs.append((name, spec, make, cond, run, q, context, qs[q]))
    print(f"{len(jobs)} answers to collect ({a.models} x {a.conditions} x {a.runs} run(s) x {len(qids)} questions)")
    done = 0
    with cf.ThreadPoolExecutor(max_workers=a.workers) as ex:
        for rec in ex.map(lambda j: ask(*j), jobs):
            done += 1
            tag = "OK " if rec["correct"] else "   "
            print(f"[{done}/{len(jobs)}] {tag} {rec['condition']:6} {rec['model_name']:7} run{rec['run']} {rec['question']} "
                  f"{rec['stopped']} ({rec['seconds']} s, {rec['tool_calls']} calls, ${rec.get('cost_usd', 0):.3f})")


if __name__ == "__main__":
    main()
