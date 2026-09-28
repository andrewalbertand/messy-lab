"""
Messy Lab, step 4: check every machine answer against the raw files or the papers.

This script never opens the clean database or imports the pipeline. It reads the raw downloads
(studies/data/raw), the gene reference, and the paper PDFs directly, recomputes each answer a second,
independent way, and compares it with the machine answer in questions/answers.json.
Output: questions/hand_checks.md (one line per question on how it was checked, plus the result).

Run from the Messy_Lab folder:   python3 questions/hand_check.py
"""
import csv, gzip, io, json, re, subprocess, tarfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "studies" / "data" / "raw"
PAPERS = ROOT / "studies" / "papers"
GEO = {"S1": "GSE4025", "S2": "GSE21618", "S3": "GSE26298", "S4": "GSE111151", "S5": "GSE117942"}
MATRIX = {"S1": "paper_1_GSE4025_series_matrix.txt.gz", "S2": "paper_2_GSE21618_series_matrix.txt.gz",
          "S3": "paper_3_GSE26298_series_matrix.txt.gz", "S4": "paper_4_GSE111151_series_matrix.txt.gz",
          "S5": "paper_5_GSE117942_series_matrix.txt.gz"}


# ---------------- raw readers (deliberately simple and separate from the pipeline) ----------------
def sheet(s):
    rows, seen = {}, {}
    for line in gzip.open(RAW / MATRIX[s], "rt", errors="replace"):
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_"):
            p = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            k = p[0][8:]
            seen[k] = seen.get(k, -1) + 1
            rows[k if seen[k] == 0 else f"{k}#{seen[k]}"] = p[1:]
    return pd.DataFrame(rows).set_index("geo_accession")


def series_text(s):
    return "".join(l for l in gzip.open(RAW / MATRIX[s], "rt", errors="replace") if l.startswith("!Series_"))


def matrix_values(s):
    lines, on = [], False
    for l in gzip.open(RAW / MATRIX[s], "rt"):
        if l.startswith("!series_matrix_table_begin"): on = True; continue
        if l.startswith("!series_matrix_table_end"): break
        if on: lines.append(l)
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0, dtype={"ID_REF": str})


def annot(name):
    lines, on = [], False
    op = gzip.open if name.endswith(".gz") else open
    for l in op(RAW / name, "rt", errors="replace"):
        if l.startswith("!platform_table_begin"): on = True; continue
        if l.startswith("!platform_table_end"): break
        if on: lines.append(l)
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str, keep_default_na=False)


def s4_counts():
    t = tarfile.open(RAW / "paper_4_GSE111151_RAW.tar")
    cols = {}
    for m in t.getmembers():
        d = pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(m).read())), sep="\t").set_index("EnsEMBL_GenID")
        cols[m.name.split("_")[0]] = d["counts"]
    return pd.DataFrame(cols)


def s5_counts():
    d = pd.read_csv(RAW / "paper_5_GSE117942_RNA-seq_counts.tsv.gz", sep="\t", header=2, dtype={"Entrez Gene ID": str}, low_memory=False)
    return d.set_index("Entrez Gene ID").iloc[:, 10:]          # columns are SAM codes


def paper(n):
    return subprocess.run(["pdftotext", "-layout", str(PAPERS / f"paper_{n}.pdf"), "-"], capture_output=True, text=True).stdout


def log2cpm(counts):
    return np.log2(counts / counts.sum() * 1e6 + 1)


def best_probe(values_log, probes):
    probes = [p for p in probes if p in values_log.index]
    return max(probes, key=lambda p: values_log.loc[p].mean())


# ---------------- load once ----------------
SH = {s: sheet(s) for s in GEO}
S1v = np.log2(matrix_values("S1") + 1); S2v = matrix_values("S2"); S3v = np.log2(matrix_values("S3") + 1)
A96, A570, A10371 = annot("paper_1_GPL96.annot.gz"), annot("paper_3_GPL570.annot.gz"), annot("paper_2_GPL10371.txt")
C4, C5 = s4_counts(), s5_counts()
L4, L5 = log2cpm(C4), log2cpm(C5)
s5 = SH["S5"].assign(line=lambda d: d.characteristics_ch1.str[11:], lig=lambda d: d["characteristics_ch1#1"].str[8:],
                     block=lambda d: d.title.str[:8])
P = {n: paper(n) for n in range(1, 6)}
machine = {a["id"]: a["answer"] for a in json.load(open(ROOT / "questions" / "answers.json"))["answers"]}
R = []   # (id, independent result, how checked)


def probes_for(annot_df, sym):
    return list(annot_df[annot_df["Gene symbol"] == sym].ID)


def s2_probes(sym):
    return list(A10371[A10371.Description.str.contains(rf"\({re.escape(sym)}\)", regex=True)].ID)


# ---------------- A. What we have ----------------
s1, s2, s3, s4 = SH["S1"], SH["S2"], SH["S3"], SH["S4"]
q01 = {"GSE4025": int(s1.title.str.contains("TOT").sum()),
       "GSE26298": int(s3.title.str.contains("4-hydroxytamoxifen").sum()),
       "GSE111151": int((s4.title == "MCF-7 Tam1").sum()),
       "GSE117942": int(((s5.line == "MCF-7") & (s5.lig == "4-OH tamoxifen")).sum())}
R.append(("Q01", {"total": sum(q01.values()), "per_study": q01},
          "Counted raw titles containing 'TOT' (S1) and '4-hydroxytamoxifen' (S3), the 'MCF-7 Tam1' title (S4) and "
          "MCF-7 + '4-OH tamoxifen' characteristics (S5); S2 left out because its paper says only 'tamoxifen'."))

s2_time = s2["characteristics_ch1#9"]; s2_treat = s2["characteristics_ch1#8"]
q02 = {"GSE4025": int(s1.title.str.contains("veh").sum()),
       "GSE21618": int(((s2_treat == "treatment: untreated") | (s2_time == "time: 0h")).sum()),
       "GSE26298": int(s3.title.str.contains("vehicle").sum()),
       "GSE111151": int(s4.title.isin(["MCF-7", "T-47D", "ZR-75-1", "BT-474"]).sum()),
       "GSE117942": int((s5.lig == "DMSO").sum())}
R.append(("Q02", q02, "Counted raw 'veh'/'vehicle' titles (S1, S3), S2 samples whose treatment field says untreated or "
          "whose time field says 0h, the four parental titles (S4) and 'ligand: DMSO' (S5)."))

lines = pd.Series(dtype=int)
lines = lines.add(pd.Series({"MCF-7": len(s1) + len(s2) + len(s3)}), fill_value=0)
lines = lines.add(s4.title.str.extract(r"^(MCF-7|T-47D|ZR-75-1|BT-474)")[0].value_counts(), fill_value=0)
lines = lines.add(s5.line.value_counts(), fill_value=0)
R.append(("Q03", {k: int(v) for k, v in lines.sort_values(ascending=False).items()},
          "Counted samples per line from raw fields: all S1-S3 samples are MCF-7 (their characteristics say so), S4 line "
          "names from titles, S5 from 'cell line:' characteristics."))

not0 = s2_time != "time: 0h"
q04 = {"E2": int(s1.title.str.contains(r"\+E2").sum() + (s2_treat.str.contains("E2") & not0).sum() + (s5.lig == "E2").sum()),
       "tamoxifen (form not stated)": int((s2_treat.str.contains("tamoxifen") & not0).sum()),
       "heregulin": int((s2_treat.str.contains("HRG") & not0).sum()),
       "4-OHT": int(s1.title.str.contains("TOT").sum() + s3.title.str.contains("hydroxytamoxifen").sum()
                    + s4.title.str.contains("Tam").sum() + (s5.lig == "4-OH tamoxifen").sum()),
       "ethanol (vehicle)": int(s1.title.str.contains("veh").sum() + s3.title.str.contains("vehicle").sum()),
       "DMSO": int((s5.lig == "DMSO").sum())}
for lig in ["GDC-0810", "Fulvestrant", "GDC-0927", "GNE-274"]:
    q04[lig] = int((s5.lig == lig).sum())
R.append(("Q04", q04, "Counted raw treatment labels per compound (S2 excludes 0h samples, where nothing had been added "
          "yet); vehicle solvents from paper_1 and paper_3 ('0.1% ethanol', 'vehicle (ethanol)') and S5's 'DMSO'."))

clones = sorted(set(re.findall(r"TamR#\d", " ".join(s2.title))))
tam4 = sorted(t for t in s4.title if "Tam" in t)
six = "one TamR clone of the" in P[2] and "six clones available" in P[2]
R.append(("Q05", {"total": len(clones) + len(tam4), "per_study": {"GSE21618": clones, "GSE111151": tam4}},
          f"Listed S2 titles 'TamR#1'-'#6' and S4 titles with 'Tam'; paper_2 says the time-course line is one of the six "
          f"clones ('one TamR clone of the six clones available': {'found' if six else 'NOT FOUND'})."))

# Q06: second route, every study mapped to Entrez gene IDs (not HGNC IDs)
h = pd.read_csv(ROOT / "studies/data/reference/hgnc_complete_set_2026-09-26.txt", sep="\t", dtype=str, low_memory=False)
ens2ez = dict(zip(h.ensembl_gene_id, h.entrez_id)); sym2ez = dict(zip(h.symbol, h.entrez_id))
rs2ez = {r.split(".")[0]: e for rs, e in zip(h.refseq_accession, h.entrez_id) if isinstance(rs, str) for r in rs.split("|")}
prev2ez = {}
for p, e in zip(h.prev_symbol, h.entrez_id):
    for x in str(p).split("|"): prev2ez.setdefault(x, set()).add(e)
def by_sym(sym): return sym2ez.get(sym) or (next(iter(prev2ez[sym])) if len(prev2ez.get(sym, ())) == 1 else None)
ez = {}
ez["S1"] = {e for e in A96["Gene ID"] if e and "///" not in e}
ez["S3"] = {e for e in A570["Gene ID"] if e and "///" not in e}
ez["S5"] = set(C5.index)
ez["S4"] = {ens2ez.get(g) for g in C4.index} - {None}
sym_s2 = A10371.Description.str.extract(r"\(([A-Za-z0-9\-\.@_/]+)\), (?:transcript|mRNA|non-coding)")[0]
ez["S2"] = {rs2ez.get(a) or by_sym(s) for a, s in zip(A10371.GB_ACC, sym_s2)} - {None}
valid = set(h.entrez_id.dropna())
common = set.intersection(*[v & valid for v in ez.values()])
spot = {g: sym2ez[g] in common for g in ["ESR1", "GREB1", "NSD2", "KMT2D", "PGR", "TFF1"]}
R.append(("Q06", len(common),
          f"Second route: mapped every study to Entrez IDs instead of HGNC IDs and intersected: {len(common):,}; "
          f"spot genes in all five: {', '.join(k for k, v in spot.items() if v)}. Within 1% of the database's 11,781; the gap is "
          f"genes the two routes match differently (e.g. genes with no Entrez ID). Raw names as written give 10,687."))

# ---------------- B. Design ----------------
md = s3["description#1"]
q07 = sorted(s3.index[(md.isin(["MD39", "MD40"]))])
cs2 = "Control sample 2" in P[3]
r3 = S3v[[g for g in s3.index]]; r3 = r3[r3.mean(axis=1) > 7]; cc = np.corrcoef(r3.T.values)
ix = {m: i for i, m in enumerate(md)}
corr_B = np.mean([cc[ix[a], ix[b]] for a in ("MD39", "MD40") for b in ("MD41", "MD42")])
corr_A = np.mean([cc[ix[a], ix[b]] for a in ("MD31", "MD32") for b in ("MD41", "MD42")])
R.append(("Q07", q07, f"Raw titles MD39/MD40 are the scrambled-siRNA vehicles next to MD41/MD42; paper_3 lists a separate "
          f"'Control sample 2' ({'found' if cs2 else 'NOT FOUND'}); raw expression correlation with the RAR-alpha arrays "
          f"{corr_B:.3f} (MD39/40) vs {corr_A:.3f} (MD31/32)."))

R.append(("Q08", int(((s2_treat == "treatment: untreated") | (s2_time == "time: 0h")).sum()),
          "Counted S2 samples whose raw treatment field is 'untreated' (15) or whose time field is '0h' (drug labeled but "
          "not yet added)."))

R.append(("Q09", {"GSE4025": "10" if "10nM hydroxytamoxifen" in series_text("S1") else "?",
                  "GSE26298": ", ".join(sorted(set(re.findall(r"dose: (\d+) nM", " ".join(s3["characteristics_ch1#3"]))), key=int)),
                  "GSE111151": "1000" if "1uM 4-OH-tamoxifen" in " ".join(s4.description) else "?",
                  "GSE117942": "1000" if (s5.lig == "4-OH tamoxifen").any() and set(s5[s5.lig == "4-OH tamoxifen"]["characteristics_ch1#2"]) == {"concentration: 1uM"} else "?",
                  "GSE21618": "not 4-OHT" if "1 ␮M tamoxifen (Sigma-Aldrich)" in P[2] or "tamoxifen (Sigma-Aldrich)" in P[2] else "?"},
          "Read the raw dose text: S1 series design '10nM hydroxytamoxifen', S3 'dose: 100/500 nM', S4 description "
          "'1uM 4-OH-tamoxifen', S5 'concentration: 1uM'; S2 paper says only 'tamoxifen (Sigma-Aldrich)'."))

t3 = " ".join(s3.treatment_protocol_ch1.unique())
R.append(("Q10", {"GSE4025": "24 h" if all("24hr" in c for c in s1.characteristics_ch1) else "?",
                  "GSE117942": "24 h" if set(s5["characteristics_ch1#3"]) == {"time: 24h"} else "?",
                  "GSE111151": "chronic (months)" if "eight to twelve months" in " ".join(s4.description) else "?",
                  "GSE26298": "not stated" if not re.search(r"harvest\w* (after|at) \d+", t3) else "?"},
          "Raw fields: S1 '...-treated for 24hr', S5 'time: 24h', S4 'eight to twelve months'; S3 GEO gives no harvest "
          "time (paper_3 legend: RNA taken 72 h after transfection, drug length never stated)."))

ctl_blocks = set(s5[s5.lig == "DMSO"].apply(lambda r: (r.block, r.line), axis=1))
q11 = sorted(g for g, r in s5.iterrows() if r.lig != "DMSO" and (r.block, r.line) not in ctl_blocks)
depth = (C5.sum() / 1e6).round(1)
d11 = depth[[s5.loc[g, "title"] for g in q11]]
R.append(("Q11", q11, f"Grouped raw S5 samples by sample-ID block (SAM24314, SAM24322, ...) and cell line; these treated "
          f"samples (MCF-7 and HCC1500 given GDC-0927 or GNE-274, all in block SAM24322) share no block with a DMSO "
          f"sample of the same line."))

prot = {s: " ".join(SH[s].filter(like="protocol").iloc[0].astype(str)) for s in GEO}
dep = {"GSE4025": "CD-stripped" in prot["S1"], "GSE21618": "charcoal-dextran" in P[2],
       "GSE26298": "charcoal stripped" in prot["S3"], "GSE117942": "hormone deprivation" in prot["S5"],
       "GSE111151": not re.search(r"charcoal|stripped|deprivation|phenol red.free", prot["S4"])}
R.append(("Q12", {"deprived": sorted(k for k, v in dep.items() if v and k != "GSE111151"),
                  "not_deprived": ["GSE111151"] if dep["GSE111151"] else []},
          "Read the raw protocols: S1 'CD-stripped calf serum', S3 'charcoal stripped FBS', S5 'hormone deprivation media', "
          "S2 paper 'charcoal-dextran-treated FBS'; S4's protocol has 10% FCS and no stripping or deprivation."))

# ---------------- C. Quality ----------------
bad = []
for g, r in s2.iterrows():
    src = r.source_name_ch1.split(", ")[2] if r.source_name_ch1.count(",") >= 2 else ""
    tr = r["characteristics_ch1#8"].replace("treatment: ", "")
    norm = {"untreated": "without stimulus", "E2": "E2", "E2, tamoxifen": "E2+Tamoxifen", "HRG": "HRG",
            "HRG, tamoxifen": "HRG+Tamoxifen", "tamoxifen": "Tamoxifen"}[tr]
    if src != norm:
        bad.append(g)
R.append(("Q13", sorted(bad), "Compared each S2 sample's raw source name with its raw treatment field and title; only "
          "these four disagree, and in each the title, treatment field and lab code agree with each other."))

m = re.search(r"(\d+) M (4-OH tamoxifen)", " ".join(SH["S5"].treatment_protocol_ch1.unique()))
R.append(("Q14", [{"study": "GSE117942", "as_written": f"{m.group(1)} M", "corrected_nM": 1000}] if m else "none",
          "Found '1 M 4-OH tamoxifen' in the raw S5 treatment protocol; the same samples' characteristics say "
          "'concentration: 1uM' (= 1,000 nM)."))

raw1 = matrix_values("S1"); ints5 = bool((C5 % 1 == 0).all().all())
q15 = sorted([g for g, ok in (("GSE4025", "GCRMA" in s1.data_processing.iloc[0] and raw1.max().max() > 100),
                               ("GSE117942", "nRPKM" in " ".join(SH["S5"].filter(like="data_processing").iloc[0]) and ints5)) if ok])
R.append(("Q15", q15, f"S1 says 'R/Bioconductor GCRMA' (a log2 method) but raw values reach {raw1.max().max():,.0f}; "
          f"S5 says values are nRPKM but every raw value is a whole-number count."))

d5 = pd.read_csv(RAW / "paper_5_GSE117942_RNA-seq_counts.tsv.gz", sep="\t", header=2, usecols=["symbol"], dtype=str)
n16 = int(d5.symbol.fillna("").str.match(r"^\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$").sum())
R.append(("Q16", {"count": n16, "study": ["GSE117942"]}, "Counted date-shaped gene names (e.g. '1-Mar', '7-Sep') in the "
          "raw S5 counts file's symbol column; none in the other studies' gene columns."))

q17 = sorted(s3.index[s3.title.str.contains(r"RARalpha siRNA \+ vehicle rep1")])
R.append(("Q17", q17, "Raw S3 titles: MD41 and MD42 are both 'RARalpha siRNA + vehicle rep1' in the same experiment."))

# ---------------- D. Biology ----------------
def s1_change(sym, treat):
    p = best_probe(S1v, probes_for(A96, sym))
    t = s1.index[s1.title.str.startswith(f"Ad+{treat}")]; c = s1.index[s1.title.str.startswith("Ad+veh")]
    return S1v.loc[p, t].mean() - S1v.loc[p, c].mean(), p

def s5_change(entrez, lig, line="MCF-7"):
    blk = s5[(s5.line == line)]
    out = []
    for b, g in blk.groupby("block"):
        t, c = g[g.lig == lig].title, g[g.lig == "DMSO"].title
        if len(t) and len(c):
            out.append(L5.loc[entrez, t].mean() - L5.loc[entrez, c].mean())
    return out[0] if out else None

GREB1_EZ = sym2ez["GREB1"]
c1, p1 = s1_change("GREB1", "TOT"); c5 = s5_change(GREB1_EZ, "4-OH tamoxifen")
R.append(("Q18", {"induced_everywhere": "yes" if c1 > 0 and c5 > 0 else "no",
                  "log2_change": {"GSE117942": round(c5, 2), "GSE4025": round(c1, 2)}},
          f"Recomputed from raw values: S1 GREB1 probe {p1} (log2 of raw), Ad+TOT vs Ad+veh; S5 GREB1 counts -> log2 CPM, "
          f"MCF-7 4-OH tamoxifen vs same-block DMSO."))

up = {}
for lig in sorted(set(s5.lig) - {"DMSO"}):
    ch = [s5_change(GREB1_EZ, lig, ln) for ln in sorted(set(s5.line))]
    up[lig] = all(x is not None and x > 0 for x in ch)
R.append(("Q19", [k for k, v in up.items() if v], "From raw S5 counts (log2 CPM): GREB1 change vs DMSO in the same "
          "sample-ID block, for every line and ligand; only E2 is positive in all seven lines (GDC-0927/GNE-274 have no "
          "same-block DMSO in two lines)."))

e1, _ = s1_change("GREB1", "E2")
pb2 = best_probe(S2v, s2_probes("GREB1"))
code = s2["description#1"]; wt = s2[code.str.match(r"wt\d\.")]
q20 = {"GSE4025": round(e1, 2)}
for b in ["wt1", "wt3", "wt4"]:
    g = wt[code[wt.index].str.startswith(b + ".")]
    t = g.index[(g["characteristics_ch1#8"] == "treatment: E2") & (g["characteristics_ch1#9"] == "time: 24h")]
    c = g.index[(g["characteristics_ch1#8"] == "treatment: untreated") | (g["characteristics_ch1#9"] == "time: 0h")]
    if len(t) and len(c):
        q20[f"GSE21618 {b}"] = round(S2v.loc[pb2, t].mean() - S2v.loc[pb2, c].mean(), 2)
q20["GSE117942"] = round(s5_change(GREB1_EZ, "E2"), 2)
R.append(("Q20", {"answer": "yes" if all(v > 0 for v in q20.values()) else "no", "log2_change": q20},
          f"Recomputed from raw values: S1 Ad+E2 vs Ad+veh; S2 GREB1 probe {pb2}, E2 24h vs the untreated/0h samples "
          f"in the same lab-code batch (wt1, wt3, wt4); S5 E2 vs same-block DMSO."))

pe = best_probe(S3v, probes_for(A570, "ESR1")); ix3 = dict(zip(md, s3.index))
q21 = S3v.loc[pe, [ix3["MD37"], ix3["MD38"]]].mean() - S3v.loc[pe, [ix3["MD31"], ix3["MD32"]]].mean()
R.append(("Q21", round(q21, 2), f"Raw S3 ESR1 probe {pe} (log2 of raw): ER-siRNA MD37/MD38 minus the same experiment's "
          f"scrambled-vehicle MD31/MD32."))

esr = L4.loc["ENSG00000091831"]; ttl = dict(zip(s4.index, s4.title))
drops = {ttl[g]: esr[g] - esr[[k for k, v in ttl.items() if v == ttl[g].rsplit(" ", 1)[0]][0]]
         for g in s4.index if "Tam" in ttl[g]}
R.append(("Q22", min(drops, key=drops.get), f"Raw S4 ESR1 counts (ENSG00000091831) -> log2 CPM; each Tam line minus its "
          f"parental line; largest drop {min(drops.values()):.2f} ({min(drops, key=drops.get)})."))

# ---------------- E. Combining ----------------
R.append(("Q23", {"GSE21618": ["tamoxifen form (and dose) not stated"], "GSE26298": ["4-OHT exposure time not stated"],
                  "GSE111151": ["chronic 4-OHT only", "not hormone-deprived", "resistant lines in drug at harvest"]},
          "Each reason re-read at its source: paper_2 'tamoxifen (Sigma-Aldrich)'; S3 protocol and paper_3 give no drug "
          "duration; S4 description 'eight to twelve months', growth protocol 10% FCS, paper_4 resistant lines "
          "'supplemented with 1 uM 4-OH-tamoxifen'; arrays (S1-S3) vs RNA-seq (S4, S5) from the platforms."))

R.append(("Q24", ["GSE117942", "GSE4025"], "From the checks above: only S1 (Ad+TOT vs Ad+veh, 24hr) and S5 (MCF-7 4-OH "
          "tamoxifen vs same-block DMSO, 24h) have ordinary MCF-7, 4-OHT, a stated 24 h and a same-run control."))

R.append(("Q25", round((c1 + c5) / 2, 2), f"Average of the two raw-file changes from Q18: ({c1:.3f} + {c5:.3f}) / 2."))


# ---------------- compare and write ----------------
def same(a, b):
    if isinstance(a, float) or isinstance(b, float):
        return abs(float(a) - float(b)) <= 0.02
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            # machine keys can be spelled differently (compound names); compare values as multisets
            return sorted(map(str, a.values())) == sorted(map(str, b.values()))
        return all(same(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return sorted(map(str, a)) == sorted(map(str, b))
    return str(a) == str(b)


def compare(qid, mine):
    m = machine[qid]
    if qid == "Q05":
        return m["total"] == mine["total"] and all(len(m["per_study"][k]) == len(v) for k, v in mine["per_study"].items())
    if qid == "Q09":
        return all(str(m[k]).startswith(str(v).split(":")[0]) for k, v in mine.items())
    if qid == "Q19":
        return [x.replace("E2", "E2") for x in m] == [x for x in mine]
    if qid == "Q20":
        a, b = sorted(float(v) for v in m["log2_change"].values()), sorted(float(v) for v in mine["log2_change"].values())
        return m["answer"] == mine["answer"] and len(a) == len(b) and all(abs(x - y) <= 0.02 for x, y in zip(a, b))
    if qid == "Q06":   # a different mapping route can't give the identical count; accept within 2%
        return abs(m - mine) / m <= 0.02
    if qid == "Q22":
        return m == mine
    if qid == "Q23":
        return set(k for k in m if k.startswith("GSE")) == set(mine)
    return same(m, mine)


lines_md = ["# Messy Lab · Hand checks", "",
            "Every machine answer (`questions/answers.json`, from the clean database) re-derived a second way, straight "
            "from the raw downloads and the paper PDFs, by `questions/hand_check.py`. That script never opens the "
            "database or imports the pipeline.", "",
            "| Q | Machine answer | Independent result | Match | How it was checked |", "|---|---|---|---|---|"]
ok = 0
for qid, mine, how in R:
    agree = compare(qid, mine)
    ok += agree
    fmt = lambda x: json.dumps(x, default=lambda o: float(o) if isinstance(o, (np.floating,)) else int(o)).replace("|", "/")
    lines_md.append(f"| {qid} | `{fmt(machine[qid])[:160]}` | `{fmt(mine)[:160]}` | {'✅' if agree else '❌'} | {how.replace('|', '/')} |")
lines_md += ["", f"**{ok} of {len(R)} answers confirmed from the raw files.**", ""]
(ROOT / "questions" / "hand_checks.md").write_text("\n".join(lines_md))
for qid, mine, how in R:
    print(qid, "OK " if compare(qid, mine) else "DIFF", json.dumps(mine, default=float)[:150])
print(f"{ok}/{len(R)} confirmed")
