#!/usr/bin/env python3
"""Turn \\cite{key} markers into numbered [n] citations against an alphabetical
reference list, and emit that list (SIAM style; --order appearance for PRE).

    python3 tools/number_citations.py introduction.src.md introduction.md
    python3 tools/number_citations.py intro.src.md section2.src.md -o manuscript.md

Multiple sources are concatenated in order and share one numbering sequence,
so sections of the same manuscript never restart at [1].

Keeps manuscript numbering mechanical: hand-numbering breaks the moment a
citation is inserted mid-draft.
"""
import re, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
BIB = os.path.join(HERE, "..", "references.bib")
args = sys.argv[1:]
ORDER = "alpha"          # SIAM: alphabetical list, numeric labels
if "--order" in args:
    i = args.index("--order")
    ORDER = args[i + 1]
    del args[i:i + 2]
if "-o" in args:
    i = args.index("-o")
    SRCS, OUT = args[:i], args[i + 1]
else:
    SRCS, OUT = args[:1], args[1]

raw = open(BIB, encoding="utf-8").read()
entries = {}
for m in re.finditer(r"@(?:article|inproceedings|incollection|book)\{([^,]+),([\s\S]*?)\n\}", raw):
    f = {}
    for fm in re.finditer(r"(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}", m.group(2)):
        f[fm.group(1).lower()] = fm.group(2)
    entries[m.group(1).strip()] = f

ACC = {"{\\~n}": "ñ", "{\\'a}": "á", "{\\'e}": "é", "{\\'i}": "í", "{\\'o}": "ó",
       "{\\'u}": "ú", '{\\"u}': "ü", '{\\"o}': "ö", "{\\'A}": "Á", "{\\~a}": "ã", "{\\`e}": "è"}

def detex(s):
    if not s:
        return ""
    for k, v in ACC.items():
        s = s.replace(k, v)
    s = re.sub(r"\{\\'([a-zA-Z])\}", r"\1", s)
    s = re.sub(r'\{\\"([a-zA-Z])\}', r"\1", s)
    s = re.sub(r"\{\\~([a-zA-Z])\}", r"\1", s)
    return re.sub(r"\s+", " ", re.sub(r"[{}\\]", "", s)).strip()

def authors(a):
    out = []
    for n in detex(a).split(" and "):
        p = [x.strip() for x in n.split(",")]
        if len(p) == 2:
            fam, giv = p
            out.append(" ".join(w[0] + "." for w in giv.replace("-", " ").split() if w) + " " + fam)
        else:
            out.append(n.strip())
    out = out[:8] + (["et al."] if len(out) > 8 else [])
    # SIAM joins the last two authors with "and".
    return out[0] if len(out) == 1 else ", ".join(out[:-1]) + (
        ", and " if len(out) > 2 else " and ") + out[-1]

def family(a):
    """Family name of the first author, for alphabetical ordering."""
    first = detex(a).split(" and ")[0]
    return (first.split(",")[0] if "," in first else first.split()[-1]).lower()

# SIAM orders the reference list alphabetically and labels it numerically, so the
# numbers are only known after every citation has been seen: collect first, then
# number, then substitute.
CITE = re.compile(r"\\cite\{([^}]+)\}")
body = "\n\n".join(open(f, encoding="utf-8").read().rstrip() for f in SRCS)
cited = []
for m in CITE.finditer(body):
    for k in (k.strip() for k in m.group(1).split(",")):
        if k not in entries:
            sys.exit(f"MISSING BIB KEY: {k}")
        if k not in cited:
            cited.append(k)

if ORDER == "appearance":
    order = cited
else:
    order = sorted(cited, key=lambda k: (family(entries[k]["author"]),
                                         entries[k].get("year", ""),
                                         detex(entries[k].get("title", ""))))
num = {k: i + 1 for i, k in enumerate(order)}

body = CITE.sub(lambda m: "[" + ", ".join(
    str(n) for n in sorted(num[k.strip()] for k in m.group(1).split(","))) + "]", body)

lines = ["", "---", "", "## 参考文献", ""]
for i, k in enumerate(order, 1):
    e = entries[k]
    bits = [authors(e["author"])]
    if e.get("title"):
        bits.append(f"*{detex(e['title'])}*")
    tail = detex(e.get("journal") or e.get("booktitle") or "")
    if e.get("volume"):
        tail += f", {e['volume']}"
    tail += f" ({e['year']})"
    if e.get("pages"):
        pg = detex(e["pages"]).replace("--", "–")
        tail += (", pp. " if "–" in pg else ", p. ") + pg
    bits.append(tail)
    doi = e.get("doi", "")
    lines.append(f"{i}. " + ", ".join(bits) + (f", https://doi.org/{doi}." if doi else "."))

open(OUT, "w", encoding="utf-8").write(body.rstrip() + "\n" + "\n".join(lines) + "\n")
print(f"{len(SRCS)} source file(s) -> {OUT}: {len(order)} references, ordered by {ORDER}")
