#!/usr/bin/env python3
"""Turn \\cite{key} markers into numbered [n] citations in order of first
appearance, and emit the matching numbered reference list (PRE style).

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
    return ", ".join(out[:8]) + (" et al." if len(out) > 8 else "")

order = []

def repl(m):
    nums = []
    for k in (k.strip() for k in m.group(1).split(",")):
        if k not in entries:
            sys.exit(f"MISSING BIB KEY: {k}")
        if k not in order:
            order.append(k)
        nums.append(order.index(k) + 1)
    return "[" + ", ".join(str(n) for n in sorted(nums)) + "]"

body = re.sub(r"\\cite\{([^}]+)\}", repl,
              "\n\n".join(open(f, encoding="utf-8").read().rstrip() for f in SRCS))

lines = ["", "---", "", "## 参考文献", ""]
for i, k in enumerate(order, 1):
    e = entries[k]
    bits = [f"{authors(e['author'])}, *{detex(e.get('journal') or e.get('booktitle') or '')}*"]
    if e.get("volume"):
        bits.append(f"**{e['volume']}**")
    if e.get("pages"):
        bits.append(detex(e["pages"]).replace("--", "–"))
    lines.append(f"{i}. " + ", ".join(bits) + f" ({e['year']}). doi:{e.get('doi', '')}")

open(OUT, "w", encoding="utf-8").write(body.rstrip() + "\n" + "\n".join(lines) + "\n")
print(f"{len(SRCS)} source file(s) -> {OUT}: {len(order)} references cited, numbered in order of appearance")
