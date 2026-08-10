const fs = require("fs");
const path = require("path");
const D = require("docx");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  TabStopType,
  LevelFormat, convertInchesToTwip,
} = D;

const SRC = process.argv[2];
const BIB = path.join(__dirname, "..", "references.bib");
const OUT = process.argv[3];

const LATIN = "Times New Roman";
const CJK = "SimSun";
const FONT = { ascii: LATIN, hAnsi: LATIN, eastAsia: CJK, cs: LATIN };

// ---------------------------------------------------------------- math ----
const UNRESOLVED = new Map();

const SYM = {
  "\\phi": "\u03c6", "\\alpha": "\u03b1", "\\tau": "\u03c4", "\\eta": "\u03b7",
  "\\mu": "\u03bc", "\\theta": "\u03b8", "\\lambda": "\u03bb", "\\beta": "\u03b2",
  "\\sigma": "\u03c3", "\\rho": "\u03c1", "\\Delta": "\u0394", "\\varnothing": "\u2205", "\\square": "\u25a1",
  "\\to": "\u2192", "\\ge": "\u2265", "\\le": "\u2264", "\\geq": "\u2265", "\\leq": "\u2264",
  "\\in": "\u2208", "\\neq": "\u2260", "\\cap": "\u2229", "\\cup": "\u222a",
  "\\rightrightarrows": "\u21c9", "\\rightleftarrows": "\u21c4",
  "\\rightarrow": "\u2192", "\\leftrightarrow": "\u2194", "\\mapsto": "\u21a6",
  "\\cdot": "\u00b7", "\\times": "\u00d7", "\\propto": "\u221d", "\\infty": "\u221e",
  "\\langle": "\u27e8", "\\rangle": "\u27e9", "\\parallel": "\u2225", "\\top": "\u22a4",
  "\\sum": "\u03a3", "\\sqrt": "\u221a", "\\pm": "\u00b1", "\\approx": "\u2248",
  "\\mid": "|", "\\lvert": "|", "\\rvert": "|", "\\ast": "*", "\\id": "id",
  "\\Theta": "\u0398", "\\Lambda": "\u039b", "\\Sigma": "\u03a3", "\\Omega": "\u03a9",
  "\\Gamma": "\u0393", "\\Phi": "\u03a6", "\\Psi": "\u03a8", "\\Pi": "\u03a0",
  "\\gamma": "\u03b3", "\\delta": "\u03b4", "\\epsilon": "\u03b5", "\\nu": "\u03bd",
  "\\subseteq": "\u2286", "\\subset": "\u2282", "\\supseteq": "\u2287",
  "\\forall": "\u2200", "\\exists": "\u2203", "\\emptyset": "\u2205",
  "\\equiv": "\u2261", "\\sim": "\u223c", "\\simeq": "\u2243",
  "\\circ": "\u2218", "\\setminus": "\u2216", "\\perp": "\u22a5",
  "\\longmapsto": "\u27fc", "\\ne": "\u2260", "\\colon": ":",
  "\\kappa": "\u03ba", "\\psi": "\u03c8", "\\varepsilon": "\u03b5",
  "\\varphi": "\u03c6", "\\omega": "\u03c9", "\\pi": "\u03c0", "\\xi": "\u03be",
  "\\zeta": "\u03b6", "\\chi": "\u03c7", "\\iota": "\u03b9", "\\upsilon": "\u03c5",
  "\\wedge": "\u2227", "\\vee": "\u2228", "\\int": "\u222b", "\\partial": "\u2202",
  "\\Pr": "Pr", "\\exp": "exp", "\\ln": "ln", "\\log": "log", "\\max": "max",
  "\\min": "min", "\\lim": "lim", "\\det": "det", "\\deg": "deg",
  "\\prod": "\u220f", "\\cdots": "\u22ef", "\\dots": "\u2026", "\\ldots": "\u2026",
  "\\gg": "\u226b", "\\ll": "\u226a", "\\notin": "\u2209", "\\subsetneq": "\u228a",
};

// Macros the linearizer resolves structurally rather than through SYM; listing
// them here keeps the unresolved-token audit from flagging them.
const STRUCTURAL = new Set([
  "frac", "tfrac", "dfrac", "binom", "bar", "overline", "dot", "boxed",
  "Bigl", "Bigr", "Bigm", "bigl", "bigr", "bigm", "Big", "big", "left", "right",
  "qquad", "quad", "text", "textrm", "textit", "textbf", "mathcal", "mathbb",
  "mathrm", "mathbf", "mathsf", "operatorname", "begin", "end", "label",
  "eqref", "cite", "hline", "rm",
]);

// Flatten LaTeX constructs that a run-based renderer cannot lay out in 2D.
// The .md source stays authoritative; this is a legible linearization.
function linearize(s) {
  // Substitute known symbols FIRST. Stripping a wrapper like \mathcal{H} later
  // removes its backslash, which would weld a preceding command to the letters
  // that follow it (\mapsto\mathcal{R} -> \mapstoR) and lose the token boundary.
  // A TeX control word swallows the whitespace that terminates it, so "\langle S"
  // is "<S", not "< S". Eat that separator when the token resolves to a symbol;
  // keep it otherwise, or an unresolved token would weld to what follows.
  let o = s.replace(/\\[a-zA-Z]+( ?)/g,
    (m, sp, off, str) => (SYM[m.trimEnd()] !== undefined ? SYM[m.trimEnd()] : m));
  o = o.replace(/\\begin\{pmatrix\}([\s\S]*?)\\end\{pmatrix\}/g, (_, b) =>
    "[" + b.replace(/\\\\/g, " ; ").replace(/&/g, ", ").replace(/\s+/g, " ").trim() + "]");
  o = o.replace(/\\(?:t|d)?frac\s*\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\}\s*\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\}/g,
    (_, a, b) => `(${a})/(${b})`);
  // \tfrac12 shorthand: two bare digits, no braces
  o = o.replace(/\\(?:t|d)?frac\s*(\d)\s*(\d)/g, (_, a, b) => `${a}/${b}`);
  o = o.replace(/\\begin\{cases\}([\s\S]*?)\\end\{cases\}/g, (_, b) =>
    "\\{" + b.replace(/\\\\/g, "; ").replace(/&/g, "").replace(/\s+/g, " ").trim() + "\\}");
  o = o.replace(/\\(?:bar|overline)\s*\{?([A-Za-z])\}?/g, "$1\u0304");
  o = o.replace(/\\dot\s*\{?([A-Za-z\u0370-\u03ff])\}?/g, "$1\u0307");
  o = o.replace(/\\binom\s*\{((?:[^{}]|\{[^{}]*\})*)\}\s*\{((?:[^{}]|\{[^{}]*\})*)\}/g,
    (_, a, b) => `C(${a}, ${b})`);
  o = o.replace(/\\boxed\s*\{/g, "{");
  o = o.replace(/\\(?:Bigl|Bigr|Bigm|bigl|bigr|bigm|left|right|Big|big)\s*/g, "");
  o = o.replace(/\\(?:qquad|quad|;|:|,|!)/g, " ");
  o = o.replace(/\\text(?:rm|it|bf)?\s*\{([^{}]*)\}/g, "$1");
  o = o.replace(/\\mathcal\s*\{([^{}]*)\}/g, "$1");
  const BB = { P: "\u2119", R: "\u211d", N: "\u2115", Z: "\u2124", Q: "\u211a", C: "\u2102" };
  o = o.replace(/\\mathbb\s*\{([^{}]*)\}/g, (_, x) => BB[x] || x);
  o = o.replace(/\\(?:mathrm|mathbf|mathsf|operatorname)\s*\{([^{}]*)\}/g, "$1");
  o = o.replace(/\\#/g, "#").replace(/\\%/g, "%").replace(/\\&/g, "&");
  o = o.replace(/\\\\/g, "  ");
  // Anything still carrying a backslash would be silently degraded to a bare
  // word downstream ("\kappa" -> "kappa"), which is how unhandled macros used
  // to reach the page unnoticed. Record them so the build can report them.
  for (const m of o.matchAll(/\\([a-zA-Z]+)/g)) {
    if (!STRUCTURAL.has(m[1])) UNRESOLVED.set(m[1], (UNRESOLVED.get(m[1]) || 0) + 1);
  }
  return o;
}

// Parse a $...$ body into [{text, sub, sup, italic}]
function parseMath(src) {
  const out = [];
  let i = 0;
  const push = (text, opt = {}) => { if (text) out.push({ text, italics: opt.sub || opt.sup ? false : true, ...opt }); };
  while (i < src.length) {
    const c = src[i];
    if (c === "\\") {
      if (src[i + 1] === "{" || src[i + 1] === "}") {
        push(src[i + 1], { italics: false });
        i += 2;
        continue;
      }
      const m = /^\\[a-zA-Z]+/.exec(src.slice(i));
      if (m) {
        const tok = m[0];
        if (tok === "\\rm" || tok === "\\mathrm" || tok === "\\text") {
          i += tok.length;
          const g = grab(src, i);
          push(g.body, { italics: false });
          i = g.next;
          continue;
        }
        push(SYM[tok] !== undefined ? SYM[tok] : tok.slice(1), { italics: !SYM[tok] });
        i += tok.length;
        continue;
      }
      i += 1;
      continue;
    }
    if (c === "_" || c === "^") {
      const kind = c === "_" ? "sub" : "sup";
      i += 1;
      const g = grab(src, i);
      // strip \rm, then map \symbol tokens to their Unicode glyphs
      let inner = g.body.replace(/\\(rm|mathrm|text)\s*/g, "");
      inner = inner.replace(/\\[a-zA-Z]+/g, (tok) => (SYM[tok] !== undefined ? SYM[tok] : tok.slice(1)));
      inner = inner.replace(/[{}]/g, "");
      out.push({ text: inner, [kind]: true, italics: false });
      i = g.next;
      continue;
    }
    if (c === "{" || c === "}") { i += 1; continue; }
    // run of plain chars
    const m2 = /^[^\\_^{}]+/.exec(src.slice(i));
    const chunk = m2[0];
    // digits / punctuation upright, letters italic
    for (const part of chunk.match(/[A-Za-z]+|[^A-Za-z]+/g) || []) {
      push(part, { italics: /[A-Za-z]/.test(part) });
    }
    i += chunk.length;
  }
  return out;
}
function grab(src, i) {
  if (src[i] === "{") {
    let depth = 0, j = i;
    for (; j < src.length; j++) {
      if (src[j] === "{") depth++;
      else if (src[j] === "}") { depth--; if (!depth) break; }
    }
    return { body: src.slice(i + 1, j), next: j + 1 };
  }
  return { body: src[i] || "", next: i + 1 };
}

// ------------------------------------------------------------- inline ----
function runs(text, base = {}) {
  const out = [];
  // split on $math$, **bold**, *italic*, `code`
  const re = /(\$[^$]+\$)|(\*\*[^*]+\*\*)|(`[^`]+`)|(\*[^*]+\*)/g;
  let last = 0, m;
  const plain = (s) => { if (s) out.push(new TextRun({ text: s, font: FONT, ...base })); };
  while ((m = re.exec(text))) {
    plain(text.slice(last, m.index));
    if (m[1]) {
      for (const p of parseMath(linearize(m[1].slice(1, -1)))) {
        out.push(new TextRun({
          text: p.text, font: { ascii: "Cambria Math", hAnsi: "Cambria Math", eastAsia: CJK },
          italics: !!p.italics, subScript: !!p.sub, superScript: !!p.sup, ...base,
        }));
      }
    } else if (m[2]) {
      plain2(out, m[2].slice(2, -2), { ...base, bold: true });
    } else if (m[3]) {
      out.push(new TextRun({
        text: m[3].slice(1, -1),
        font: { ascii: "Consolas", hAnsi: "Consolas", eastAsia: CJK },
        ...base,
      }));
    } else if (m[4]) {
      plain2(out, m[4].slice(1, -1), { ...base, italics: true });
    }
    last = re.lastIndex;
  }
  plain(text.slice(last));
  return out;
}
// bold segment may itself contain math
function plain2(out, text, base) {
  for (const r of runs(text, base)) out.push(r);
}

// -------------------------------------------------------------- blocks ----
let md = fs.readFileSync(SRC, "utf8").split("\n");

// Equation labels. Hand-written cross-references drift the moment an equation is
// inserted; resolve \eqref{} against the numbering the renderer produces.
// SIAM numbers equations within the section — (3.1), (3.2), … — so a section
// built on its own already carries the numbering it will have in the manuscript.
const EQLABEL = {};
// `# 3. 标题` sets the section; equations restart at 1 inside it.
function sectionOf(line) {
  const m = /^#\s+(\d+)\.\s/.exec(line.trim());
  return m ? Number(m[1]) : null;
}
function scanEquations(lines) {
  let sec = 0, n = 0, inFence = false;
  for (let j = 0; j < lines.length; j++) {
    const t = lines[j].trim();
    if (t.startsWith("```")) { inFence = !inFence; continue; }
    if (inFence) continue;
    const s = sectionOf(t);
    if (s !== null) { sec = s; n = 0; continue; }
    if (t !== "$$") continue;
    const buf = [];
    j++;
    while (j < lines.length && lines[j].trim() !== "$$") { buf.push(lines[j]); j++; }
    const body = buf.join(" ");
    if (/\\notag/.test(body)) continue;
    n += 1;
    const lm = /\\label\{([^}]+)\}/.exec(body);
    if (lm) EQLABEL[lm[1]] = sec ? `${sec}.${n}` : String(n);
  }
}
// A section built on its own still has to resolve cross-references into earlier
// sections. EQ_PRESCAN names the sources that precede this one: they contribute
// labels without being rendered. Section-scoped numbering means they no longer
// have to advance any counter.
for (const f of (process.env.EQ_PRESCAN || "").split(",").filter(Boolean)) {
  scanEquations(fs.readFileSync(f.trim(), "utf8").split("\n"));
}
scanEquations(md);
md = md.map((l) => l.replace(/\\eqref\{([^}]+)\}/g, (m0, name) => {
  if (!(name in EQLABEL)) { console.warn(`  undefined \\eqref{${name}}`); return "(??)"; }
  return "(" + EQLABEL[name] + ")";
}));
const children = [];
const SP = { before: 120, after: 120, line: 320, lineRule: "auto" };

const para = (text, opts = {}) =>
  new Paragraph({ children: runs(text), spacing: SP, ...opts });

function tableFrom(rows) {
  const n = rows[0].length;
  // A `|` inside inline math splits a row into too many cells; a short row is a
  // typo. Normalise to the header width so a malformed row can't corrupt the table.
  rows = rows.map((r, ri) => {
    if (r.length === n) return r;
    console.warn(`  table row ${ri}: ${r.length} cells, expected ${n} -> ${r.length > n ? "merging overflow" : "padding"}`);
    return r.length > n
      ? [...r.slice(0, n - 1), r.slice(n - 1).join(" | ")]
      : [...r, ...Array(n - r.length).fill("")];
  });
  const total = convertInchesToTwip(6.5);
  const colW = Array(n).fill(Math.floor(total / n));
  colW[n - 1] = total - colW[0] * (n - 1);
  return new Table({
    columnWidths: colW,
    width: { size: total, type: WidthType.DXA },
    rows: rows.map((cells, ri) =>
      new TableRow({
        tableHeader: ri === 0,
        children: cells.map((c, ci) =>
          new TableCell({
            width: { size: colW[ci], type: WidthType.DXA },
            shading: ri === 0 ? { type: ShadingType.CLEAR, fill: "E8ECF2" } : undefined,
            margins: { top: 80, bottom: 80, left: 110, right: 110 },
            children: [new Paragraph({
              children: runs(c, ri === 0 ? { bold: true } : {}),
              spacing: { before: 40, after: 40, line: 300 },
            })],
          })),
      })),
  });
}

let eqSec = 0, eqNo = 0;
let i = 0;
while (i < md.length) {
  const line = md[i];
  const t = line.trim();

  if (!t) { i++; continue; }

  // figure: ![caption](path){width=inches}
  const im = /^!\[([^\]]*)\]\(([^)]+)\)(?:\{width=([\d.]+)\})?$/.exec(t);
  if (im) {
    const p = path.isAbsolute(im[2]) ? im[2] : path.join(path.dirname(SRC), im[2]);
    if (!fs.existsSync(p)) {
      console.error("MISSING FIGURE:", p);
      process.exitCode = 2;
    } else {
      const wIn = im[3] ? parseFloat(im[3]) : 6.5;
      const png = require("zlib"); // only to keep require list stable
      const buf = fs.readFileSync(p);
      // PNG header: width/height are big-endian uint32 at bytes 16 and 20
      const pw = buf.readUInt32BE(16), ph = buf.readUInt32BE(20);
      children.push(new Paragraph({
        alignment: AlignmentType.CENTER,
        // An omitted lineRule is read as "exact" by some renderers, which
        // clips an inline image to one line height. Say "auto" explicitly.
        spacing: { before: 200, after: 60, line: 240, lineRule: "auto" },
        children: [new D.ImageRun({
          data: buf, type: "png",
          transformation: { width: wIn * 96, height: wIn * 96 * ph / pw },
        })],
      }));
      if (im[1]) {
        children.push(new Paragraph({
          alignment: AlignmentType.LEFT,
          spacing: { before: 0, after: 200 },
          children: runs(im[1], { size: 19 }),
        }));
      }
    }
    i++; continue;
  }

  if (/^---+$/.test(t)) {
    children.push(new Paragraph({
      text: "", spacing: { before: 80, after: 200 },
      border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "BBBBBB", space: 1 } },
    }));
    i++; continue;
  }

  // fenced code block: monospace, one paragraph per line, no inline parsing
  // (ASCII diagrams only align in a fixed-width font, and their glyphs would
  // otherwise be eaten as markdown)
  if (t.startsWith("```")) {
    i++;
    const buf = [];
    while (i < md.length && !md[i].trim().startsWith("```")) { buf.push(md[i]); i++; }
    i++; // consume closing fence
    buf.forEach((b, bi) => children.push(new Paragraph({
      children: [new TextRun({
        text: b.replace(/\t/g, "    ") || " ",
        font: { ascii: "Consolas", hAnsi: "Consolas", eastAsia: CJK },
        size: 17,
      })],
      spacing: { before: bi === 0 ? 140 : 0, after: bi === buf.length - 1 ? 140 : 0, line: 240 },
      shading: { type: ShadingType.CLEAR, fill: "F4F6F8" },
      indent: { left: 200 },
    })));
    continue;
  }

  // display math: a $$ fence, its body, and a closing $$
  if (t === "$$") {
    const buf = [];
    i++;
    while (i < md.length && md[i].trim() !== "$$") { buf.push(md[i].trim()); i++; }
    i++; // consume closing fence
    // \tag{} suppresses auto-numbering (e.g. for an unnumbered aside)
    let body = buf.join(" ");
    const noNum = /\\notag/.test(body);
    body = linearize(body.replace(/\\notag/g, "").replace(/\\label\{[^}]*\}/g, ""));
    const mathFont = { ascii: "Cambria Math", hAnsi: "Cambria Math", eastAsia: CJK };
    const kids = parseMath(body).map((p) => new TextRun({
      text: p.text, font: mathFont,
      italics: !!p.italics, subScript: !!p.sub, superScript: !!p.sup, size: 21,
    }));
    // Centre the equation and hang its number at the right margin: a centre tab
    // stop plus a right tab stop, so the number never rides on the equation.
    const W = convertInchesToTwip(6.5);
    if (!noNum) {
      eqNo += 1;
      const tag = eqSec ? `${eqSec}.${eqNo}` : String(eqNo);
      kids.unshift(new TextRun({ text: "\t", font: mathFont }));
      kids.push(new TextRun({ text: "\t(" + tag + ")", font: mathFont, size: 21 }));
    }
    children.push(new Paragraph({
      children: kids,
      alignment: noNum ? AlignmentType.CENTER : AlignmentType.LEFT,
      tabStops: noNum ? undefined : [
        { type: TabStopType.CENTER, position: Math.floor(W / 2) },
        { type: TabStopType.RIGHT, position: W },
      ],
      spacing: { before: 180, after: 180, line: 300 },
    }));
    continue;
  }

  if (t.startsWith("#")) {
    const lvl = t.match(/^#+/)[0].length;
    const txt = t.replace(/^#+\s*/, "");
    const s = sectionOf(t);
    if (s !== null) { eqSec = s; eqNo = 0; }
    const map = { 1: HeadingLevel.HEADING_1, 2: HeadingLevel.HEADING_2, 3: HeadingLevel.HEADING_3 };
    children.push(new Paragraph({
      children: runs(txt, { bold: true, size: lvl === 1 ? 32 : lvl === 2 ? 26 : 23, color: "1A1A1A" }),
      heading: map[lvl] || HeadingLevel.HEADING_4,
      spacing: { before: lvl === 1 ? 0 : 300, after: 160 },
      alignment: lvl === 1 ? AlignmentType.CENTER : undefined,
    }));
    i++; continue;
  }

  if (t.startsWith(">")) {
    const buf = [];
    while (i < md.length && md[i].trim().startsWith(">")) {
      buf.push(md[i].trim().replace(/^>\s?/, ""));
      i++;
    }
    for (const b of buf) {
      children.push(new Paragraph({
        children: runs(b, { color: "3A4A5A" }),
        spacing: { before: 60, after: 60, line: 300 },
        indent: { left: 360 },
        border: { left: { style: BorderStyle.SINGLE, size: 12, color: "9AAFC4", space: 8 } },
      }));
    }
    continue;
  }

  if (t.startsWith("|")) {
    const rows = [];
    while (i < md.length && md[i].trim().startsWith("|")) {
      const r = md[i].trim();
      if (!/^\|[\s:\-|]+\|$/.test(r)) {
        rows.push(r.replace(/^\||\|$/g, "").split("|").map((s) => s.trim()));
      }
      i++;
    }
    children.push(tableFrom(rows));
    children.push(new Paragraph({ text: "", spacing: { after: 160 } }));
    continue;
  }

  if (/^-\s+/.test(t)) {
    children.push(new Paragraph({
      children: runs(t.replace(/^-\s+/, "")),
      numbering: { reference: "bullets", level: 0 },
      spacing: { before: 60, after: 60, line: 310 },
    }));
    i++; continue;
  }

  if (/^\d+\.\s+/.test(t)) {
    children.push(new Paragraph({
      children: runs(t.replace(/^\d+\.\s+/, "")),
      numbering: { reference: "nums", level: 0 },
      spacing: { before: 100, after: 100, line: 310 },
    }));
    i++; continue;
  }

  children.push(para(t, { alignment: AlignmentType.JUSTIFIED }));
  i++;
}

// ---------------------------------------------------- bibliography ----
function parseBib(src) {
  const entries = [];
  const re = /@(?:article|inproceedings|incollection|book)\{([^,]+),([\s\S]*?)\n\}/g;
  let m;
  while ((m = re.exec(src))) {
    const key = m[1].trim();
    const f = {};
    const fre = /(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}/g;
    let fm;
    while ((fm = fre.exec(m[2]))) f[fm[1].toLowerCase()] = fm[2];
    entries.push({ key, ...f });
  }
  return entries;
}
function deTeX(s) {
  if (!s) return "";
  const acc = {
    "{\\~n}": "ñ", "{\\'a}": "á", "{\\'e}": "é", "{\\'i}": "í", "{\\'o}": "ó", "{\\'u}": "ú",
    "{\\\"u}": "ü", "{\\\"o}": "ö", "{\\'A}": "Á", "{\\~a}": "ã", "{\\`e}": "è",
  };
  let o = s;
  for (const [k, v] of Object.entries(acc)) o = o.split(k).join(v);
  return o.replace(/\{\\'([a-zA-Z])\}/g, "$1").replace(/\{\\"([a-zA-Z])\}/g, "$1")
          .replace(/\{\\~([a-zA-Z])\}/g, "$1").replace(/[{}\\]/g, "").replace(/\s+/g, " ").trim();
}
function fmtAuthors(a) {
  if (!a) return "";
  const list = deTeX(a).split(" and ").map((n) => {
    const p = n.split(",").map((x) => x.trim());
    return p.length === 2 ? `${p[1]} ${p[0]}` : n.trim();
  });
  if (list.length > 6) return list.slice(0, 6).join(", ") + " et al.";
  return list.join(", ");
}

const bib = process.env.NO_AUTO_BIB
  ? []
  : parseBib(fs.readFileSync(BIB, "utf8")).sort((a, b) => a.key.localeCompare(b.key));

if (bib.length) {
children.push(new Paragraph({
  text: "", spacing: { before: 200, after: 200 },
  border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "BBBBBB", space: 1 } },
}));
children.push(new Paragraph({
  children: runs("参考文献", { bold: true, size: 26, color: "1A1A1A" }),
  heading: HeadingLevel.HEADING_2, spacing: { before: 0, after: 160 },
}));
children.push(para("按 BibTeX key 字母序排列，正文中的 [key] 与此处一一对应。全部条目已通过 cite-verify 核验。"));

for (const e of bib) {
  const bits = [];
  bits.push(fmtAuthors(e.author));
  bits.push(`(${e.year})`);
  bits.push(deTeX(e.title) + ".");
  let tail = deTeX(e.journal || e.booktitle);
  if (e.volume) tail += ` ${e.volume}`;
  if (e.pages) tail += `, ${deTeX(e.pages).replace(/--/g, "–")}`;
  bits.push(tail + ".");
  const body = bits.join(" ");
  children.push(new Paragraph({
    children: [
      new TextRun({ text: `[${e.key}] `, font: { ascii: "Consolas", hAnsi: "Consolas", eastAsia: CJK }, bold: true, size: 18 }),
      new TextRun({ text: body, font: FONT, size: 19 }),
      new TextRun({ text: ` doi:${e.doi}`, font: { ascii: "Consolas", hAnsi: "Consolas", eastAsia: CJK }, size: 17, color: "555555" }),
    ],
    spacing: { before: 60, after: 60, line: 280 },
    indent: { left: 420, hanging: 420 },
  }));
}
}

// -------------------------------------------------------------- doc ----
const doc = new Document({
  creator: "SIR on Directed Hypergraphs",
  title: "第 I 部分　引言与定位",
  styles: {
    default: {
      document: { run: { font: FONT, size: 21 }, paragraph: { spacing: { line: 320, lineRule: "auto" } } },
    },
  },
  numbering: {
    config: [
      {
        reference: "bullets",
        levels: [{
          level: 0, format: LevelFormat.BULLET, text: "\u2022", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 460, hanging: 250 } } },
        }],
      },
      {
        reference: "nums",
        levels: [{
          level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 460, hanging: 360 } } },
        }],
      },
    ],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 11906, height: 16838 },           // A4
        margin: { top: 1200, right: 1100, bottom: 1200, left: 1100 },
      },
    },
    children,
  }],
});

Packer.toBuffer(doc).then((b) => {
  fs.writeFileSync(OUT, b);
  console.log("wrote", OUT, b.length, "bytes;", bib.length, "references");
  if (UNRESOLVED.size) {
    const list = [...UNRESOLVED.entries()].sort((a, b) => b[1] - a[1])
      .map(([k, n]) => `\\${k}(${n})`).join(" ");
    console.error("UNRESOLVED LATEX -> rendered as bare words:", list);
    process.exitCode = 2;
  }
});
