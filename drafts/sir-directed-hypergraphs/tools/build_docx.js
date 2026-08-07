const fs = require("fs");
const path = require("path");
const D = require("docx");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  LevelFormat, convertInchesToTwip,
} = D;

const SRC = "/home/user/Claude-skills/drafts/sir-directed-hypergraphs/part1-introduction.md";
const BIB = "/home/user/Claude-skills/drafts/sir-directed-hypergraphs/references.bib";
const OUT = "/home/user/Claude-skills/drafts/sir-directed-hypergraphs/part1-introduction.docx";

const LATIN = "Times New Roman";
const CJK = "SimSun";
const FONT = { ascii: LATIN, hAnsi: LATIN, eastAsia: CJK, cs: LATIN };

// ---------------------------------------------------------------- math ----
const SYM = {
  "\\phi": "\u03c6", "\\alpha": "\u03b1", "\\tau": "\u03c4", "\\eta": "\u03b7",
  "\\mu": "\u03bc", "\\theta": "\u03b8", "\\lambda": "\u03bb", "\\beta": "\u03b2",
  "\\to": "\u2192", "\\ge": "\u2265", "\\le": "\u2264", "\\in": "\u2208",
  "\\rightrightarrows": "\u21c9", "\\rightleftarrows": "\u21c4",
  "\\rightarrow": "\u2192", "\\cdot": "\u00b7", "\\times": "\u00d7",
};

// Parse a $...$ body into [{text, sub, sup, italic}]
function parseMath(src) {
  const out = [];
  let i = 0;
  const push = (text, opt = {}) => { if (text) out.push({ text, italics: opt.sub || opt.sup ? false : true, ...opt }); };
  while (i < src.length) {
    const c = src[i];
    if (c === "\\") {
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
  // split on $math$, **bold**, `code`
  const re = /(\$[^$]+\$)|(\*\*[^*]+\*\*)|(`[^`]+`)/g;
  let last = 0, m;
  const plain = (s) => { if (s) out.push(new TextRun({ text: s, font: FONT, ...base })); };
  while ((m = re.exec(text))) {
    plain(text.slice(last, m.index));
    if (m[1]) {
      for (const p of parseMath(m[1].slice(1, -1))) {
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
const md = fs.readFileSync(SRC, "utf8").split("\n");
const children = [];
const SP = { before: 120, after: 120, line: 320 };

const para = (text, opts = {}) =>
  new Paragraph({ children: runs(text), spacing: SP, ...opts });

function tableFrom(rows) {
  const header = rows[0];
  const n = header.length;
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

let i = 0;
while (i < md.length) {
  const line = md[i];
  const t = line.trim();

  if (!t) { i++; continue; }

  if (/^---+$/.test(t)) {
    children.push(new Paragraph({
      text: "", spacing: { before: 80, after: 200 },
      border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "BBBBBB", space: 1 } },
    }));
    i++; continue;
  }

  if (t.startsWith("#")) {
    const lvl = t.match(/^#+/)[0].length;
    const txt = t.replace(/^#+\s*/, "");
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
  const re = /@article\{([^,]+),([\s\S]*?)\n\}/g;
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

const bib = parseBib(fs.readFileSync(BIB, "utf8")).sort((a, b) => a.key.localeCompare(b.key));

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
  let tail = deTeX(e.journal);
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

// -------------------------------------------------------------- doc ----
const doc = new Document({
  creator: "SIR on Directed Hypergraphs",
  title: "第 I 部分　引言与定位",
  styles: {
    default: {
      document: { run: { font: FONT, size: 21 }, paragraph: { spacing: { line: 320 } } },
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
});
