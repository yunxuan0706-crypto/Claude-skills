// Combined manuscript: Part 1 (message-passing group closure + time evolution)
// followed by Part 2 (subcritical determination of the threshold + the rho12
// scan). Both sections are built by their own generators, which now export
// their content arrays; this file concatenates them under one title and page
// setup so the whole thing is one document. Regenerate the two figures and run
// the section generators first if their inputs changed.
const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, AlignmentType } = require("docx");

const part1 = require("./make_doc_cn_sec1.js").children;
const part2 = require("./make_doc_cn3.js").children;
const part3 = require("./make_doc_cn_sec3.js").children;

const CJK = "SimSun", HEI = "SimHei", INK = "000000", BODY = 21;

// paper title (added for the merged document; edit freely)
const title = new Paragraph({
  alignment: AlignmentType.CENTER,
  spacing: { before: 0, after: 340 },
  children: [new TextRun({
    text: "多重超图上的 SIR：消息传递群体闭合与爆发阈值",
    font: HEI, bold: true, size: 32, color: INK })],
});

const doc = new Document({
  styles: { default: { document: { run: { font: CJK, size: BODY, color: INK } } } },
  sections: [{
    properties: { page: { margin: { top: 1400, bottom: 1400, left: 1440, right: 1440 } } },
    children: [title, ...part1, ...part2, ...part3],
  }],
});
Packer.toBuffer(doc).then(b => {
  fs.writeFileSync("combined_cn.docx", b);
  console.log("wrote combined_cn.docx", b.length,
              "(", part1.length, "+", part2.length, "+", part3.length,
              "paragraphs + title )");
});
