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

const CJK = "SimSun", HEI = "SimHei", INK = "000000", BODY = 21, SMALL = 18;

const title = new Paragraph({
  alignment: AlignmentType.CENTER,
  spacing: { before: 0, after: 340 },
  children: [new TextRun({
    text: "多重超图上的 SIR：消息传递群体闭合与爆发阈值",
    font: HEI, bold: true, size: 32, color: INK })],
});

const R = (t, o = {}) => new TextRun({
  text: t, font: o.f || CJK, size: o.size || BODY, color: INK,
  bold: o.bold, italics: o.italics,
});
const refH = new Paragraph({ spacing: { before: 400, after: 200 }, indent: { firstLine: 0 },
  children: [new TextRun({ text: "参考文献", font: HEI, bold: true, size: 26, color: INK })] });

const refs = [
  "[1] Newman M E J. Spread of epidemic disease on networks. Phys Rev E, 2002, 66: 016128.",
  "[2] Newman M E J. Random graphs with clustering. Phys Rev Lett, 2009, 103: 058701.",
  "[3] Iacopini I, Petri G, Barrat A, Latora V. Simplicial models of social contagion. Nat Commun, 2019, 10: 2485.",
  "[4] de Arruda G F, Petri G, Moreno Y. Social contagion models on hypergraphs. Phys Rev Research, 2020, 2: 023032.",
  "[5] Boccaletti S, Bianconi G, Criado R, et al. The structure and dynamics of multilayer networks. Phys Rep, 2014, 544: 1–122.",
  "[6] De Domenico M, Granell C, Porter M A, Arenas A. The physics of spreading processes in multilayer networks. Nat Phys, 2016, 12: 901–906.",
  "[7] Gómez S, Arenas A, Borge-Holthoefer J, Meloni S, Moreno Y. Discrete-time Markov chain approach to contact-based disease spreading in complex networks. Europhys Lett, 2010, 89: 38009.",
  "[8] Mata A S, Ferreira S C. Pair quenched mean-field theory for the susceptible-infected-susceptible model on complex networks. Europhys Lett, 2013, 103: 48003.",
  "[9] Ding L, Zhu Y. Epidemic dynamics in higher-order networks with directed edges. Commun Nonlinear Sci Numer Simul, 2026.",
  "[10] Wang W, Tang M, Stanley H E, Braunstein L A. Unification of theoretical approaches for epidemic spreading on complex networks. Rep Prog Phys, 2017, 80: 036603.",
  "[11] Allard A, Althouse B M, Scarpino S V, Hébert-Dufresne L. The role of directionality, heterogeneity, and correlations in epidemic risk and spread. SIAM Rev, 2023, 65: 471–492.",
  "[12] Mézard M, Parisi G. The Bethe lattice spin glass revisited. Eur Phys J B, 2001, 20: 217–233.",
  "[13] Karrer B, Newman M E J. Message passing approach for general epidemic models. Phys Rev E, 2010, 82: 016101.",
  "[14] Diekmann O, Heesterbeek J A P, Metz J A J. On the definition and the computation of the basic reproduction ratio R₀ in models for infectious diseases in heterogeneous populations. J Math Biol, 1990, 28: 365–382.",
  "[15] Diekmann O, Heesterbeek J A P, Roberts M G. The construction of next-generation matrices for compartmental epidemic models. J R Soc Interface, 2010, 7: 873–885.",
  "[16] Gillespie D T. Exact stochastic simulation of coupled chemical reactions. J Phys Chem, 1977, 81: 2340–2361.",
  "[17] Fisher M E, Barber M N. Scaling theory for finite-size effects in the critical region. Phys Rev Lett, 1972, 28: 1516–1519.",
  "[18] Richardson L F. The approximate arithmetical solution by finite differences of physical problems. Phil Trans R Soc A, 1911, 210: 307–357.",
  "[19] Miller J C. Percolation and epidemics in random clustered networks. Phys Rev E, 2009, 80: 020901(R).",
];
const refParagraphs = refs.map(r => new Paragraph({
  children: [new TextRun({ text: r, font: CJK, size: SMALL, color: INK })],
  alignment: AlignmentType.BOTH,
  spacing: { after: 60, line: 280 },
  indent: { firstLine: 0, left: 420, hanging: 420 },
}));

const doc = new Document({
  styles: { default: { document: { run: { font: CJK, size: BODY, color: INK } } } },
  sections: [{
    properties: { page: { margin: { top: 1400, bottom: 1400, left: 1440, right: 1440 } } },
    children: [title, ...part1, ...part2, ...part3, refH, ...refParagraphs],
  }],
});
Packer.toBuffer(doc).then(b => {
  fs.writeFileSync("combined_cn.docx", b);
  console.log("wrote combined_cn.docx", b.length,
              "(", part1.length, "+", part2.length, "+", part3.length,
              "paragraphs + title )");
});
