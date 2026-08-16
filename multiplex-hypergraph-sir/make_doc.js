const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, BorderStyle, ShadingType, ImageRun,
} = require("docx");

const INK = "20201E", ACCENT = "1C6E64", HEAD = "2F5D57", GREY = "6B6A63";
const HDR = "E7F1EF", ZEBRA = "F6F5F1";
const FONT = "Calibri";

// ---------- helpers ----------
const R = (t, o = {}) => new TextRun({ text: t, font: FONT, size: o.size || 20, color: o.color || INK, bold: o.bold, italics: o.italics, subScript: o.sub, superScript: o.sup });
const P = (children, o = {}) => new Paragraph({ children: Array.isArray(children) ? children : [children], spacing: { after: o.after == null ? 120 : o.after, before: o.before || 0, line: 276 }, alignment: o.align });
const H = (t, lvl) => new Paragraph({ heading: lvl, spacing: { before: 240, after: 100 }, children: [new TextRun({ text: t, font: FONT, bold: true, color: HEAD, size: lvl === HeadingLevel.HEADING_1 ? 28 : 24 })] });

function cell(runs, { w, header, zebra, align } = {}) {
  return new TableCell({
    width: { size: w, type: WidthType.DXA },
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    shading: header ? { type: ShadingType.CLEAR, fill: HDR } : (zebra ? { type: ShadingType.CLEAR, fill: ZEBRA } : undefined),
    children: [new Paragraph({ alignment: align || AlignmentType.LEFT, spacing: { after: 0, line: 252 },
      children: (Array.isArray(runs) ? runs : [runs]).map(t => typeof t === "string" ? new TextRun({ text: t, font: FONT, size: 18, bold: header, color: header ? HEAD : INK }) : t) })],
  });
}
function table(widths, rows) {
  const total = widths.reduce((a, b) => a + b, 0);
  return new Table({
    columnWidths: widths, width: { size: total, type: WidthType.DXA },
    borders: {
      top: { style: BorderStyle.SINGLE, size: 4, color: "BFC6C4" },
      bottom: { style: BorderStyle.SINGLE, size: 4, color: "BFC6C4" },
      left: { style: BorderStyle.NONE }, right: { style: BorderStyle.NONE },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 2, color: "E2E2DC" },
      insideVertical: { style: BorderStyle.NONE },
    },
    rows: rows.map((cells, i) => new TableRow({
      tableHeader: i === 0, cantSplit: true,
      children: cells.map((c, j) => cell(c, { w: widths[j], header: i === 0, zebra: i > 0 && i % 2 === 0, align: j === 0 ? AlignmentType.LEFT : AlignmentType.CENTER })),
    })),
  });
}
const rule = () => new Paragraph({ spacing: { before: 60, after: 120 }, border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "D8D8D2" } }, children: [] });

// math run shorthands
const sub = (b, s) => [R(b), R(s, { sub: true })];
const rho = "ρ", lam = "λ", eps = "ε", sig = "σ", chi = "χ", inf = "∞", le = "≤", times = "×", to = "→", deg = "°";

// ---------- content ----------
const children = [];

children.push(new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 40 },
  children: [new TextRun({ text: "由亚临界终态外推测定 " + lam + "c ", font: FONT, bold: true, size: 34, color: INK }),
             new TextRun({ text: "（数据、结果与图 2 说明）", font: FONT, bold: true, size: 26, color: HEAD })] }));
children.push(new Paragraph({ spacing: { after: 160 }, children: [new TextRun({ text: "多重超图上的 SIR 传播 · 群体闭合理论的第三条独立阈值复核", font: FONT, italics: true, size: 20, color: GREY })] }));
children.push(rule());

// 1. 方法
children.push(H("1  方法概述", HeadingLevel.HEADING_1));
children.push(P([R("阈值 "), ...sub(lam, "c"), R(" 的第三条独立复核不经由线性化，而由亚临界终态的规模发散给出。取一个小的初始感染比例 "), R(eps), R("（每个节点独立以概率 "), R(eps), R(" 置为感染态），在 "), R(lam), R(" < "), ...sub(lam, "c"), R(" 处将群体闭合方程组 (2.2)、(2.4)–(2.6) 积分至终态，由 I("), R(inf), R(")=0 得最终规模 R("), R(inf), R(") = 1 − (1−"), R(eps), R(")·"), R("Ψ"), R("("), R("Φ"), R("("), R(inf), R(")).")]));
children.push(P([R("亚临界簇的平均总规模为 "), R(chi), R("("), R(lam), R(") = "), R("\u{1D7D9}"), R("ᵀ", { sup: true }), R("(I−N)"), R("−1", { sup: true }), ...sub("v", "0"), R("，在 "), R(rho), R("(N)"), R(to), R("1"), R("−", { sup: true }), R(" 时按 (1−"), R(rho), R(")"), R("−1", { sup: true }), R(" 发散；因而其倒数 "), R(eps), R("/R("), R(inf), R(") = 1/"), R(chi), R(" 随 "), R(lam), R(" 线性趋零，零点即 "), ...sub(lam, "c"), R("。有限种子的偏置由 "), R(eps), R(to), R("0 的 Richardson 外推消去。")]));
children.push(P([R("该路径与 "), R(rho), R("(N)=1 的特征值计算"), R("不共享代码", { bold: true }), R("：一侧是非线性方程组积分至终态再作线性拟合，另一侧是次代矩阵 N 的谱半径二分。两者仅共享模型定义，故构成一次真正独立的比对。")]));

// 2. 理论正确性校验
children.push(H("2  理论正确性校验", HeadingLevel.HEADING_1));
children.push(P([R("在引入结果之前先确立理论自身的可信度。下表汇总解析锚点与闭合方程组的内部自检，全部通过。", {})]));
children.push(table([3100, 2400, 3060, 1000], [
  ["校验项", "期望", "计算结果", "状态"],
  [[R("C(2,"+lam+",1) = "+lam+"/(1+"+lam+")", { size: 18 })], "精确", [R(lam+"=0.13/0.35/1/2 均到 10", { size: 18 }), R("−12", { sup: true, size: 18 })], "通过"],
  [[R("C(3,4,5; "+lam+"=1) 超出 (m−1)T", { size: 18 })], "11 / 22 / 31 %", "11.11 / 21.53 / 31.05 %", "通过"],
  [[R("单种子激活 "+"θ"+"≥2 的群体", { size: 18 })], "C = 0", "0", "通过"],
  [[R("Anchor A: P={(2,2),(3,3)}, m=(3,3), "+lam+"=0.13", { size: 17 })], [R(rho+"(N) = 1.013", { size: 18 })], [R("N", { size: 18 }), R("11", { sub: true, size: 18 }), R("=0.3860, "+rho+"=1.0132", { size: 18 })], "通过"],
  [[R("成对退化 5-正则单层", { size: 18 })], [...sub(lam, "c").map(r => r), R(" = 1/3", { size: 18 })], "0.3333333333", "通过"],
  [[R("Anchor B: "+"θ"+"=2 通道出列", { size: 18 })], "阈值不变", [R("C=0，精确相等", { size: 18 })], "通过"],
  [[R("Anchor B: 大纲 0.410619", { size: 18 })], [R(rho+"值 × (1+2.5%)", { size: 18 })], [R("0.400567 × 1.0251 = +2.51%", { size: 18 })], "通过"],
  [[R("求和规则 (2.6)", { size: 18 })], "≈ 0", [R("残差 ~ 10", { size: 18 }), R("−17", { sup: true, size: 18 })], "通过"],
  [[R("全易感恒等式 (2.7)", { size: 18 })], "≈ 0", [R("残差 ~ 10", { size: 18 }), R("−12", { sup: true, size: 18 })], "通过"],
  [[R("无病定态 Jacobian 谱横标", { size: 18 })], [R("0  at  ", { size: 18 }), ...sub(lam, "c")], "0（机器精度）", "通过"],
]));
children.push(P([R("要点：大纲 2.9 节的 0.410619 并非 ", { size: 18 }), R(rho, { size: 18 }), R("(N)=1 的值，而是有限 t", { size: 18 }), R("max", { sub: true, size: 18 }), R(" 的终态二分阈值——本文 ", { size: 18 }), R(rho, { size: 18 }), R("(N)=1 给 0.400567，二者比值恰为 1.0251，正是 3.1 节记录的“积分窗口偏高 2.5–2.6%”。", { size: 18, italics: true })], { after: 60 }));

// 3. 主结果
children.push(H("3  主结果：八组构型的 " + lam + "c 外推", HeadingLevel.HEADING_1));
children.push(P([R("窗口取 "), R(lam), R(" ∈ [0.90, 0.995]·"), ...sub(lam, "c"), R(" 的 10 个点，"), R(eps), R(to), R("0 后作线性外推，x 截距即外推 "), ...sub(lam, "c"), R("。跨越 "), ...sub(lam, "c"), R(" ∈ [0.08, 0.40] 的八组构型（单层高阶超图、成对退化网络、双层多重超图）结果如下。")]));
const rows3 = [
  ["构型", [R("λᴄ [ρ(N)=1]", { size: 18 })], [R("λᴄ [外推]", { size: 18 })], "相对偏差", [R("σ", { size: 18 })], [R("偏离 (σ)", { size: 18 })]],
  ["2-layer  m=(4,4)", "0.080302", "0.080354", "+0.065%", [R("3.4×10", { size: 18 }), R("−5", { sup: true, size: 18 })], "+1.51"],
  ["2-layer  m=(3,4)", "0.099309", "0.099370", "+0.062%", [R("4.4×10", { size: 18 }), R("−5", { sup: true, size: 18 })], "+1.38"],
  ["2-layer  m=(3,3)", "0.128162", "0.128236", "+0.058%", [R("5.9×10", { size: 18 }), R("−5", { sup: true, size: 18 })], "+1.26"],
  ["3-reg.  m=4", "0.174483", "0.174536", "+0.031%", [R("9.8×10", { size: 18 }), R("−5", { sup: true, size: 18 })], "+0.55"],
  ["2-layer  m=(2,3)", "0.185417", "0.185512", "+0.052%", [R("9.3×10", { size: 18 }), R("−5", { sup: true, size: 18 })], "+1.03"],
  ["6-reg.  pairwise", "0.250000", "0.250096", "+0.039%", [R("1.3×10", { size: 18 }), R("−4", { sup: true, size: 18 })], "+0.76"],
  ["5-reg.  pairwise", "0.333333", "0.333440", "+0.032%", [R("1.9×10", { size: 18 }), R("−4", { sup: true, size: 18 })], "+0.57"],
  ["1-layer  m=3", "0.400567", "0.400702", "+0.034%", [R("2.9×10", { size: 18 }), R("−4", { sup: true, size: 18 })], "+0.47"],
];
children.push(table([2260, 1560, 1560, 1300, 1340, 1140], rows3));
children.push(P([R("最大相对偏差 6.5", { bold: true }), R(times, { bold: true }), R("10", { bold: true }), R("−4", { sup: true, bold: true }), R("，八组全部落在拟合不确定度的 1.51", { bold: true }), R(sig, { bold: true }), R(" 之内。", { bold: true }), R(" 残余为已知的有限窗口效应——1/"), R(chi), R(" 仅在 "), ...sub(lam, "c"), R(" 邻域严格线性，窗口凸性使外推零点略偏高，其符号随窗口位置翻转、量级与 "), R(sig), R(" 相当，故不可与真实偏离区分。")], { after: 60 }));

// 4. 独立交叉验证
children.push(H("4  独立交叉验证", HeadingLevel.HEADING_1));
children.push(table([3050, 3050, 3460], [
  ["检验", "方法", "结果"],
  [[R("闭合终态 = 分支公式", { size: 18 })], [R("χ_ODE vs 1+\u{1D7D9}ᵀ(I−Nᵀ)⁻¹g", { size: 17 })], [R("ε→0 收敛：3.7×10", { size: 18 }), R("−3", { sup: true, size: 18 }), R("→ 1.1×10", { size: 18 }), R("−5", { sup: true, size: 18 })]],
  [[R("级联 C 独立复核", { size: 18 })], [R("单群体直接 Gillespie，10⁶ 抽样", { size: 18 })], [R("≤ ~1σ 吻合（含 m=2 解析、θ≥2）", { size: 18 })]],
  [[R("终态 t_max 收敛", { size: 18 })], [R("t_max ∈ [2×10⁴, 1.2×10⁵]", { size: 18 })], [R("χ 不变；残余激活质量 ≤ 10", { size: 18 }), R("−20", { sup: true, size: 18 })]],
  [[R("ρ(λ) 单调性", { size: 18 })], [R("8 点扫描", { size: 18 })], [R("严格递增 (0.21→2.72)，二分适定", { size: 18 })]],
  [[R("ODE 过阈行为", { size: 18 })], [R("ε=0.02，跨越 λᴄ", { size: 18 })], [R("亚临界 R(∞)~O(ε)，超临界 O(1)", { size: 18 })]],
]));

// 5. 图 2 说明 + 图
children.push(H("5  图 2 说明", HeadingLevel.HEADING_1));
const img = fs.readFileSync("figure2_lambda_c.png");
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 60, after: 60 },
  children: [new ImageRun({ type: "png", data: img, transformation: { width: 648, height: 209 } })] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 160 },
  children: [new TextRun({ text: "图 2  由亚临界终态外推测定爆发阈值 " + lam + "c，并与 " + rho + "(N)=1 独立比对。", font: FONT, italics: true, size: 18, color: GREY })] }));

const panel = (letter, runs) => children.push(P([R("（" + letter + "）", { bold: true, color: ACCENT }), ...runs], { after: 100 }));
panel("a", [R("方法。逆平均簇规模 "), R(eps), R("/R("), R(inf), R(") 随归一化控制参数 "), R(lam), R("/"), ...sub(lam, "c"), R(" 线性趋零（示例三组构型：双层 m=(3,4)、成对 6-正则、单层 m=3）。散点为 "), R(eps), R(to), R("0 外推后的闭合终态，实线为近阈值窗口 [0.90, 0.995]·"), ...sub(lam, "c"), R(" 的线性拟合；三线一同汇聚于空心圆标记处，即 "), R(rho), R("(N)=1 给出的 "), ...sub(lam, "c"), R("（归一化后为 "), R(lam), R("/"), ...sub(lam, "c"), R("=1）。")]);
panel("b", [R("一致性。八组构型的外推值与 "), R(rho), R("(N)=1 值落在对角线（identity）上，最大相对偏差 6"), R(times), R("10"), R("−4", { sup: true }), R("。误差棒为回归标准误 "), R(sig), R("（在此标度下不可见，反映一致程度之高）。")]);
panel("c", [R("偏离的 "), R(sig), R(" 数。纵轴 ("), ...sub(lam, "c"), R("^extr − "), ...sub(lam, "c"), R("^"), R(rho), R(")/"), R(sig), R(" 直接给出每组构型偏离 "), R(rho), R("(N)=1 的 "), R(sig), R(" 数；八组全部落在 1.5"), R(sig), R(" 之内，阴影为 ±1"), R(sig), R("、±2"), R(sig), R(" 带。全部为正是有限窗口凸性所致的微小系统性，量级 "), R(le), R(" 0.065%。")]);

children.push(rule());
children.push(P([R("结论：", { bold: true, color: HEAD }), R("由亚临界终态外推得到的 "), ...sub(lam, "c"), R(" 与 "), R(rho), R("(N)=1 在八组构型上一致到 6"), R(times), R("10"), R("−4", { sup: true }), R("、全部 "), R(le), R(" 1.5"), R(sig), R("；这是对爆发阈值的一次全非线性、独立的确认，与数值 Jacobian 复核（检验线性化谱）互补。")], { after: 0 }));

// ---------- build ----------
const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: 20, color: INK } } } },
  sections: [{
    properties: { page: { margin: { top: 1100, bottom: 1100, left: 1200, right: 1200 } } },
    children,
  }],
});
Packer.toBuffer(doc).then(buf => { fs.writeFileSync("lambda_c_report.docx", buf); console.log("wrote lambda_c_report.docx", buf.length, "bytes"); });
