// Section 3, in the layout of Sections 1 and 2: numbered subsections,
// section-numbered displayed equations on centre/right tab stops, bold run-in
// labels, NO tables, figures only. Covers the structural dependence of the
// threshold (figures 5 and 6) and the falsifiable overlap test (figure 7).
// The figure-7 numbers are read from figure7_data.json at build time.
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, AlignmentType, TabStopType,
  BorderStyle, ImageRun,
} = require("docx");

const INK = "000000";
const CJK = "SimSun", HEI = "SimHei", EQF = "Cambria Math";
const BODY = 21, SMALL = 18;
const CENTER = 4500, RIGHTT = 9000;

const R = (t, o = {}) => new TextRun({
  text: t, font: o.f || CJK, size: o.size || BODY, color: INK,
  bold: o.bold, italics: o.italics, subScript: o.sub, superScript: o.sup,
});
function E(s, size) {
  size = size || BODY; const runs = []; let i = 0, buf = "";
  const flush = () => { if (buf) { runs.push(new TextRun({ text: buf, font: EQF, size, color: INK })); buf = ""; } };
  const strip = (t) => t.replace(/[_^]\{([^{}]*)\}/g, "$1");
  while (i < s.length) {
    const c = s[i];
    if ((c === "_" || c === "^") && s[i + 1] === "{") {
      flush();
      let d = 1, j = i + 2;
      while (j < s.length && d > 0) { if (s[j] === "{") d++; else if (s[j] === "}") d--; j++; }
      let inner = s.slice(i + 2, j - 1);
      while (/[_^]\{/.test(inner)) inner = strip(inner);
      runs.push(new TextRun({ text: inner, font: EQF, size, color: INK,
        subScript: c === "_", superScript: c === "^" }));
      i = j;
    } else { buf += c; i++; }
  }
  flush(); return runs;
}
const P = (runs, o = {}) => new Paragraph({
  children: Array.isArray(runs) ? runs : [runs],
  alignment: AlignmentType.BOTH,
  spacing: { after: o.after == null ? 120 : o.after, before: o.before || 0, line: 300 },
  indent: { firstLine: o.ind === false ? 0 : 420 },
});
const eq = (s, num) => new Paragraph({
  spacing: { before: 90, after: 130, line: 288 },
  tabStops: [{ type: TabStopType.CENTER, position: CENTER }, { type: TabStopType.RIGHT, position: RIGHTT }],
  children: [new TextRun({ text: "\t" }), ...E(s), new TextRun({ text: "\t" }),
             new TextRun({ text: "(" + num + ")", font: EQF, size: BODY, color: INK })],
});
const H1 = (t) => new Paragraph({ spacing: { before: 300, after: 150 }, indent: { firstLine: 0 },
  children: [new TextRun({ text: t, font: HEI, bold: true, size: 26, color: INK })] });
const H2 = (t) => new Paragraph({ spacing: { before: 250, after: 120 }, indent: { firstLine: 0 },
  children: [new TextRun({ text: t, font: HEI, bold: true, size: 22, color: INK })] });
const lead = (label, runs) => P([R(label, { bold: true }), ...runs]);
const figure = (file, w, h) => new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { before: 160, after: 70 }, indent: { firstLine: 0 },
  children: [new ImageRun({ type: "png", data: fs.readFileSync(file), transformation: { width: w, height: h } })] });
const caption = (runs) => new Paragraph({
  alignment: AlignmentType.BOTH, spacing: { after: 180, line: 264 }, indent: { firstLine: 0 },
  children: runs });

// ---- figure-7 numbers ----
const d7 = JSON.parse(fs.readFileSync("figure7_data.json", "utf8"));
const f4 = (x) => Number(x).toFixed(4);
const f5 = (x) => Number(x).toFixed(5);
const lc0 = d7.lc[0], lcEnd = d7.lc[d7.lc.length - 1];
const se0 = d7.se[0], seEnd = d7.se[d7.se.length - 1];
const rise = ((lcEnd / lc0 - 1) * 100).toFixed(0);
const minus = (s) => String(s).replace(/-/g, "−");   // typographic minus
const pull0 = minus(((lc0 - d7.lc_theory) / se0).toFixed(1));
const pullEnd = ((lcEnd - d7.lc_theory) / seEnd).toFixed(0);

// worst residual pull over the six weighted linear fits, so the prose cannot
// claim a fit quality the data does not have
function maxPullAll(d) {
  let worst = 0;
  for (let i = 0; i < d.o.length; i++) {
    const L = d.lams_all[i], y = d.y_all[i], e = d.ye_all[i];
    let sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (let j = 0; j < L.length; j++) {
      const w = 1 / (e[j] * e[j]);
      sw += w; sx += w * L[j]; sy += w * y[j];
      sxx += w * L[j] * L[j]; sxy += w * L[j] * y[j];
    }
    const b = (sw * sxy - sx * sy) / (sw * sxx - sx * sx);
    const a = (sy - b * sx) / sw;
    for (let j = 0; j < L.length; j++)
      worst = Math.max(worst, Math.abs((y[j] - (a + b * L[j])) / e[j]));
  }
  return worst;
}
const MAXPULL = maxPullAll(d7).toFixed(1);
const bias0 = Math.abs((lc0 / d7.lc_theory - 1) * 100).toFixed(1);
const bias1 = Math.abs((lcEnd / d7.lc_o1 - 1) * 100).toFixed(1);
const ratioMeas = lcEnd / lc0;
const ratioPred = d7.lc_o1 / d7.lc_theory;
const ratioGap = (Math.abs(ratioMeas / ratioPred - 1) * 100).toFixed(1);

const k = [];

k.push(H1("3. 阈值的结构依赖与理论的边界"));
k.push(P([R("阈值公式 (2.7) 把结构压缩成两样东西：联合层度分布的二阶矩，以及各层的群体基数。压缩得越狠，可问的问题就越具体——哪些结构自由度真的移动阈值，哪些被压掉了。本节沿三条方向扫描前者（图 5、图 6），再用一条被压掉的自由度检验公式本身（图 7）。")], { ind: false }));

// ---------------- 3.1 ----------------
k.push(H2("3.1. 协同区：两条弱渠道合并后的超临界"));
k.push(P([R("此前各层的速率由一个标量 "), ...E("λ"), R(" 经固定权重 "), ...E("w_{a}"), R(" 统一定标。放开这一约束，令两层各自的速率 "), ...E("λ_{1}"), R(" 与 "), ...E("λ_{2}"), R(" 独立变化，阈值条件 "), ...E("ρ(N)=1"), R(" 就不再是一个点，而是 "), ...E("(λ_{1},λ_{2})"), R(" 平面上的一条曲线。与之对照的是两条单层条件 "), ...E("N_{aa}=1"), R("，即各层若被孤立出来自行点燃所需的速率。三条线围出一个楔形")], { ind: false }));
k.push(eq("𝒮={(λ_{1},λ_{2}): N_{11}<1, N_{22}<1, ρ(N)>1}", "3.1"));
k.push(P([R("楔内每一点都对应同一个局面：任一渠道单独存在都不足以引发爆发，两者并存却足够。图 5(a) 以 "), ...E("P={(2,2),(3,3)}"), R("、"), ...E("m=(3,3)"), R(" 画出该楔形。单层各自的阈值为 "), ...E("λ_{c}^{(a)}=0.400567"), R("，两层对称合并后降到 0.128162，低 68.0%；§2.5 用作解析锚点的 "), ...E("λ_{1}=λ_{2}=0.13"), R("（"), ...E("N_{11}=N_{22}=0.3860"), R("、"), ...E("ρ=1.0132"), R("）正落在楔内紧贴边界处，图中以菱形标出。楔形之大值得留意：它并非临界线附近的一条窄缝，而是占据了整个 "), ...E("[0,λ_{c}^{(a)}]^{2}"), R(" 方块中临界曲线以上的全部区域。")]));
k.push(P([R("楔形的宽窄由层间耦合决定。图 5(b) 取图 4 的 "), ...E("{2,4}"), R(" 族：边缘分布固定意味着两条单层线被钉死，"), ...E("ρ_{12}"), R(" 只能移动临界曲线本身。曲线随 "), ...E("ρ_{12}"), R(" 由 −1 增到 +1 而整体下压，对称阈值由 0.106203 降到 0.092956，与图 4 的两个端点逐位相同。"), R("协同不是可有可无的修饰，而是多重结构的默认状态", { bold: true }), R("：只要两层都不是完全孤立，合并阈值就严格低于任一单层阈值——这一点在 §2.1 已由 "), ...E("ρ(N)≥max_{a}N_{aa}"), R(" 给出，图 5 把它在整个平面上兑现。")]));
k.push(figure("figure5_synergy.png", 520, 214));
k.push(caption([
  R("图 5", { bold: true, size: SMALL }),
  R("　协同区相图。(a) ", { size: SMALL }), ...E("P={(2,2),(3,3)}", SMALL), R("、", { size: SMALL }), ...E("m=(3,3)", SMALL),
  R(" 下的 ", { size: SMALL }), ...E("(λ_{1},λ_{2})", SMALL),
  R(" 平面：蓝线为 ", { size: SMALL }), ...E("ρ(N)=1", SMALL),
  R("，两条灰虚线为单层条件 ", { size: SMALL }), ...E("N_{aa}=1", SMALL),
  R("，橙色区域即 (3.1) 的楔形 ", { size: SMALL }), ...E("𝒮", SMALL),
  R("——两层各自亚临界而合并后超临界；菱形为 §2.5 的解析锚点 ", { size: SMALL }), ...E("λ_{1}=λ_{2}=0.13", SMALL),
  R("，空心圆为对称合并阈值 0.128162，细灰线为对角线。(b) 图 4 的 ", { size: SMALL }), ...E("{2,4}", SMALL),
  R(" 族：边缘固定使两条单层线不动（", { size: SMALL }), ...E("λ_{c}^{(a)}=0.2493", SMALL),
  R("），仅临界曲线随 ", { size: SMALL }), ...E("ρ_{12}", SMALL),
  R(" 移动，阴影为 ", { size: SMALL }), ...E("ρ_{12}=−1", SMALL), R(" 与 ", { size: SMALL }), ...E("+1", SMALL),
  R(" 两条曲线之间的展宽。", { size: SMALL }),
]));

// ---------------- 3.2 ----------------
k.push(H2("3.2. 渠道配置与群体粒度"));
k.push(lead("最劣配置。", [R("固定总传播预算 "), ...E("Σ_{a}w_{a}"), R("，在两层之间转移权重，问哪种分配最危险。答案在四组构型上一致："), R("最劣配置总是内点", { bold: true }), R("，从不是把预算全部押在最强的那条渠道上。理由可以从 (2.7) 直接读出——纯配置令 "), ...E("λ_{b}=0"), R(" 的那一列整列归零，"), R("把 N 的非对角元一并关掉", { bold: true }), R("，而跨层项恰是多重结构的全部收益所在。图 6(a) 以各构型自身最优纯配置为单位作归一：对称的 "), ...E("m=(3,3)"), R("、"), ...E("k=(3,3)"), R(" 在 "), ...E("w_{1}=0.5"), R(" 处取到最劣，阈值比最优纯配置低 29.3%；非对称的 "), ...E("m=(2,5)"), R("、"), ...E("k=(4,2)"), R(" 最劣点移到 "), ...E("w_{1}=0.32"), R("，低 13.8%。即便层间差异大到最劣点几乎贴住某一端（"), ...E("m=(3,5)"), R(" 时 "), ...E("w_{1}=0.079"), R("，仅低 0.45%），它仍严格落在内部。此问题为多层所特有：单层没有可供分配的自由度。")]));
k.push(lead("群体粒度。", [R("「少数大群体」与「多数小群体」孰更危险，在问清"), R("固定什么", { italics: true }), R("之前没有答案，而两种自然的固定方式给出相反的结论。固定每节点的群体数 "), ...E("k=4"), R("（图 6b），阈值从 "), ...E("m=2"), R(" 的 0.500000 一路降到 "), ...E("m=8"), R(" 的 0.044149，大群体危险得多，与 "), ...E("C"), R(" 关于 "), ...E("m"), R(" 超线性增长的预期一致。改为固定接触预算 "), ...E("k(m−1)=12"), R("，即每节点可及的邻居数相同（图 6c），结论"), R("反转", { bold: true }), R("：阈值从 "), ...E("m=2"), R(" 的 0.100000 升到 "), ...E("m=7"), R(" 的 0.148224，"), R("多数小群体反而更危险", { bold: true }), R("。")]));
k.push(P([R("反转的来源是单层阈值条件本身")], { ind: false }));
k.push(eq("(k−1)·C(m,λ_{c},θ)=1", "3.2"));
k.push(P([R("群体层面的分支因子是余度与级联的乘积。增大 "), ...E("m"), R(" 确实换来 "), ...E("C"), R(" 的超线性增益——固定接触预算时 "), ...E("C(λ_{c})"), R(" 由 0.0909 一路升到 1.0000——但同时付出余度 "), ...E("k−1"), R(" 由 11 坍缩到 1 的代价，后者取胜。极端情形最能说明问题："), ...E("k=1"), R(" 时 "), ...E("k−1=0"), R("，无论群体多大，该层单独都永远无法传播，因为被感染者没有第二个群体可以传出去。"), R("这一条否证了「大群体更危险」的朴素预期", { bold: true }), R("：该预期只在固定 "), ...E("k"), R(" 时成立，一旦按可及邻居数作公平比较便告失效。")]));
k.push(figure("figure6_design.png", 520, 158));
k.push(caption([
  R("图 6", { bold: true, size: SMALL }),
  R("　两个多层特有的结构问题。(a) 固定总预算 ", { size: SMALL }), ...E("w_{1}+w_{2}", SMALL),
  R(" 在两层间转移，纵轴以各构型自身最优纯配置归一，圆点为最劣（阈值最低）配置——四组构型的最劣点"),
  R("均为内点", { size: SMALL, bold: true }),
  R("。(b) 固定每节点群体数 ", { size: SMALL }), ...E("k=4", SMALL),
  R(" 时单层阈值随 ", { size: SMALL }), ...E("m", SMALL),
  R(" 下降（纵轴对数）。(c) 改固定接触预算 ", { size: SMALL }), ...E("k(m−1)=12", SMALL),
  R(" 后结论反转，阈值随 ", { size: SMALL }), ...E("m", SMALL),
  R(" 上升；各点旁标注相应的 ", { size: SMALL }), ...E("k", SMALL),
  R("，其中 ", { size: SMALL }), ...E("m=13、k=1", SMALL),
  R(" 因 ", { size: SMALL }), ...E("k−1=0", SMALL), R(" 而永不传播，未予绘出。", { size: SMALL }),
]));

// ---------------- 3.3 ----------------
k.push(H2("3.3. 层间重叠：一条可证伪的推论"));
k.push(P([R("以上三条扫描都在公式看得见的自由度上进行。更能说明理论边界的，是一条它"), R("看不见", { italics: true }), R("的自由度。(2.7) 只读取联合层度分布与群体基数，因此两个 "), ...E("P(k)"), R(" 与 "), ...E("m"), R(" 相同的多重超图必被赋予同一个 "), ...E("λ_{c}"), R("，无论两层的群体在节点上如何相对摆放。定义层间重叠")], { ind: false }));
k.push(eq("o=#{节点对{u,v}: 同处某层1群体 ∧ 同处某层2群体}/#{节点对{u,v}: 同处某层1群体}", "3.3"));
k.push(P([R("一个被重复使用的节点对恰好闭合一条 4-循环，"), R("正是局域树状假设所丢弃的东西", { bold: true }), R("。于是得到一条可证伪的推论：固定 "), ...E("P(k)"), R("、"), ...E("m"), R(" 与 "), ...E("θ"), R("，只改 "), ...E("o"), R("，则 "), ...E("λ_{c}"), R(" 必须不动。若测得阈值随 "), ...E("o"), R(" 移动，移动的幅度就是树状闭合失效的定量刻度。")]));
k.push(lead("构造。", [R("取每节点层度 "), ...E("(2,2)"), R("、"), ...E("m=(3,3)"), R("。层 1 为位形模型；层 2 把层 1 中随机选出的一部分群体"), R("逐字复制", { bold: true }), R("，其余成员槽按残余度随机连线。由于每个节点的层 1 度恒为 2，被复制群体给它的度数不会超标，残余度非负，最终每个节点的层 2 度仍恰为 2。全族因而共享同一个 "), ...E("P(k)"), R(" 与同一个 "), ...E("m"), R("，"), ...E("o"), R(" 是唯一变动的量（脚本 overlap.py，度数与群体基数逐节点、逐群体核验无误）。")]));
k.push(lead("两条参照线。", [R("理论给出与 "), ...E("o"), R(" 无关的水平线 "), ...E("3C(3,λ_{c})=1"), R("，即 "), ...E("λ_{c}=" + f5(d7.lc_theory)), R("。另一端 "), ...E("o=1"), R(" 可解析求出：层 2 成为层 1 的副本后，每个"), R("物理", { italics: true }), R("群体同时承载两层的传播、速率合为 "), ...E("2λ"), R("，而节点只剩 2 个物理群体、余度为 1，阈值条件退化为")]));
k.push(eq("1·C(3,2λ_{c})=1  ⟹  λ_{c}=" + f5(d7.lc_o1), "3.4"));
k.push(P([R("两条参照线相差 "), ...E("" + ((d7.lc_o1 / d7.lc_theory - 1) * 100).toFixed(0) + "%"), R("——若树状闭合在重叠下依然成立，测点应贴住下面那条；若重叠确实起作用，测点应从下面那条爬向上面那条。")], { ind: false }));
k.push(lead("结果。", [R("在六个重叠水平上按 §2.2 的亚临界外推测定阈值（"), ...E("N=2×10^{4}"), R("，窗口取 "), ...E("[0.60,0.84]λ_{c}"), R(" 并三次重定心，六条拟合的残差均在 " + MAXPULL + "σ 以内），所得 "), ...E("λ_{c}"), R(" 由 "), ...E("o=0"), R(" 的 " + f5(lc0) + " 单调升到 "), ...E("o=1"), R(" 的 " + f5(lcEnd) + "，"), R("整整高出 " + rise + "%", { bold: true }), R("，而理论预言的是一条水平线。上升并非匀速：小重叠下抬升平缓，逼近 "), ...E("o=1"), R(" 时急剧加速。")]));
k.push(lead("偏置与稳健量。", [R("测量本身的系统偏置须予说明。"), ...E("o=0"), R(" 的控制点为 " + f5(lc0) + "，比理论值低 " + bias0 + "%（" + pull0 + "σ）。为把窗口退到安全区间，本节的拟合区比 §2.3 更远离阈值，线性律的有限窗口曲率因而更强；"), ...E("o=1"), R(" 端相对 (3.4) 同样偏低 " + bias1 + "%。两端偏置"), R("同号同量级", { bold: true }), R("，故稳健的量不是绝对阈值而是"), R("比值", { bold: true }), R("：实测 "), ...E("λ_{c}(1)/λ_{c}(0)=" + ratioMeas.toFixed(3)), R("，而 (3.4) 与水平线之比的解析值为 " + ratioPred.toFixed(3) + "，两者仅差 " + ratioGap + "%。"), R("独立测得的比值与独立推出的解析比值就此相互印证", { bold: true }), R("，说明抬升的幅度本身是真实的，而非窗口选择的产物。")]));
k.push(P([R("推论因此被干净地证伪。"), ...E("o=1"), R(" 处偏离水平线达 " + pullEnd + "σ，而 "), R("层间重叠不是可以忽略的高阶修正，而是比层间参与相关更强的阈值调控量", { bold: true }), R("——图 4 中 "), ...E("ρ_{12}"), R(" 扫遍全程只移动阈值 12.5%，此处 "), ...E("o"), R(" 却移动了 " + rise + "%，相差近一个数量级。")], { ind: false }));
k.push(P([R("失效的方向也符合机制。重叠把传播投向已经可及的节点，"), R("同一条边被两层重复计入", { bold: true }), R("，而 (2.7) 按互不相交的分支计数，于是系统性地高估了分支能力、低估了阈值。图 7(c) 给出小 "), ...E("λ"), R(" 下的展开：理论按 "), ...E("3C(3,λ)≈6λ"), R(" 计数，"), ...E("o=1"), R(" 的实际过程只有 "), ...E("C(3,2λ)≈4λ"), R("，比值 3:2 即高估的量级。"), R("树状闭合的适用边界由此定量给出", { bold: true }), R("：它要求层间群体近乎不相交，位形模型系综自动满足（"), ...E("o=O(1/N)"), R("），而任何具有实质层间重叠的真实系统都会落在公式之外。")]));
k.push(figure("figure7_overlap.png", 520, 168));
k.push(caption([
  R("图 7", { bold: true, size: SMALL }),
  R("　层间重叠的可证伪检验。全族共享同一 ", { size: SMALL }), ...E("P(k)", SMALL),
  R("（每节点层度 ", { size: SMALL }), ...E("(2,2)", SMALL), R("）与同一 ", { size: SMALL }), ...E("m=(3,3)", SMALL),
  R("，仅 (3.3) 的重叠 ", { size: SMALL }), ...E("o", SMALL),
  R(" 变化。(a) 测得的 ", { size: SMALL }), ...E("λ_{c}", SMALL),
  R(" 随 ", { size: SMALL }), ...E("o", SMALL),
  R(" 单调上升，误差棒为回归标准误；下方水平线为 (2.7) 的预言（与 ", { size: SMALL }), ...E("o", SMALL),
  R(" 无关），上方水平线为 (3.4) 的 ", { size: SMALL }), ...E("o=1", SMALL),
  R(" 解析值。(b) 各 ", { size: SMALL }), ...E("o", SMALL),
  R(" 下的亚临界外推 ", { size: SMALL }), ...E("1/χ→0", SMALL),
  R("，空心圈为 x 截距，按 ", { size: SMALL }), ...E("o", SMALL),
  R(" 着色。(c) 机制：(2.7) 所计的群体层面分支因子 ", { size: SMALL }), ...E("3C(3,λ)", SMALL),
  R(" 与 ", { size: SMALL }), ...E("o=1", SMALL), R(" 实际过程的 ", { size: SMALL }), ...E("C(3,2λ)", SMALL),
  R("，两者穿过 1 的位置即各自的阈值；前者系统性偏高，重叠越大高估越甚。", { size: SMALL }),
]));

module.exports = { children: k };

if (require.main === module) {
  const doc = new Document({
    styles: { default: { document: { run: { font: CJK, size: BODY, color: INK } } } },
    sections: [{ properties: { page: { margin: { top: 1400, bottom: 1400, left: 1440, right: 1440 } } }, children: k }],
  });
  Packer.toBuffer(doc).then(b => { fs.writeFileSync("section3_cn.docx", b); console.log("wrote section3_cn.docx", b.length); });
}
