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

// PDG-style scaling: where residual curvature leaves chi2/dof > 1, the
// delta-method error on the x-intercept is optimistic and is inflated by
// sqrt(chi2/dof) before any significance is quoted
function scaledSE(d) {
  const out = [], fac = [];
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
    let chi2 = 0;
    for (let j = 0; j < L.length; j++) {
      const r = (y[j] - (a + b * L[j])) / e[j];
      chi2 += r * r;
    }
    const f = Math.max(1, Math.sqrt(chi2 / (L.length - 2)));
    fac.push(f); out.push(d.se[i] * f);
  }
  return { se: out, fac };
}
const SC = scaledSE(d7);
const seS0 = SC.se[0], seSEnd = SC.se[SC.se.length - 1];
const facMax = Math.max(...SC.fac).toFixed(1);
const pull0s = minus(((lc0 - d7.lc_theory) / seS0).toFixed(1));
const pullEnds = ((lcEnd - d7.lc_theory) / seSEnd).toFixed(0);
const ratioSE = (lcEnd / lc0) * Math.sqrt(Math.pow(seSEnd / lcEnd, 2) + Math.pow(seS0 / lc0, 2));
const bias0 = Math.abs((lc0 / d7.lc_theory - 1) * 100).toFixed(1);
const bias1 = Math.abs((lcEnd / d7.lc_o1 - 1) * 100).toFixed(1);
const ratioMeas = lcEnd / lc0;
const ratioPred = d7.lc_o1 / d7.lc_theory;
const ratioGap = (Math.abs(ratioMeas / ratioPred - 1) * 100).toFixed(1);

const k = [];

k.push(H1("3. 阈值的结构依赖与理论的边界"));
k.push(P([R("次代矩阵 (2.7) 将网络结构压缩为两个充分统计量 (sufficient statistics)：联合层度分布的二阶矩 "), ...E("⟨k^{(a)}k^{(b)}⟩"), R("，以及各层的超边基数 "), ...E("m_{a}"), R("。这种低维约化 (dimensional reduction) 赋予理论以预测力，但同时限定了其适用范围——被约化丢弃的结构自由度构成理论的盲区。本节从三个方向系统考察阈值的结构依赖：(i) 将阈值推广至 "), ...E("(λ_{1},λ_{2})"), R(" 平面上的相图，揭示多层协同效应 (synergy)（图 4）；(ii) 在理论可见的自由度内扫描渠道配置与群体粒度效应（图 5a、b）；(iii) 利用理论不可见的自由度——层间重叠 (inter-layer overlap)——对树状闭合假设本身进行可证伪检验 (falsifiable test)（图 5c、d）。")], { ind: false }));

// ---------------- 3.1 ----------------
k.push(H2("3.1. 协同区：多层耦合的超临界相变"));
k.push(P([R("多层网络中的协同效应 (synergistic effect) 是跨层耦合研究的核心问题之一 [5,6]：两条各自亚临界的传播渠道合并后是否足以引发宏观爆发？此前各层的速率由一个标量 "), ...E("λ"), R(" 经固定权重 "), ...E("w_{a}"), R(" 统一定标。放开约束，令两层各自的速率 "), ...E("λ_{1}"), R(" 与 "), ...E("λ_{2}"), R(" 独立变化，阈值条件 "), ...E("ρ(N)=1"), R(" 升维为 "), ...E("(λ_{1},λ_{2})"), R(" 参数平面上的一条临界曲线 (critical curve)。与之对照的是两条单层条件 "), ...E("N_{aa}=1"), R("，即各层孤立时各自的渗流阈值 [1]。三条线围出协同区")], { ind: false }));
k.push(eq("𝒮={(λ_{1},λ_{2}): N_{11}<1, N_{22}<1, ρ(N)>1}", "3.1"));
k.push(P([R("楔内每一点都对应同一个局面：任一渠道单独存在都不足以引发爆发，两者并存却足够。图 4(a) 以 "), ...E("P={(2,2),(3,3)}"), R("、"), ...E("m=(3,3)"), R(" 的谱半径热力图画出该楔形，白色等值线为 "), ...E("ρ(N)"), R(" 的等高线。单层各自的阈值为 "), ...E("λ_{c}^{(a)}=0.400567"), R("，两层对称合并后降到 0.128162，低 68.0%；§2.5 用作解析锚点的 "), ...E("λ_{1}=λ_{2}=0.13"), R("（"), ...E("N_{11}=N_{22}=0.3860"), R("、"), ...E("ρ=1.0132"), R("）正落在楔内紧贴边界处，图中以菱形标出。楔形并非临界线附近的窄缝，而是占据了整个 "), ...E("[0,λ_{c}^{(a)}]^{2}"), R(" 方块中临界曲线以上的全部区域。")]));
k.push(P([R("楔形的宽窄由层间耦合决定。取 §2.4 的 "), ...E("{2,4}"), R(" 族：边缘分布固定意味着两条单层线被钉死，"), ...E("ρ_{12}"), R(" 只能移动临界曲线本身。曲线随 "), ...E("ρ_{12}"), R(" 由 −1 增到 +1 而整体内移（图 4c），对称阈值由 0.106203 降到 0.092956，与 §2.4 的两个端点逐位相同；图 4(d) 把移动量画成"), ...E("Δλ_{2,c}(λ_{1},ρ_{12})"), R(" 的完整曲面。"), R("协同是多重结构的默认态", { bold: true }), R("——这一结论与多层网络传播的一般理论 [5,6] 一致，但本文将其从成对边推广至高阶超边：只要两层非完全孤立，合并阈值严格低于任一单层，§2.1 由 Perron–Frobenius 不等式 "), ...E("ρ(N)≥max_{a}N_{aa}"), R(" 给出严格下界。图 4(b) 将 "), ...E("λ_{c}(ρ_{12})"), R(" 的精确解析曲线与位形模型上的 Gillespie 精确仿真 [16] 自助分布 (bootstrap distribution) 对照，六点一致至 1.9σ。")]));
// Fig 4 (the phase portrait) now sits in 2.4, where the rho12 scan first
// discusses panels (b)-(d); 3.1 refers back to it.

// ---------------- 3.2 ----------------
k.push(H2("3.2. 渠道配置与群体粒度"));
k.push(lead("最劣渠道配置 (worst-case channel allocation)。", [R("固定总传播预算 "), ...E("Σ_{a}w_{a}"), R("，在各层之间转移权重，问哪种分配使爆发阈值最低（传播最危险）。这是多层网络特有的优化问题——单层网络不存在可供分配的自由度。答案在四组构型上一致："), R("最劣配置总是内解 (interior solution)", { bold: true }), R("，从不是边界解（把预算全部押在最强的那条渠道上）。原因可从 (2.7) 的矩阵结构直接读出——纯配置令 "), ...E("λ_{b}=0"), R(" 的那一层整列归零，"), R("消去了 N 的全部非对角元", { bold: true }), R("，而跨层项恰是多重结构降低阈值的唯一来源。图 5(a) 以各构型自身的最优纯配置为基准归一：对称构型 "), ...E("m=(3,3)"), R("、"), ...E("k=(3,3)"), R(" 在 "), ...E("w_{1}=0.5"), R(" 处取到最劣，阈值比最优纯配置低 29.3%；非对称构型 "), ...E("m=(2,5)"), R("、"), ...E("k=(4,2)"), R(" 最劣点移到 "), ...E("w_{1}=0.32"), R("，低 13.8%。即便层间差异使最劣点几乎贴住某一端（"), ...E("m=(3,5)"), R(" 时 "), ...E("w_{1}=0.079"), R("，仅低 0.45%），它仍严格落在内部。")]));
k.push(lead("群体粒度反转 (granularity reversal)。", [R("高阶相互作用中一个基本问题是：「少数大群体」与「多数小群体」孰更有利于传播 [3,4]？答案取决于归一化方式——在问清"), R("固定什么", { italics: true }), R("之前没有唯一答案。两种自然的归一化给出定性相反的结论（图 5b 叠示二者）。固定每节点的群体参与度 (membership) "), ...E("k=4"), R("，阈值从 "), ...E("m=2"), R(" 的 0.500000 单调降至 "), ...E("m=8"), R(" 的 0.044149——大群体确实更危险，与群内级联 "), ...E("C"), R(" 关于 "), ...E("m"), R(" 的超线性增长一致。改为固定接触预算 (contact budget) "), ...E("k(m−1)=12"), R("，即每节点可及的邻居总数不变，结论"), R("反转", { bold: true }), R("：阈值从 "), ...E("m=2"), R(" 的 0.100000 升到 "), ...E("m=7"), R(" 的 0.148224，"), R("多数小群体反而更危险", { bold: true }), R("。")]));
k.push(P([R("反转的机制可从单层阈值条件的分支因子分解中读出")], { ind: false }));
k.push(eq("(k−1)·C(m,λ_{c},θ)=1", "3.2"));
k.push(P([R("群体层面的分支因子是余度与级联的乘积。增大 "), ...E("m"), R(" 确实换来 "), ...E("C"), R(" 的超线性增益——固定接触预算时 "), ...E("C(λ_{c})"), R(" 由 0.0909 一路升到 1.0000——但同时付出余度 "), ...E("k−1"), R(" 由 11 坍缩到 1 的代价，后者取胜。极端情形最能说明问题："), ...E("k=1"), R(" 时 "), ...E("k−1=0"), R("，无论群体多大，该层单独都永远无法传播，因为被感染者没有第二个群体可以传出去。"), R("余度坍缩 (excess degree collapse) 否证了「大群体更危险」的朴素预期", { bold: true }), R("：该预期只在固定参与度 "), ...E("k"), R(" 时成立，一旦在接触预算可比条件下进行公平比较便不再成立（"), ...E("k=1"), R(" 时永不传播，未绘于图 5b）。该反转提示：在实际传播控制中，干预策略的评估必须明确归一化基准 [6]。")]));

// ---------------- 3.3 ----------------
k.push(H2("3.3. 层间重叠：树状闭合的可证伪检验"));
k.push(P([R("以上扫描都在次代矩阵 (2.7) 可见的自由度上进行。任何理论的可信度不仅取决于它在已知极限下的正确性，更取决于它在何处失效——理论边界的定量刻画是Popper 意义上的可证伪性 (falsifiability) 的体现。(2.7) 的输入仅为联合层度分布 "), ...E("P(k)"), R(" 与超边基数 "), ...E("m"), R("，因此两个 "), ...E("P(k)"), R(" 与 "), ...E("m"), R(" 相同的多重超图必被赋予同一个 "), ...E("λ_{c}"), R("，无论两层的超边在节点上如何相对排列。定义层间重叠 (inter-layer overlap)")], { ind: false }));
k.push(eq("o=#{节点对{u,v}: 同处某层1群体 ∧ 同处某层2群体}/#{节点对{u,v}: 同处某层1群体}", "3.3"));
k.push(P([R("一个被重复使用的节点对在两层间闭合一条长度为 4 的环路 (4-cycle)，"), R("恰为局域树状假设所忽略的短环路效应", { bold: true }), R(" [1,12]。由此给出一条严格的可证伪推论：固定 "), ...E("P(k)"), R("、"), ...E("m"), R(" 与 "), ...E("θ"), R("，只改 "), ...E("o"), R("，则 "), ...E("λ_{c}"), R(" 必须不动。若实测阈值随 "), ...E("o"), R(" 系统性移动，移动幅度即树状闭合失效的定量标尺 (quantitative gauge of breakdown)。")]));
k.push(lead("构造。", [R("取每节点层度 "), ...E("(2,2)"), R("、"), ...E("m=(3,3)"), R("。层 1 为位形模型；层 2 把层 1 中随机选出的一部分群体"), R("逐字复制", { bold: true }), R("，其余成员槽按残余度随机连线。由于每个节点的层 1 度恒为 2，被复制群体给它的度数不会超标，残余度非负，最终每个节点的层 2 度仍恰为 2。全族因而共享同一个 "), ...E("P(k)"), R(" 与同一个 "), ...E("m"), R("，"), ...E("o"), R(" 是唯一变动的量（脚本 overlap.py，度数与群体基数逐节点、逐群体核验无误）。")]));
k.push(lead("两条参照线。", [R("理论给出与 "), ...E("o"), R(" 无关的水平线 "), ...E("3C(3,λ_{c})=1"), R("，即 "), ...E("λ_{c}=" + f5(d7.lc_theory)), R("。另一端 "), ...E("o=1"), R(" 可解析求出：层 2 成为层 1 的副本后，每个"), R("物理", { italics: true }), R("群体同时承载两层的传播、速率合为 "), ...E("2λ"), R("，而节点只剩 2 个物理群体、余度为 1，阈值条件退化为")]));
k.push(eq("1·C(3,2λ_{c})=1  ⟹  λ_{c}=" + f5(d7.lc_o1), "3.4"));
k.push(P([R("两条参照线相差 "), ...E("" + ((d7.lc_o1 / d7.lc_theory - 1) * 100).toFixed(0) + "%"), R("——若树状闭合在重叠下依然成立，测点应贴住下面那条；若重叠确实起作用，测点应从下面那条爬向上面那条。")], { ind: false }));
k.push(lead("结果。", [R("在六个重叠水平上按 §2.2 的亚临界外推测定阈值（"), ...E("N=2×10^{4}"), R("，窗口取 "), ...E("[0.60,0.84]λ_{c}"), R(" 并三次重定心，六条拟合的残差均在 " + MAXPULL + "σ 以内），所得 "), ...E("λ_{c}"), R(" 由 "), ...E("o=0"), R(" 的 " + f5(lc0) + " 单调升到 "), ...E("o=1"), R(" 的 " + f5(lcEnd) + "，"), R("整整高出 " + rise + "%", { bold: true }), R("，而理论预言的是一条水平线。上升并非匀速：小重叠下抬升平缓，逼近 "), ...E("o=1"), R(" 时急剧加速。")]));
k.push(lead("误差的标定。", [R("截距的标准误由 delta 方法给出，其前提是线性律在整个窗口内严格成立。凡残余曲率未被消尽之处，"), ...E("χ^{2}/dof"), R(" 便大于 1，该标准误偏乐观，按惯例应以 "), ...E("√(χ^{2}/dof)"), R(" 放大。六条拟合中四条的 "), ...E("χ^{2}/dof"), R(" 在 1 附近而无需放大，"), ...E("o=1"), R(" 一条为 4.7（放大 " + facMax + " 倍）——它的窗口跨越最宽的绝对 "), ...E("λ"), R(" 区间，保留的曲率也最多。以下所有显著性均取放大后的误差。")]));
k.push(lead("偏置与稳健量。", [R("测量本身的系统偏置须予说明。"), ...E("o=0"), R(" 的控制点为 " + f5(lc0) + "，比理论值低 " + bias0 + "%（" + pull0s + "σ）。为把窗口退到安全区间，本节的拟合区比 §2.3 更远离阈值，线性律的有限窗口曲率因而更强；"), ...E("o=1"), R(" 端相对 (3.4) 同样偏低 " + bias1 + "%。两端偏置"), R("同号同量级", { bold: true }), R("，故稳健的量不是绝对阈值而是"), R("比值", { bold: true }), R("：实测")]));
k.push(eq("λ_{c}(1)/λ_{c}(0)=" + ratioMeas.toFixed(3) + "±" + ratioSE.toFixed(3) + "，  解析值 " + ratioPred.toFixed(3), "3.5"));
k.push(P([R("两者相差 " + Math.abs((ratioMeas - ratioPred) / ratioSE).toFixed(1) + "σ，"), R("在噪声之内完全相容", { bold: true }), R("。独立测得的比值与独立推出的解析比值相互印证：抬升幅度确为真实，并非窗口效应。")], { ind: false }));
k.push(P([R("推论因此被干净地证伪。"), ...E("o=1"), R(" 处偏离水平线达 " + pullEnds + "σ，逐点递增也都显著（相邻两点之差最小为 4.5σ，全部同号）。"), R("层间重叠并非高阶修正，其调控力度远超层间参与相关", { bold: true }), R("：图 4 中 "), ...E("ρ_{12}"), R(" 全程只移动阈值 12.5%，"), ...E("o"), R(" 却移动了 " + rise + "%，相差近一个数量级。")], { ind: false }));
k.push(P([R("失效的方向也符合机制。重叠把传播投向已经可及的节点，"), R("同一条边被两层重复计入", { bold: true }), R("，而 (2.7) 按互不相交的分支计数，于是系统性地高估了分支能力、低估了阈值。图 5(d) 给出小 "), ...E("λ"), R(" 下的展开：理论按 "), ...E("3C(3,λ)≈6λ"), R(" 计数，"), ...E("o=1"), R(" 的实际过程只有 "), ...E("C(3,2λ)≈4λ"), R("，比值 3:2 即高估的量级。"), R("树状闭合的适用边界由此定量给出", { bold: true }), R("：它要求层间超边近乎不相交，位形模型系综自动满足（"), ...E("o=O(1/N)"), R("），而任何具有实质层间重叠的真实网络都会落在公式之外。这一发现为将消息传递方法推广至聚类网络 (clustered networks) [19] 和具有显著层间关联的真实多层系统 [5,6] 提出了具体的修正方向。")]));
k.push(figure("figure4_structure.png", 462, 348));
k.push(caption([
  R("图 5", { bold: true, size: SMALL }),
  R("　结构依赖与理论的边界。(a) 渠道配置：固定总预算 ", { size: SMALL }), ...E("w_{1}+w_{2}", SMALL),
  R(" 在两层间转移，纵轴以各构型最优纯配置归一，圆点为最劣配置——三组构型的最劣点", { size: SMALL }),
  R("均为内点", { size: SMALL, bold: true }),
  R("。(b) 群体粒度：固定每节点群体数 ", { size: SMALL }), ...E("k=4", SMALL),
  R(" 时单层阈值随 ", { size: SMALL }), ...E("m", SMALL),
  R(" 下降，固定接触预算 ", { size: SMALL }), ...E("k(m−1)=12", SMALL),
  R(" 时却随 ", { size: SMALL }), ...E("m", SMALL),
  R(" 上升，两种归一给出", { size: SMALL }), R("相反", { size: SMALL, bold: true }),
  R("结论（纵轴对数）。(c) 层间重叠 ", { size: SMALL }), ...E("o", SMALL),
  R(" 的可证伪检验：固定 ", { size: SMALL }), ...E("P(k)", SMALL), R(" 与 ", { size: SMALL }), ...E("m", SMALL),
  R("，(2.7) 预言 ", { size: SMALL }), ...E("λ_{c}", SMALL),
  R(" 与 ", { size: SMALL }), ...E("o", SMALL),
  R(" 无关（下方蓝线），实测却单调升到 ", { size: SMALL }), ...E("o=1", SMALL),
  R(" 的解析值（上方红虚线），误差棒按 ", { size: SMALL }), ...E("√(χ^{2}/dof)", SMALL),
  R(" 标定。(d) 机制：(2.7) 所计分支因子 ", { size: SMALL }), ...E("3C(3,λ)", SMALL),
  R(" 高于 ", { size: SMALL }), ...E("o=1", SMALL), R(" 实际的 ", { size: SMALL }), ...E("C(3,2λ)", SMALL),
  R("，两者穿过 1 处即各自阈值。", { size: SMALL }),
]));

module.exports = { children: k };

if (require.main === module) {
  const doc = new Document({
    styles: { default: { document: { run: { font: CJK, size: BODY, color: INK } } } },
    sections: [{ properties: { page: { margin: { top: 1400, bottom: 1400, left: 1440, right: 1440 } } }, children: k }],
  });
  Packer.toBuffer(doc).then(b => { fs.writeFileSync("section3_cn.docx", b); console.log("wrote section3_cn.docx", b.length); });
}
