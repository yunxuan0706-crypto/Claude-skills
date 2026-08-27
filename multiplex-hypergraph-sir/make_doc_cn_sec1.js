// Section 1, in the same layout as Sections 2 and 3: numbered subsections,
// section-numbered displayed equations on centre/right tab stops, bold run-in
// labels, NO tables, figures only. Content reproduces the Part-1 manuscript
// (model, group closure, time evolution, verification); the chapters 一/二/三
// are flattened into subsections 1.1-1.6 and equations renumbered (1.x) so the
// three parts share one scheme. The mechanism schematic is Fig 1, the
// time-evolution plot Fig 2.
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, AlignmentType, TabStopType,
  BorderStyle, ImageRun,
} = require("docx");

const INK = "000000";
const CJK = "SimSun", HEI = "SimHei", EQF = "Cambria Math";
const BODY = 21, SMALL = 18;
const CENTER = 4500, RIGHTT = 9000;   // A4 text block is 9026 twips

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
// two-line displayed equation for one that does not fit a single line
const eqBreak = (s1, s2, num) => [
  new Paragraph({ spacing: { before: 90, after: 0, line: 288 },
    tabStops: [{ type: TabStopType.CENTER, position: CENTER }],
    children: [new TextRun({ text: "\t" }), ...E(s1)] }),
  new Paragraph({ spacing: { before: 0, after: 130, line: 288 },
    tabStops: [{ type: TabStopType.CENTER, position: CENTER }, { type: TabStopType.RIGHT, position: RIGHTT }],
    children: [new TextRun({ text: "\t" }), ...E(s2), new TextRun({ text: "\t" }),
               new TextRun({ text: "(" + num + ")", font: EQF, size: BODY, color: INK })] }),
];
const H1 = (t) => new Paragraph({ spacing: { before: 300, after: 150 }, indent: { firstLine: 0 },
  children: [new TextRun({ text: t, font: HEI, bold: true, size: 26, color: INK })] });
const H2 = (t) => new Paragraph({ spacing: { before: 250, after: 120 }, indent: { firstLine: 0 },
  children: [new TextRun({ text: t, font: HEI, bold: true, size: 22, color: INK })] });
const lead = (label, runs) => P([R(label, { bold: true }), ...runs]);

const k = [];

k.push(H1("1. 模型、群体闭合与含时演化"));

// ---------------- 1.1 ----------------
k.push(H2("1.1. 多重超图与记号"));
k.push(P([R("本模型的三个要素示于图 1：多重超图和层度 (a)、单个群体内部的级联 (b)，以及群体层面的分支 (c)。")], { ind: false }));
// figure 1: mechanism schematic
{
  const img = fs.readFileSync("figure1_mechanism.png");
  k.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 160, after: 70 }, indent: { firstLine: 0 },
    children: [new ImageRun({ type: "png", data: img, transformation: { width: 592, height: 197 } })] }));
  k.push(new Paragraph({ alignment: AlignmentType.BOTH, spacing: { after: 180, line: 264 }, indent: { firstLine: 0 },
    children: [
      R("图 1", { bold: true, size: SMALL }),
      R("　多重超图 SIR 群体闭合的机制示意。(a) 多重超图和节点的层度；(b) 单个群体内部的级联，", { size: SMALL }),
      ...E("(i,s)", SMALL), R(" 为递推所用的计数，其规模 ", { size: SMALL }), ...E("C", SMALL),
      R(" 超过朴素的 ", { size: SMALL }), ...E("(m−1)T", SMALL),
      R("；(c) 群体层面的分支，来路群体（层 ", { size: SMALL }), ...E("a", SMALL),
      R("）因余度相减 ", { size: SMALL }), ...E("−δ_{ab}", SMALL), R(" 出列，故 ", { size: SMALL }),
      ...E("N_{ab}=C(X_{ab}−δ_{ab})", SMALL), R("，减法仅在同层 ", { size: SMALL }), ...E("a=b", SMALL),
      R(" 生效。", { size: SMALL }),
    ] }));
}
k.push(P([R("节点集为 "), ...E("V"), R("，"), ...E("N=|V|"), R("。层 "), ...E("a=1,…,M"), R("，层 a 是 V 上的超图 "), ...E("H_{a}=(V,E_{a})"), R("，主线取同质基数 "), ...E("|e|=m_{a}"), R("。节点的层度 "), ...E("k^{(a)}(v)=|{e∈E_{a}:v∈e}|"), R("，写成向量 "), ...E("k(v)"), R("。引入两个生成函数：")], { ind: false }));
k.push(eq("Ψ(x)=Σ_{k}P(k)∏_{c}x_{c}^{k_c}", "1.1"));
k.push(eq("ψ_{a}(x)=(1/⟨k^{(a)}⟩)Σ_{k}k_{a}P(k)∏_{c}x_{c}^{k_c−δ_ac}", "1.2"));
k.push(P([...E("Ψ"), R(" 描述节点侧，"), ...E("ψ_{a}"), R(" 描述沿层 a 的一个成员槽到达某节点时该节点的余度分布，指数中的 "), ...E("−δ_{ac}"), R(" 即余度相减：来路那个群体不再计入。联合层度分布 "), ...E("P(k)"), R(" 是本模型的结构输入，层间效应只经由 "), ...E("Ψ"), R(" 和 "), ...E("ψ_{a}"), R(" 的交叉项进入动力学。")], { ind: false }));

// ---------------- 1.2 ----------------
k.push(H2("1.2. SIR 动力学"));
k.push(P([R("节点状态 "), ...E("x_{v}∈{S,I,R}"), R(" 跨层共享——个体在任一渠道被感染即在所有渠道中都为感染态。记层 a 的群体 e 中已感染成员数为 "), ...E("n_{e}(t)=|e∩I(t)|"), R("。当 "), ...E("n_{e}≥θ_{a}"), R(" 时群体激活，向 e 中每个易感成员各以速率 "), ...E("β_{a}"), R(" 独立传递，否则静默；感染者以速率 "), ...E("μ"), R(" 康复。易感节点 u 的总风险率为")], { ind: false }));
k.push(eq("Λ_{u}(t)=Σ_{a}Σ_{e∈E_a: u∈e}β_{a}1[n_{e}(t)≥θ_{a}]", "1.3"));
k.push(P([R("以 "), ...E("μ"), R(" 归一时间，"), ...E("λ_{a}=β_{a}/μ"), R("，主线取 "), ...E("θ_{a}=1"), R("。初始时刻每个节点独立地以概率 "), ...E("ε"), R(" 处于感染态、其余易感，"), ...E("ε"), R(" 即均匀初始感染比例。和有向超图不同，无向群体中每个成员既是传染源也是传染靶：一个成员被感染后群体激活，可感染同群其他成员，后者被感染又延长群体的激活时长，形成群内级联，后续分析均围绕它展开。")], { ind: false }));

// ---------------- 1.3 ----------------
k.push(H2("1.3. 空腔构造与群体计数方程组"));
k.push(P([R("固定一个检验节点 u，令其永久易感（空腔节点）。取一个 u 所属的层 a 群体 e，定义边基变量")], { ind: false }));
k.push(new Paragraph({
  spacing: { before: 90, after: 130, line: 288 },
  tabStops: [{ type: TabStopType.CENTER, position: CENTER }, { type: TabStopType.RIGHT, position: RIGHTT }],
  children: [new TextRun({ text: "\t" }), ...E("Φ_{a}(t)=Pr[e"), R(" 到 "), ...E("t"), R(" 时刻尚未向 "), ...E("u"), R(" 传递"), ...E("]"),
             new TextRun({ text: "\t" }), new TextRun({ text: "(1.4)", font: EQF, size: BODY, color: INK })],
}));
k.push(P([R("在局域树状假设下（位形模型系综 "), ...E("N→∞"), R(" 时局域收敛于超树，回路概率 "), ...E("O(1/N)"), R("），u 所属各群体的上游分支互不相交，节点侧因而闭合为")], { ind: false }));
k.push(eq("S(t)=(1−ε)Ψ(Φ_{1}(t),…,Φ_{M}(t))", "1.5"));
k.push(P([R("以 e 除 u 外的 "), ...E("m_{a}−1"), R(" 个成员的计数 "), ...E("(s,i,r)"), R("（"), ...E("s+i+r=m_{a}−1"), R("）为仓室，定义未传递子概率 "), ...E("x^{(a)}_{sir}(t)"), R("，则 "), ...E("Φ_{a}=Σx^{(a)}_{sir}"), R("，激活且未传递的部分为")], { ind: false }));
k.push(eq("Φ^{A}_{a}(t)=Σ_{i≥θ_a}x^{(a)}_{sir}(t)", "1.6"));
k.push(P([R("一个成员经由其他群体到 t 仍易感的概率为 "), ...E("(1−ε)ψ_{a}(Φ)"), R("，故其外部风险率")], { ind: false }));
k.push(eq("h_{a}(t)=−(d/dt)ln[(1−ε)ψ_{a}(Φ)]=Σ_{c}β_{c}Φ^{A}_{c}(∂_{c}ψ_{a}(Φ))/(ψ_{a}(Φ))", "1.7"));
k.push(P([R("一个易感成员的总感染率是外部和群内之和 "), ...E("h_{a}+β_{a}1[i≥θ_{a}]"), R("；群体向 u 的传递以速率 "), ...E("β_{a}1[i≥θ_{a}]"), R(" 把概率移出「未传递」类。据此")], { ind: false }));
k.push(...eqBreak(
  "ẋ^{(a)}_{sir}=(h_{a}+β_{a}1[i−1≥θ_{a}])(s+1)x^{(a)}_{s+1,i−1,r}−(h_{a}+β_{a}1[i≥θ_{a}])s x^{(a)}_{sir}",
  "+μ[(i+1)x^{(a)}_{s,i+1,r−1}−i x^{(a)}_{sir}]−β_{a}1[i≥θ_{a}]x^{(a)}_{sir}", "1.8"));
k.push(P([R("（越界下标记为零。）对 (1.8) 求和，"), ...E("h_{a}"), R(" 项、群内项和 "), ...E("μ"), R(" 项各自守恒相消，只余")], { ind: false }));
k.push(eq("Φ̇_{a}=−β_{a}Φ^{A}_{a}", "1.9"));
k.push(P([R("和 (1.7) 自洽。初始条件 "), ...E("x^{(a)}_{s,i,0}(0)=C(m_{a}−1,i)(1−ε)^{s}ε^{i}"), R("，其余为零，"), ...E("Φ_{a}(0)=1"), R("。闭合方程组只含各层的群体计数，共 "), ...E("Σ_{a}C(m_{a}+1,2)"), R(" 维，和系统规模 N 无关，"), ...E("M=2"), R("、"), ...E("m=(3,4)"), R(" 时为 16 个；节点侧的 S 由 (1.5) 代数给出，终态规模 "), ...E("R(∞)=1−(1−ε)Ψ(Φ(∞))"), R(" 直接读出，二者都不增维。欲绘出完整时间演化中的 "), ...E("I(t)"), R(" 和 "), ...E("R(t)"), R("（图 2），再附一个求积 "), ...E("Ṙ=μI"), R("（"), ...E("I=1−S−R"), R("）分离即可，它不改变闭合本身的维数。")], { ind: false }));

// ---------------- 1.4 ----------------
k.push(H2("1.4. 爆发阈值"));
k.push(P([R("按「感染经由哪一层传来」分型，经由层 a 感染的节点其层 b 群体数的期望为 "), ...E("⟨k^{(a)}k^{(b)}⟩/⟨k^{(a)}⟩−δ_{ab}"), R("，每个被点燃的群体产出 "), ...E("C(m_{b},λ_{b},θ_{b})"), R(" 个新感染成员（C 为孤立群体内的群内级联期望，由一个小型连续时间 Markov 链的递推给出，"), ...E("m=2"), R(" 时退化为成对传递概率）。故多型分支矩阵")], { ind: false }));
k.push(eq("N_{ab}=C(m_{b},λ_{b},θ_{b})[⟨k^{(a)}k^{(b)}⟩/⟨k^{(a)}⟩−δ_{ab}]", "1.10"));
k.push(P([R("爆发阈值由谱半径条件 "), ...E("ρ(N)=1"), R(" 确定。图 2(a) 的控制参数 "), ...E("λ"), R(" 即以此阈值 "), ...E("λ_{c}"), R(" 定标。")], { ind: false }));

// ---------------- 1.5 ----------------
k.push(H2("1.5. 含时演化：闭合与精确仿真"));
k.push(P([R("闭合方程组不仅给出阈值，也复现整个传播过程。图 2(a) 给出 "), ...E("M=2"), R("、"), ...E("m=(3,4)"), R("、"), ...E("N=6000"), R("、"), ...E("λ=1.6λ_{c}"), R(" 下感染比例的完整时间演化，联合层度分布 "), ...E("P(k)"), R(" 在 "), ...E("(2,2)"), R(" 和 "), ...E("(3,3)"), R(" 上各占一半，"), ...E("λ_{c}=0.0993"), R("，"), ...E("ε=0.02"), R("。群体闭合 (1.8) 和精确仿真在 "), ...E("S(t)"), R("、"), ...E("I(t)"), R(" 上的最大偏差分别为 0.0019 和 0.0006，都小于图中线宽；感染峰高相差约 0.3%，峰位 "), ...E("t≈4.4"), R("。含时演化的验证图并入图 2（见第 2 节），与阈值的两条独立测定并列。")], { ind: false }));

// ---------------- 1.6 ----------------
k.push(H2("1.6. 数值校验"));
k.push(P([R("仿真为 (1.3) 所定义连续时间 Markov 过程的精确 Gillespie 抽样，取 "), ...E("μ=1"), R("。")], { ind: false }));
k.push(lead("内部自洽。", [R("方程组自带两条不需要仿真的自检。其一为求和规则：对 (1.8) 逐态求和应严格给出 (1.9)，残差处于机器精度（"), ...E("10^{−20}"), R(" 量级），说明它是方程组的代数恒等式。其二为全体成员易感的恒等式 "), ...E("x^{(a)}_{m_a−1,0,0}=[(1−ε)ψ_{a}(Φ)]^{m_a−1}"), R("，其残差随积分步长按 "), ...E("dt^{4}"), R(" 收敛（dt 从 0.02 逐次减半，残差依次是 "), ...E("1.3×10^{−10}"), R("、"), ...E("8.0×10^{−12}"), R("、"), ...E("5.0×10^{−13}"), R("、"), ...E("3.1×10^{−14}"), R("，相邻比约 16），和四阶 Runge–Kutta 的阶一致——该恒等式是方程组的严格推论、残差纯属积分误差。")]));
k.push(lead("阈值自洽。", [R("阈值公式 (1.10) 和方程组自身在无病定态的线性失稳独立比对：把 (1.8) 作数值 Jacobian、对最大实部特征值二分求阈值，和 "), ...E("ρ(N)=1"), R(" 在六组构型（含反相关度分布、非对称速率和群体规模、三层，以及成对退化情形 "), ...E("λ_{c}=1/3"), R("）上的相对差都在 "), ...E("10^{−9}"), R(" 量级，受数值失稳判据的分辨率所限（Jacobian 取中心差分）。两条路径不共享代码；若以 "), ...E("(m−1)T"), R(" 取代群内级联 C，该一致性立即消失。")]));
k.push(lead("含时复现。", [R("如图 2(a) 所示，群体闭合和精确仿真在整个时间演化上逐点吻合，"), ...E("S(t)"), R(" 和 "), ...E("I(t)"), R(" 的最大偏差都小于线宽（0.0019 和 0.0006），感染峰高相差约 0.3%、峰位都在 "), ...E("t≈4.4"), R("。仿真为 "), ...E("N=6000"), R("、400 次独立实现（脚本 figure2_evolution.py，含 Fenwick 树加速的精确 Gillespie）。残余偏差随系统规模增大单调趋零，属有限尺寸效应，增大系统即可消除。")]));
k.push(lead("独立复算。", [R("上述结果都由和原推导不共享路径的方式复核："), ...E("m=2"), R("（成对）退化情形下方程组可解析积出终态不动点，与 RK4 积分的 "), ...E("R(∞)"), R(" 相差 "), ...E("3×10^{−14}"), R("；群内级联 C 由孤立群体的直接蒙特卡洛复核，全部十五组参数（"), ...E("m=2,…,6"), R(" 配 "), ...E("λ=0.5,1,2"), R("）的标准化偏差都不超过 "), ...E("2.5σ"), R("，和「无系统偏差」相容；阈值本身亦由亚临界终态外推独立测得，和 (1.10) 的解析值相容。")]));

module.exports = { children: k };

if (require.main === module) {
  const doc = new Document({
    styles: { default: { document: { run: { font: CJK, size: BODY, color: INK } } } },
    sections: [{ properties: { page: { margin: { top: 1400, bottom: 1400, left: 1440, right: 1440 } } }, children: k }],
  });
  Packer.toBuffer(doc).then(b => { fs.writeFileSync("section1_cn.docx", b); console.log("wrote section1_cn.docx", b.length); });
}
