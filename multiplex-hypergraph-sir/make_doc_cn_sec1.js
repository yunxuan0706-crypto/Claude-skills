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
k.push(P([R("本模型有三个要素，都画在了图 1 里：多重超图和层度 (a)、单个群体内部的级联 (b)，还有群体层面的分支 (c)。")], { ind: false }));
// figure 1: mechanism schematic
{
  const img = fs.readFileSync("figure1_mechanism.png");
  k.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 160, after: 70 }, indent: { firstLine: 0 },
    children: [new ImageRun({ type: "png", data: img, transformation: { width: 592, height: 197 } })] }));
  k.push(new Paragraph({ alignment: AlignmentType.BOTH, spacing: { after: 180, line: 264 }, indent: { firstLine: 0 },
    children: [
      R("图 1", { bold: true, size: SMALL }),
      R("　多重超图 SIR 群体闭合的机制示意。(a) 多重超图和节点的层度；(b) 单个群体内部的级联，", { size: SMALL }),
      ...E("(i,s)", SMALL), R(" 是递推所用的计数，它的规模 ", { size: SMALL }), ...E("C", SMALL),
      R(" 超过朴素的 ", { size: SMALL }), ...E("(m−1)T", SMALL),
      R("；(c) 群体层面的分支，来路群体（层 ", { size: SMALL }), ...E("a", SMALL),
      R("）因余度相减 ", { size: SMALL }), ...E("−δ_{ab}", SMALL), R(" 出列，所以 ", { size: SMALL }),
      ...E("N_{ab}=C(X_{ab}−δ_{ab})", SMALL), R("，这个减法只在同层 ", { size: SMALL }), ...E("a=b", SMALL),
      R(" 时才生效。", { size: SMALL }),
    ] }));
}
k.push(P([R("节点集为 "), ...E("V"), R("，"), ...E("N=|V|"), R("。层 "), ...E("a=1,…,M"), R("，每一层 a 都是 V 上的一个超图 "), ...E("H_{a}=(V,E_{a})"), R("，主线取同质基数 "), ...E("|e|=m_{a}"), R("。节点的层度 "), ...E("k^{(a)}(v)=|{e∈E_{a}:v∈e}|"), R("，写成向量就是 "), ...E("k(v)"), R("。这里引入两个生成函数：")], { ind: false }));
k.push(eq("Ψ(x)=Σ_{k}P(k)∏_{c}x_{c}^{k_c}", "1.1"));
k.push(eq("ψ_{a}(x)=(1/⟨k^{(a)}⟩)Σ_{k}k_{a}P(k)∏_{c}x_{c}^{k_c−δ_ac}", "1.2"));
k.push(P([...E("Ψ"), R(" 管节点侧，"), ...E("ψ_{a}"), R(" 管的是沿层 a 的一个成员槽到达某个节点时、这个节点的余度分布，指数里的 "), ...E("−δ_{ac}"), R(" 就是余度相减：来路的那个群体不再算进去。联合层度分布 "), ...E("P(k)"), R(" 是本模型的结构输入，层和层之间的一切效应，最后都只通过 "), ...E("Ψ"), R(" 和 "), ...E("ψ_{a}"), R(" 的交叉项进入动力学。")], { ind: false }));

// ---------------- 1.2 ----------------
k.push(H2("1.2. SIR 动力学"));
k.push(P([R("节点状态 "), ...E("x_{v}∈{S,I,R}"), R(" 是跨层共享的——一个个体只要在任一条渠道里被感染，它在所有渠道里就都算感染态。记层 a 的群体 e 里已经感染的成员数为 "), ...E("n_{e}(t)=|e∩I(t)|"), R("。当 "), ...E("n_{e}≥θ_{a}"), R(" 时群体激活，向 e 里每个易感成员各以速率 "), ...E("β_{a}"), R(" 独立传递，否则就静默着；感染者以速率 "), ...E("μ"), R(" 康复。易感节点 u 的总风险率是")], { ind: false }));
k.push(eq("Λ_{u}(t)=Σ_{a}Σ_{e∈E_a: u∈e}β_{a}1[n_{e}(t)≥θ_{a}]", "1.3"));
k.push(P([R("把时间以 "), ...E("μ"), R(" 归一，"), ...E("λ_{a}=β_{a}/μ"), R("，主线取 "), ...E("θ_{a}=1"), R("。初始时刻，每个节点各自独立地以概率 "), ...E("ε"), R(" 处在感染态、其余易感，"), ...E("ε"), R(" 就是均匀的初始感染比例。和有向超图不一样，无向群体里每个成员既是传染源、也是传染靶：一个成员被感染后群体就激活了，它能去感染同群的其他成员，后者被感染又把群体的激活时长拖长，这样就形成了群内级联。本模型里新的内容，都是从这一点来的。")], { ind: false }));

// ---------------- 1.3 ----------------
k.push(H2("1.3. 空腔构造与群体计数方程组"));
k.push(P([R("固定一个检验节点 u，让它永久保持易感（也就是空腔节点）。取一个 u 所属的层 a 群体 e，定义边基变量")], { ind: false }));
k.push(new Paragraph({
  spacing: { before: 90, after: 130, line: 288 },
  tabStops: [{ type: TabStopType.CENTER, position: CENTER }, { type: TabStopType.RIGHT, position: RIGHTT }],
  children: [new TextRun({ text: "\t" }), ...E("Φ_{a}(t)=Pr[e"), R(" 到 "), ...E("t"), R(" 时刻还没向 "), ...E("u"), R(" 传递"), ...E("]"),
             new TextRun({ text: "\t" }), new TextRun({ text: "(1.4)", font: EQF, size: BODY, color: INK })],
}));
k.push(P([R("在局域树状假设下（位形模型系综在 "), ...E("N→∞"), R(" 时局域收敛到超树，回路概率是 "), ...E("O(1/N)"), R("），u 所属的各个群体、它们的上游分支互不相交，节点侧于是就闭合成了")], { ind: false }));
k.push(eq("S(t)=(1−ε)Ψ(Φ_{1}(t),…,Φ_{M}(t))", "1.5"));
k.push(P([R("拿 e 里除 u 之外那 "), ...E("m_{a}−1"), R(" 个成员的计数 "), ...E("(s,i,r)"), R("（"), ...E("s+i+r=m_{a}−1"), R("）当仓室，定义未传递子概率 "), ...E("x^{(a)}_{sir}(t)"), R("，那么 "), ...E("Φ_{a}=Σx^{(a)}_{sir}"), R("，其中激活了、却还没传递出去的那部分是")], { ind: false }));
k.push(eq("Φ^{A}_{a}(t)=Σ_{i≥θ_a}x^{(a)}_{sir}(t)", "1.6"));
k.push(P([R("一个成员通过别的群体、到 t 时刻还保持易感的概率是 "), ...E("(1−ε)ψ_{a}(Φ)"), R("，所以它的外部风险率是")], { ind: false }));
k.push(eq("h_{a}(t)=−(d/dt)ln[(1−ε)ψ_{a}(Φ)]=Σ_{c}β_{c}Φ^{A}_{c}(∂_{c}ψ_{a}(Φ))/(ψ_{a}(Φ))", "1.7"));
k.push(P([R("一个易感成员的总感染率，是外部和群内两块之和 "), ...E("h_{a}+β_{a}1[i≥θ_{a}]"), R("；而群体向 u 的传递，以速率 "), ...E("β_{a}1[i≥θ_{a}]"), R(" 把概率移出「未传递」这一类。据此可以写出")], { ind: false }));
k.push(...eqBreak(
  "ẋ^{(a)}_{sir}=(h_{a}+β_{a}1[i−1≥θ_{a}])(s+1)x^{(a)}_{s+1,i−1,r}−(h_{a}+β_{a}1[i≥θ_{a}])s x^{(a)}_{sir}",
  "+μ[(i+1)x^{(a)}_{s,i+1,r−1}−i x^{(a)}_{sir}]−β_{a}1[i≥θ_{a}]x^{(a)}_{sir}", "1.8"));
k.push(P([R("（越界的下标都记为零。）把 (1.8) 逐项求和，"), ...E("h_{a}"), R(" 项、群内项和 "), ...E("μ"), R(" 项各自守恒、两两相消，最后只剩下")], { ind: false }));
k.push(eq("Φ̇_{a}=−β_{a}Φ^{A}_{a}", "1.9"));
k.push(P([R("这和 (1.7) 自洽。初始条件 "), ...E("x^{(a)}_{s,i,0}(0)=C(m_{a}−1,i)(1−ε)^{s}ε^{i}"), R("，其余都是零，"), ...E("Φ_{a}(0)=1"), R("。闭合方程组只含各层的群体计数，一共 "), ...E("Σ_{a}C(m_{a}+1,2)"), R(" 维，和系统规模 N 无关，"), ...E("M=2"), R("、"), ...E("m=(3,4)"), R(" 时就是 16 个；节点侧的 S 由 (1.5) 代数地给出，终态规模 "), ...E("R(∞)=1−(1−ε)Ψ(Φ(∞))"), R(" 直接读出来就行，这两样都不增加维数。要是想把完整时间演化里的 "), ...E("I(t)"), R(" 和 "), ...E("R(t)"), R(" 也画出来（图 2），只要再加一个求积 "), ...E("Ṙ=μI"), R("（"), ...E("I=1−S−R"), R("）把它们分开就够了，这一步并不改变闭合本身的维数。")], { ind: false }));

// ---------------- 1.4 ----------------
k.push(H2("1.4. 爆发阈值"));
k.push(P([R("把「感染是经由哪一层传来的」当作型，那么经由层 a 感染的节点、它的层 b 群体数的期望就是 "), ...E("⟨k^{(a)}k^{(b)}⟩/⟨k^{(a)}⟩−δ_{ab}"), R("，而每个被点燃的群体又产出 "), ...E("C(m_{b},λ_{b},θ_{b})"), R(" 个新的感染成员（C 是孤立群体内部群内级联的期望，由一个小型连续时间 Markov 链的递推算出来，"), ...E("m=2"), R(" 时就退化成了成对传递概率）。所以，多型分支矩阵是")], { ind: false }));
k.push(eq("N_{ab}=C(m_{b},λ_{b},θ_{b})[⟨k^{(a)}k^{(b)}⟩/⟨k^{(a)}⟩−δ_{ab}]", "1.10"));
k.push(P([R("爆发阈值由谱半径条件 "), ...E("ρ(N)=1"), R(" 定出来。图 2 里的控制参数 "), ...E("λ"), R("，就是拿这个阈值 "), ...E("λ_{c}"), R(" 来定标的。")], { ind: false }));

// ---------------- 1.5 ----------------
k.push(H2("1.5. 含时演化：闭合与精确仿真"));
k.push(P([R("闭合方程组不光给出阈值，也能把整个传播过程复现出来。图 2 画的是 "), ...E("M=2"), R("、"), ...E("m=(3,4)"), R("、"), ...E("N=6000"), R("、"), ...E("λ=1.6λ_{c}"), R(" 下的完整时间演化，联合层度分布 "), ...E("P(k)"), R(" 在 "), ...E("(2,2)"), R(" 和 "), ...E("(3,3)"), R(" 上各占一半，"), ...E("λ_{c}=0.0993"), R("，"), ...E("ε=0.02"), R("。群体闭合 (1.8) 和精确仿真在 "), ...E("S(t)"), R("、"), ...E("I(t)"), R(" 上的最大偏差分别是 0.0019 和 0.0006，都比图里的线宽还小；感染峰高相差大约 0.3%，峰位在 "), ...E("t≈4.4"), R("。可以说，闭合把完整的含时演化都复现了下来。")], { ind: false }));
// figure 2: time evolution
{
  const img = fs.readFileSync("figure2_evolution.png");
  k.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 160, after: 70 }, indent: { firstLine: 0 },
    children: [new ImageRun({ type: "png", data: img, transformation: { width: 560, height: 196 } })] }));
  k.push(new Paragraph({ alignment: AlignmentType.BOTH, spacing: { after: 180, line: 264 }, indent: { firstLine: 0 },
    children: [
      R("图 2", { bold: true, size: SMALL }),
      R("　多重超图 SIR 的含时演化，", { size: SMALL }), ...E("N=6000", SMALL), R("、", { size: SMALL }),
      ...E("M=2", SMALL), R("、", { size: SMALL }), ...E("m=(3,4)", SMALL), R("、", { size: SMALL }),
      ...E("λ=1.6λ_{c}", SMALL), R("、初始感染比例 ", { size: SMALL }), ...E("ε=0.02", SMALL),
      R("。实线为 400 次独立 Gillespie 实现的均值，灰带为 95% 置信区间，最宽处大约和线宽相当；虚线为群体闭合 (1.8)。(a) 易感比例 ", { size: SMALL }),
      ...E("S(t)", SMALL), R("；(b) 感染比例 ", { size: SMALL }), ...E("I(t)", SMALL), R("。", { size: SMALL }),
    ] }));
}

// ---------------- 1.6 ----------------
k.push(H2("1.6. 数值校验"));
k.push(P([R("仿真用的是 (1.3) 所定义的那个连续时间 Markov 过程的精确 Gillespie 抽样，取 "), ...E("μ=1"), R("。")], { ind: false }));
k.push(lead("内部自洽。", [R("方程组自带两条不需要仿真的自检。其一是求和规则：对 (1.8) 逐态求和，应当严格给出 (1.9)，残差落在机器精度（"), ...E("10^{−20}"), R(" 量级），这说明它就是方程组的一条代数恒等式。其二是全体成员易感时的恒等式 "), ...E("x^{(a)}_{m_a−1,0,0}=[(1−ε)ψ_{a}(Φ)]^{m_a−1}"), R("，它的残差随积分步长按 "), ...E("dt^{4}"), R(" 收敛（dt 从 0.02 逐次减半，残差依次是 "), ...E("1.3×10^{−10}"), R("、"), ...E("8.0×10^{−12}"), R("、"), ...E("5.0×10^{−13}"), R("、"), ...E("3.1×10^{−14}"), R("，相邻比大约是 16），和四阶 Runge–Kutta 的阶一致，这说明该恒等式是方程组的严格推论、残差纯粹是积分误差。")]));
k.push(lead("阈值自洽。", [R("阈值公式 (1.10) 还能和方程组自己在无病定态处的线性失稳做一次独立比对：把 (1.8) 取数值 Jacobian、对最大实部特征值二分求阈值，它和 "), ...E("ρ(N)=1"), R(" 在六组构型（含反相关度分布、非对称速率和群体规模、三层，还有成对退化情形 "), ...E("λ_{c}=1/3"), R("）上的相对差都在 "), ...E("10^{−9}"), R(" 量级，这已经受限于数值失稳判据本身的分辨率了（Jacobian 用的是中心差分）。这两条路径并不共享代码；要是把群内级联 C 换成 "), ...E("(m−1)T"), R("，这个一致性立刻就没了。")]));
k.push(lead("含时复现。", [R("像图 2 那样，群体闭合和精确仿真在整个时间演化上逐点吻合，"), ...E("S(t)"), R("、"), ...E("I(t)"), R(" 的最大偏差都比线宽还小（0.0019 和 0.0006），感染峰高相差大约 0.3%、峰位都落在 "), ...E("t≈4.4"), R("。仿真取 "), ...E("N=6000"), R("、跑了 400 次独立实现（脚本 figure2_evolution.py，里面是带 Fenwick 树加速的精确 Gillespie）。剩下的那点偏差会随系统规模变大而单调趋零，是把系统放大就能抹掉的有限尺寸效应，而不是闭合层级上的系统偏差。")]));
k.push(lead("独立复算。", [R("上面这些结果，都是用和原推导不共享路径的办法复核过的："), ...E("m=2"), R("（成对）这个退化情形下，方程组可以解析地积出终态不动点，它和 RK4 积分给出的 "), ...E("R(∞)"), R(" 只差 "), ...E("3×10^{−14}"), R("；群内级联 C 由孤立群体的直接蒙特卡洛来复核，十五组参数（"), ...E("m=2,…,6"), R(" 配 "), ...E("λ=0.5,1,2"), R("）的标准化偏差都不超过 "), ...E("2.5σ"), R("，和「没有系统偏差」相容；阈值本身也由亚临界终态外推独立测了一遍，同样和 (1.10) 的解析值相容。")]));

module.exports = { children: k };

if (require.main === module) {
  const doc = new Document({
    styles: { default: { document: { run: { font: CJK, size: BODY, color: INK } } } },
    sections: [{ properties: { page: { margin: { top: 1400, bottom: 1400, left: 1440, right: 1440 } } }, children: k }],
  });
  Packer.toBuffer(doc).then(b => { fs.writeFileSync("section1_cn.docx", b); console.log("wrote section1_cn.docx", b.length); });
}
