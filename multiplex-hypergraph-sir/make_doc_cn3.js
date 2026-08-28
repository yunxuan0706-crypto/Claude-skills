// Section in the layout of the user's directed-hypergraph draft:
// numbered subsections, section-numbered displayed equations set on
// centre/right tab stops, bold run-in labels, NO tables, figure only.
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, AlignmentType, TabStopType,
  BorderStyle, ImageRun,
} = require("docx");

const INK = "000000";
const CJK = "SimSun", HEI = "SimHei", EQF = "Cambria Math";
const BODY = 21, SMALL = 18;          // 10.5pt body, 9pt caption
const CENTER = 4500, RIGHTT = 9000;   // A4 text block is 9026 twips   // tab stops inside the text block

// plain run
const R = (t, o = {}) => new TextRun({
  text: t, font: o.f || CJK, size: o.size || BODY, color: INK,
  bold: o.bold, italics: o.italics, subScript: o.sub, superScript: o.sup,
});
// math run(s): supports _{...} and ^{...}
function E(s, size) {
  size = size || BODY; const runs = []; let i = 0, buf = "";
  const flush = () => { if (buf) { runs.push(new TextRun({ text: buf, font: EQF, size, color: INK })); buf = ""; } };
  // strip nested _{..}/^{..} markers: Word cannot nest sub/sup, so flatten them
  const strip = (t) => t.replace(/[_^]\{([^{}]*)\}/g, "$1");
  while (i < s.length) {
    const c = s[i];
    if ((c === "_" || c === "^") && s[i + 1] === "{") {
      flush();
      let d = 1, j = i + 2;                       // match the BALANCED closing brace
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
// body paragraph; ind:true = first-line indent (new topic), false = continues an equation
const P = (runs, o = {}) => new Paragraph({
  children: Array.isArray(runs) ? runs : [runs],
  alignment: AlignmentType.BOTH,
  spacing: { after: o.after == null ? 120 : o.after, before: o.before || 0, line: 300 },
  indent: { firstLine: o.ind === false ? 0 : 420 },
});
// centred equation with right-aligned number
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
// bold run-in label + text, as in the reference
const lead = (label, runs) => P([R(label, { bold: true }), ...runs]);
const figure = (file, w, h) => new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { before: 160, after: 70 }, indent: { firstLine: 0 },
  children: [new ImageRun({ type: "png", data: fs.readFileSync(file), transformation: { width: w, height: h } })] });
const caption = (runs) => new Paragraph({
  alignment: AlignmentType.BOTH, spacing: { after: 180, line: 264 }, indent: { firstLine: 0 },
  children: runs });

const k = [];

// ============================ 2. ============================
k.push(H1("2. 由亚临界终态测定爆发阈值"));

k.push(P([R("第 1 节建立的次代矩阵阈值 "), ...E("ρ(N)=1"), R(" 源于线性化分析。线性化是否遗漏了非线性效应？本节构造一条不依赖线性化的独立判据——从亚临界区 (subcritical regime) 的终态规模发散行为反向测定阈值，与谱方法互为非平凡检验。该思路与渗流理论 (percolation theory) [1] 中从簇规模分布的临界发散读出渗流阈值 (percolation threshold) 同源。")], { ind: false }));

k.push(H2("2.1. 模型与次代矩阵"));
k.push(P([R("设节点集为 V，"), ...E("N=|V|"), R("；层 "), ...E("a=1,…,M"), R("，层 a 是 V 上的均匀超图，超边基数 "), ...E("m_{a}"), R("。节点的参与由层度向量 "), ...E("k=(k^{(1)},…,k^{(M)})"), R(" 编码，"), ...E("k^{(a)}"), R(" 为节点所属的层 a 群体数，其联合分布 "), ...E("P(k)"), R(" 是模型唯一的结构输入。群体在其感染成员数达到激活阈值 "), ...E("θ_{a}"), R(" 时被激活，随即以速率 "), ...E("β_{a}=λw_{a}μ"), R(" 向其中每个易感成员独立传递；感染者以速率 μ 康复。以 μ 归一时间后，"), ...E("λ"), R(" 是唯一的标量控制参数（"), ...E("w_{a}"), R(" 为固定的渠道权重），类似于经典 SIR 中基本再生数 "), ...E("R_{0}"), R(" [14] 的角色。")], { ind: false }));

k.push(P([R("固定一个检验节点 u 并令其永久保持易感（空腔节点），取 u 所属的一个层 a 群体 e，定义边基变量")]));
k.push(eq("Φ_{a}(t) = Pr[ e 到 t 时刻尚未向 u 传递感染 ]", "2.1"));
k.push(P([R("无向群体和有向超图 [11] 的分界正在于此。有向情形下传播图 (transmission graph) [11] 的尾集 (tail set) 和头集 (head set) 不相交，头成员被感染并不延长激活，故"), R("尾成员相互独立", { italics: true }), R("是一条可用的假设。无向群体中每个成员既是传染源也是传染靶——这与 Allard 等 [11] 指出的「有向传播中入度 (in-degree) 与出度 (out-degree) 可解耦」形成对照——一个成员被感染后群体随即激活，可感染同群其他成员，后者被感染又延长群体的激活时长，剩余成员被感染的概率随之抬高。"), R("群内级联 (intra-group cascade)", { bold: true }), R(" 使逐成员的乘积闭合失效。本文的策略不是忽略群内耦合（如 MMCA [7,10] 对逐节点状态的独立性假设），而是把耦合完整吸收进群体层面的状态空间：取 e 中除 u 外 "), ...E("m_{a}−1"), R(" 个成员的计数 "), ...E("(s,i,r)"), R("（易感、感染、康复，"), ...E("s+i+r=m_{a}−1"), R("）为房室——类似于 Ding 和 Zhu [9] 的状态转移树 (state transition tree)，但直接在连续时间下建立——定义未传递子概率")], { ind: false }));
k.push(eq("x^{(a)}_{sir}(t) = Pr[ e 尚未向 u 传递 ∧ e∖{u} 处于 (s,i,r) ]", "2.2"));
k.push(P([R("于是 "), ...E("Φ_{a}=Σ x^{(a)}_{sir}"), R("，而激活但尚未传递的部分为 "), ...E("Φ^{A}_{a}=Σ_{i≥θ_{a}} x^{(a)}_{sir}"), R("。由于 u 恒为易感，不计入群体的感染计数，群体激活与否只取决于 i。")], { ind: false }));

k.push(P([R("闭合建立在两条假设上。其一为"), R("局域树状 (locally tree-like)", { bold: true }), R("：位形模型 (configuration model) [1] 系综在 "), ...E("N→∞"), R(" 时局域收敛于超树（Bethe 格近似 [12]），任一有限邻域内出现短环路的概率为 "), ...E("O(1/N)"), R("，据此 u 所属各群体的上游分支互不相交，节点侧因而因子化为")]));
k.push(eq("S(t) = (1−ε)·Ψ(Φ_{1}(t),…,Φ_{M}(t))", "2.3"));
k.push(P([R("其中 ε 为初始均匀感染比例，Ψ 为节点侧生成函数。其二为"), R("成员外部暴露的独立同分布", { bold: true }), R("：去掉 e 后，其各成员落在互不相交的分支中，故它们经由其他群体被感染的过程相互独立。该假设弱于有向情形下的"), R("尾内无相互作用", { italics: true }), R("——它只断言"), R("外部", { italics: true }), R("暴露独立，群内耦合并未被忽略，而是已由 (2.2) 的计数状态完整承载。方程组的维数为 "), ...E("Σ_{a}C(m_{a}+1,2)"), R("，与系统规模 N 无关（相较 MMCA [7,10] 的 "), ...E("O(N)"), R(" 维和配对近似 [8] 的 "), ...E("O(N^{2})"), R(" 维，复杂度不随网络增大），积分一次即给出 "), ...E("S(t)"), R(" 与终态规模 (final size)")], { ind: false }));
k.push(eq("R(∞) = 1 − (1−ε)·Ψ(Φ(∞))", "2.4"));

k.push(lead("群内级联 (intra-group cascade)。", [R("一个种子在一个孤立的 m 元群体内造成的期望感染数并非成对独立传递的简单加和 "), ...E("(m−1)T"), R("，因为群体在整个级联期间持续激活——这一非线性反馈在成对网络中不存在，是高阶相互作用的核心区别 [3,4]。记 "), ...E("u(i,s)"), R(" 为从 i 个感染、s 个易感出发此后新增感染数的条件期望，以 "), ...E("α_{i}=λs·1[i≥θ]"), R(" 记群体的总感染速率（以康复速率为单位），由全概率公式得")]));
k.push(eq("u(i,s) = [α_{i}/(α_{i}+i)]·[1+u(i+1,s−1)] + [i/(α_{i}+i)]·u(i−1,s)", "2.5"));
k.push(P([R("边界为 "), ...E("u(0,s)=u(i,0)=0"), R("，级联规模")], { ind: false }));
k.push(eq("C(m,λ,θ) = u(1, m−1)", "2.6"));
k.push(P([R("是一个可精确计算的标量。"), ...E("m=2"), R(" 时 "), ...E("C=λ/(1+λ)=T"), R("，级联消失而回到成对传递概率；"), ...E("m≥3"), R(" 时 "), ...E("C>(m−1)T"), R("，"), ...E("λ=1"), R(" 时 "), ...E("m=3,4,5"), R(" 分别超出 11.1%、21.5% 和 31.0%，并非可忽略的修正；"), ...E("θ≥2"), R(" 时单个种子无法激活群体，故 "), ...E("C=0"), R("。图 2 给出 C 的完整定量刻画及其独立校验。")], { ind: false }));
k.push(P([R("超出量的走势本身有结构（图 2b）：超出比例并非随 "), ...E("λ"), R(" 单调，而是在 "), ...E("λ≈0.7"), R(" 附近取极大后回落。极限两端都退化——"), ...E("λ→0"), R(" 时群体来不及产生二次感染，级联无从展开；"), ...E("λ→∞"), R(" 时群体内几乎必然全员感染，朴素计数也趋于饱和上限 "), ...E("m−1"), R("。级联修正因而在中等传染强度下最为显著，恰是阈值所在的区间。")]));
k.push(figure("figure_cascade.png", 500, 273));
k.push(caption([
  R("图 2", { bold: true, size: SMALL }),
  R("　群内级联 ", { size: SMALL }), ...E("C(m,λ,θ)", SMALL),
  R(" 的定量刻画与蒙特卡洛校验。(a) ", { size: SMALL }), ...E("C", SMALL),
  R(" 随 ", { size: SMALL }), ...E("λ", SMALL),
  R(" 的增长（实线，", { size: SMALL }), ...E("m=2,…,6", SMALL),
  R("）与朴素的成对独立计数 ", { size: SMALL }), ...E("(m−1)T", SMALL),
  R("（同色虚线）之比较；灰虚线为 ", { size: SMALL }), ...E("C=1", SMALL),
  R("。(b) 超出量 ", { size: SMALL }), ...E("C/[(m−1)T]−1", SMALL),
  R("，空心圆标出 ", { size: SMALL }), ...E("λ=1", SMALL),
  R(" 处 ", { size: SMALL }), ...E("m=3,4,5", SMALL),
  R(" 的 11.1%、21.5% 与 31.0%。(c) ", { size: SMALL }), ...E("C", SMALL),
  R(" 在 ", { size: SMALL }), ...E("(λ,m)", SMALL),
  R(" 平面上的分布，色标以 ", { size: SMALL }), ...E("C=1", SMALL),
  R(" 为中心发散：冷色区单个群体不足以复制自身，暖色区已足够，实线为 ", { size: SMALL }),
  ...E("C=1", SMALL),
  R(" 等值线。", { size: SMALL }), ...E("m=2", SMALL),
  R(" 一行恒为冷色——成对规模的群体无论 ", { size: SMALL }), ...E("λ", SMALL),
  R(" 多大都不能自我复制。(d) 递推 (2.5) 的状态格点 ", { size: SMALL }),
  ...E("u(i,s)", SMALL), R("（", { size: SMALL }), ...E("m=6", SMALL), R("、", { size: SMALL }), ...E("λ=1", SMALL),
  R("），格内数字即 ", { size: SMALL }), ...E("u", SMALL),
  R(" 的取值，灰色 0 为吸收边界（", { size: SMALL }), ...E("i=0", SMALL),
  R(" 熄灭、", { size: SMALL }), ...E("s=0", SMALL),
  R(" 群体耗尽）；级联规模由方框标出的角点 ", { size: SMALL }), ...E("u(1,m−1)", SMALL),
  R(" 读出。(e) 递推值与单群体直接 Gillespie 抽样的对照（十二组 ", { size: SMALL }), ...E("(m,λ)", SMALL),
  R("，每组 ", { size: SMALL }), ...E("4×10^{5}", SMALL),
  R(" 次实现），落在对角线上。(f) 同一比较以蒙特卡洛标准误为单位，最大偏离 1.43σ，阴影为 ±1σ、±2σ 带。", { size: SMALL }),
]));

k.push(lead("次代矩阵 (next-generation matrix) [15]。", [R("分支过程必须定义在群体层面而非节点层面——这与 Allard 等 [11] 在有向传播图上以边为基本单元的计数类似。每个新感染节点点燃其余各群体，每个被点燃的群体产出 C 个新感染成员，这些成员再点燃各自其余的群体。若改在节点层面统计传染链，群内级联的后继会被重复计入。按"), R("感染经由哪一层传来", { italics: true }), R("分型——即引入传播类型 (transmission type)——经由层 a 感染的节点其层 b 群体数的期望为余度 "), ...E("⟨k^{(a)}k^{(b)}⟩/⟨k^{(a)}⟩−δ_{ab}"), R(" [1]，故")]));
k.push(eq("N_{ab} = C(m_{b},λ_{b},θ_{b})·[ ⟨k^{(a)}k^{(b)}⟩/⟨k^{(a)}⟩ − δ_{ab} ],    λ_{b}=λw_{b}", "2.7"));
k.push(P([R("爆发阈值由谱半径条件 "), ...E("ρ(N)=1"), R(" 确定。标量再生数 "), ...E("R_{0}"), R(" 由此升维为矩阵——单层网络 [1]、多层网络 [5,6]、成对相互作用与高阶相互作用 [3,4] 均为特例。极限行为一一对应：所有 "), ...E("m_{a}=2"), R(" 时 "), ...E("C=T"), R("，(2.7) 回到多层网络 SIR 的已知阈值 [5]；"), ...E("M=1"), R("、"), ...E("m=2"), R(" 时传递概率阈值 "), ...E("T_{c}=λ_{c}/(1+λ_{c})"), R(" 回到 Newman [1] 的键渗流公式 "), ...E("⟨k⟩/(⟨k^{2}⟩−⟨k⟩)"), R("。另外，N 逐元非负，Perron–Frobenius 定理给出 "), ...E("ρ(N)≥max_{a}N_{aa}"), R("，而 "), ...E("N_{aa}"), R(" 恰为孤立层 a 的分支数——这意味着多重结构的引入严格降低阈值，与多层网络的一般结论 [5,6] 一致。")], { ind: false }));

// ---------------- 2.2 ----------------
k.push(H2("2.2. 亚临界簇规模与外推判据"));
k.push(P([R("(2.7) 给出的阈值是线性失稳分析 (linear stability analysis) 的产物。本小节构造一条非线性判据：始终停在亚临界区，从单一种子引发的有限簇 (finite cluster) 的平均总规模在 "), ...E("λ→λ_{c}^{−}"), R(" 处的发散行为读出 "), ...E("λ_{c}"), R("。这与渗流理论中从平均簇规模 (mean cluster size) 的幂律发散 "), ...E("χ∼|p−p_{c}|^{−γ}"), R(" 测定渗流阈值 [1] 异曲同工。以同一多型分支过程计数，多型平均总规模为")], { ind: false }));
k.push(eq("χ(λ) = 1^{T}(I−N)^{−1} v_{0}", "2.8"));
k.push(P([R("其中 "), ...E("v_{0}"), R(" 为种子的型分布。χ "), R("不可", { bold: true }), R("写成标量 "), ...E("1/(1−ρ)"), R("：后者仅在 "), ...E("M=1"), R(" 或 "), ...E("v_{0}"), R(" 恰为 Perron 向量时成立，多层情形下 "), ...E("(I−N)^{−1}"), R(" 对各型的权重并不相同。可用的性质是 (2.8) 的极点来自 "), ...E("(I−N)^{−1}"), R("，故它恰在 "), ...E("ρ(N)=1"), R(" 处按 "), ...E("(1−ρ)^{−1}"), R(" 发散。其倒数因而是一条在阈值处线性穿零的光滑函数：")], { ind: false }));
k.push(eq("ε/R(∞) = 1/χ(λ) ≃ c·(λ_{c}−λ),    λ→λ_{c}^{−}", "2.9"));
k.push(P([R("零点即 "), ...E("λ_{c}"), R("。该判据和 "), ...E("ρ(N)=1"), R(" 不共享任何中间量：一侧把非线性方程组 (2.2)–(2.4) 积分至终态、再对 (2.9) 作线性拟合，另一侧对次代矩阵 (2.7) 求谱半径并二分，二者仅共享模型定义。两条路线的一致因而是对彼此的非平凡检验，也正是图 3 所报告的。")], { ind: false }));
k.push(lead("有限种子偏置与 Richardson 外推。", [R("(2.9) 是 "), ...E("ε→0"), R(" 的渐近陈述。有限 "), ...E("ε>0"), R(" 时即便恰在 "), ...E("λ_{c}"), R(" 处 "), ...E("R(∞)"), R(" 仍为有限，"), ...E("ε/R(∞)"), R(" 因而不归零，零点被系统性地推高——这是有限种子效应 (finite seed effect)，类似有限尺度标度 (finite-size scaling) [17] 中系统规模对临界点的偏移。由于 "), ...E("R(∞)/ε"), R(" 对小 ε 的展开以线性项为主导，该偏置可由 Richardson 外推 [18] 消去，即以 "), ...E("2f(ε/2)−f(ε)"), R(" 代替 "), ...E("f(ε)=R(∞)/ε"), R("。其效果可以定量验证：在几个有代表性的构型上，不作外推时零点一律偏高 0.37%–0.44%，作外推后降到 0.031%–0.065%。偏置的系统性（各构型同号、量级相当）连同修正的效力由此一并确认。")]));

// ---------------- 2.3 ----------------
k.push(H2("2.3. 数值结果"));
k.push(P([R("在窗口 "), ...E("λ∈[0.90,0.995]λ_{c}"), R(" 上取十个点，每点按上述方式求 "), ...E("ε→0"), R(" 的 "), ...E("ε/R(∞)"), R("，作线性拟合并取其 x 截距。八组构型的阈值覆盖 "), ...E("λ_{c}∈[0.08,0.40]"), R("，含单层高阶超图、成对退化网络和双层多重超图，其 "), ...E("ρ(N)=1"), R(" 值依次为 0.080302、0.099309、0.128162、0.174483、0.185417、0.250000、0.333333 和 0.400567，而外推值依次为 0.080354、0.099370、0.128236、0.174536、0.185512、0.250096、0.333440 和 0.400702。相对偏差依次为 +0.065%、+0.062%、+0.058%、+0.031%、+0.052%、+0.039%、+0.032% 和 +0.034%，最大者为 "), ...E("6.5×10^{−4}"), R("。以外推截距的回归标准误 σ 为单位，偏离依次为 +1.51、+1.38、+1.26、+0.55、+1.03、+0.76、+0.57 和 +0.47，最大者 1.51σ，全部落在 1.6σ 之内。")], { ind: false }));
k.push(P([R("残余偏差为线性律 (2.9) 的"), R("有限窗口效应", { bold: true }), R("，而非闭合和谱方法之间的实质分歧。"), ...E("1/χ"), R(" 仅在 "), ...E("λ_{c}"), R(" 的紧邻严格线性，有限窗口上的拟合因而带一个由曲率所致的偏置；实测该偏置的符号随窗口位置翻转——窗口取 "), ...E("[0.85,0.99]λ_{c}"), R(" 时为负、取 "), ...E("[0.92,0.998]λ_{c}"), R(" 时为正——其量级和 σ 相当，故无法和真实偏离区分。")]));

// figure 2: closure validation (time evolution + subcritical determination)
const img = fs.readFileSync("figure2_validation.png");
k.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 160, after: 70 }, indent: { firstLine: 0 },
  children: [new ImageRun({ type: "png", data: img, transformation: { width: 462, height: 348 } })] }));
k.push(new Paragraph({ alignment: AlignmentType.BOTH, spacing: { after: 180, line: 264 }, indent: { firstLine: 0 },
  children: [
    R("图 3", { bold: true, size: SMALL }),
    R("　群体闭合的两重独立验证。(a) 含时演化：感染比例 ", { size: SMALL }), ...E("I(t)", SMALL),
    R(" 的闭合解与精确 Gillespie 仿真（", { size: SMALL }), ...E("N=6000", SMALL),
    R("、400 次实现、", { size: SMALL }), ...E("λ=1.6λ_{c}", SMALL),
    R("）逐点吻合，实线为仿真、虚线为闭合、灰带为 95% 置信区间。(b) 阈值的亚临界测定：逆平均簇规模 ", { size: SMALL }),
    ...E("ε/R(∞)=1/χ", SMALL), R(" 随 ", { size: SMALL }), ...E("λ/λ_{c}", SMALL),
    R(" 线性趋零，其 x 截距即 ", { size: SMALL }), ...E("λ_{c}", SMALL),
    R("（三组示例构型，虚线为向零点的延长，空心圆为 ", { size: SMALL }), ...E("ρ(N)=1", SMALL),
    R(" 预言）。(c) 八组构型的外推 ", { size: SMALL }), ...E("λ_{c}", SMALL),
    R(" 与谱值 ", { size: SMALL }), ...E("ρ(N)=1", SMALL),
    R(" 落在对角线上，最大相对偏差 ", { size: SMALL }), ...E("6.5×10^{−4}", SMALL),
    R("，误差棒为回归标准误 σ。(d) 偏离的 σ 数，最大者 1.51σ，阴影为 ±1σ、±2σ 带。", { size: SMALL }),
  ] }));

// ---------------- 2.4 (rho12 scan) ----------------
k.push(H2("2.4. 层间度相关对阈值的效应"));
k.push(P([R("多层网络中层间度相关 (inter-layer degree correlation) 是关键的结构特征 [5,6]——一个节点在不同层的活跃程度是否相关？阈值只经由度矩 "), ...E("⟨k^{(a)}k^{(b)}⟩"), R(" 进入 N，故可构造一族"), R("边缘分布固定、仅联合分布不同", { bold: true }), R("的受控构型，隔离出纯粹的层间度相关效应。设两层度各以等概率取 2 或 4（边缘均值 3、方差 1），联合分布取 "), ...E("P(2,2)=P(4,4)=(1+ρ_{12})/4"), R("、"), ...E("P(2,4)=P(4,2)=(1−ρ_{12})/4"), R("，"), ...E("ρ_{12}"), R(" 恰为两层度的 Pearson 相关系数。"), ...E("ρ_{12}"), R(" 由 −1（完全反相关）扫到 +1（完全正相关）时边缘分布不变，唯交叉矩 "), ...E("⟨k^{(1)}k^{(2)}⟩=9+ρ_{12}"), R(" 受调控。")], { ind: false }));
k.push(P([R("在 "), ...E("m=(3,3)"), R(" 下，对角元 "), ...E("N_{aa}=C(λ)(⟨k^{2}⟩/⟨k⟩−1)=(7/3)C(λ)"), R(" 和 "), ...E("ρ_{12}"), R(" 无关，非对角元 "), ...E("N_{ab}=C(λ)⟨k^{(1)}k^{(2)}⟩/⟨k⟩"), R(" 随 "), ...E("ρ_{12}"), R(" 线性增大，故 "), ...E("ρ(N)=C(λ)(16+ρ_{12})/3"), R("，阈值由 "), ...E("C(λ_{c})=3/(16+ρ_{12})"), R(" 定出——"), R("λ_c 随 ρ_12 单调下降", { bold: true }), R("。此曲线即图 4(b) 的蓝线："), ...E("λ_{c}"), R(" 从 "), ...E("ρ_{12}=−1"), R(" 处的 0.1062 降到 "), ...E("ρ_{12}=+1"), R(" 处的 0.0930，共约 12.5%。同一族构型在位形模型多重超图上的直接 Gillespie 仿真（由亚临界簇规模外推 "), ...E("1/χ→0"), R("、"), ...E("N=2×10^{4}"), R("）在六个 "), ...E("ρ_{12}"), R(" 上复现该曲线（图 4(b) 箱线），偏离最大 1.9σ、相对差不超过 1.1%。机制明确：正相关使同一节点在两层同时充当度枢纽 (degree hub)——这与 Allard 等 [11] 所讨论的「友谊悖论 (friendship paradox) 同时放大风险与范围」一致——在边缘结构不变时集中了跨层传播能力、压低阈值；反相关则分散枢纽效应、抬高阈值。这一效应在图 4 的相图中表现为临界曲线的整体平移（图 4(c)、(d)）。")]));

// Fig 4: the phase portrait. Panels (b)-(d) are what 2.4 above discusses;
// panel (a) and the synergy wedge are taken up again in 3.1.
k.push(figure("figure3_phase.png", 452, 380));
k.push(caption([
  R("图 4", { bold: true, size: SMALL }),
  R("　爆发阈值的相图。(a) ", { size: SMALL }), ...E("P={(2,2),(3,3)}", SMALL), R("、", { size: SMALL }), ...E("m=(3,3)", SMALL),
  R(" 下 ", { size: SMALL }), ...E("ρ(N)", SMALL),
  R(" 在 ", { size: SMALL }), ...E("(λ_{1},λ_{2})", SMALL),
  R(" 平面上的热力图（magma），白色细线为等值线、加粗白线为 ", { size: SMALL }), ...E("ρ(N)=1", SMALL),
  R(" 临界曲线，两条白虚线为单层条件 ", { size: SMALL }), ...E("N_{aa}=1", SMALL),
  R("；临界曲线与两条单层线围出的斜纹填充区即 (3.1) 的协同楔 ", { size: SMALL }), ...E("𝒮", SMALL),
  R("（§3.1 详述），空心圆为解析锚点 ", { size: SMALL }), ...E("λ_{1}=λ_{2}=0.13", SMALL),
  R("（", { size: SMALL }), ...E("ρ=1.0132", SMALL),
  R("）。(b) 固定边缘的 ", { size: SMALL }), ...E("{2,4}", SMALL),
  R(" 族：精确 ", { size: SMALL }), ...E("λ_{c}(ρ_{12})", SMALL),
  R("（蓝线）与直接 Gillespie 仿真的自助箱线（", { size: SMALL }), ...E("N=2×10^{4}", SMALL),
  R("）。(c) 同族在 ", { size: SMALL }), ...E("ρ_{12}=−1,0,+1", SMALL),
  R(" 的临界曲线：边缘钉死，仅曲线随相关内移。(d) ", { size: SMALL }), ...E("Δλ_{2,c}(λ_{1},ρ_{12})", SMALL),
  R(" 相对 ", { size: SMALL }), ...E("ρ_{12}=0", SMALL),
  R(" 的移动量（发散色标）。", { size: SMALL }),
]));

// ---------------- 2.5 ----------------
k.push(H2("2.5. 数值校验"));
k.push(P([R("严格的数值校验是理论框架可靠性的基石。校验按逻辑层次分为四级：(i) 闭合方程组的内部自洽性 (internal consistency)；(ii) 已知极限下的解析锚点 (analytic benchmarks)；(iii) 亚临界簇规模的独立复算；(iv) 群内级联量的蒙特卡洛复核 (Monte Carlo verification)。四级校验所使用的代码路径互不共享，构成对理论的交叉验证。")], { ind: false }));

k.push(lead("内部自洽。", [R("对 (2.2) 的方程组求和，外部风险项、群内项和康复项各自守恒相消，只余 "), ...E("Φ̇_{a}=−β_{a}Φ^{A}_{a}"), R("；数值上该求和规则的残差为 "), ...E("10^{−17}"), R(" 量级，达机器精度。全部成员易感蕴含群体从未激活、因而必未传递，故 "), ...E("x^{(a)}_{m_{a}−1,0,0}(t)=[(1−ε)ψ_{a}(Φ(t))]^{m_{a}−1}"), R(" 对一切 t 成立；该恒等式不是额外假设而是方程组的推论，其残差为 "), ...E("10^{−12}"), R(" 量级，纯属积分误差。此外，把 (2.2) 在无病定态作数值 Jacobian，其谱横标在 "), ...E("ρ(N)=1"), R(" 所给的 "), ...E("λ_{c}"), R(" 处为零，达机器精度：闭合和次代矩阵共享同一阈值，(2.9) 的外推所测者因而确为 (2.7) 的阈值。")]));

k.push(lead("解析锚点。", [R("级联递推 (2.5)–(2.6) 在 "), ...E("m=2"), R(" 时给出 "), ...E("C=λ/(1+λ)"), R("，在 "), ...E("λ=0.13、0.35、1"), R(" 与 2 处和解析值相差不超过 "), ...E("10^{−12}"), R("；"), ...E("θ≥2"), R(" 时恒为零。阈值公式 (2.7) 在成对退化的 5-正则单层上给出 "), ...E("λ_{c}=1/3"), R("，和解析值逐位相同至十位有效数字。在 "), ...E("P={(2,2),(3,3)}"), R("、"), ...E("m=(3,3)"), R("、"), ...E("λ=0.13"), R(" 处 "), ...E("N_{11}=N_{22}=0.3860"), R("，两层各自处于深度亚临界，而合并后 "), ...E("ρ=1.0132"), R(" 已越过临界——两条各自亚临界的弱渠道合并之后可以超临界。最后，"), ...E("θ=2"), R(" 的层因 "), ...E("C=0"), R(" 而在 N 中整列为零，精确出列，所得阈值和相应单层值逐位相同：一条高门槛渠道不改变爆发与否，只改变爆发规模。")]));

k.push(lead("独立复算。", [R("闭合的亚临界终态和 (2.8) 的解析值不共享代码——一侧积分非线性方程组，另一侧求一个矩阵逆。以均匀播种写出 "), ...E("χ=1+1^{T}(I−N^{T})^{−1}g"), R("、"), ...E("g_{b}=⟨k^{(b)}⟩C(m_{b},λ_{b},θ_{b})"), R("，二者的相对差随 ε 单调下降：ε 依次取 "), ...E("2×10^{−4}"), R("、"), ...E("10^{−4}"), R("、"), ...E("5×10^{−5}"), R("、"), ...E("2×10^{−5}"), R(" 与 "), ...E("10^{−5}"), R(" 时，相对差依次为 "), ...E("3.7×10^{−3}"), R("、"), ...E("1.0×10^{−3}"), R("、"), ...E("2.7×10^{−4}"), R("、"), ...E("4.5×10^{−5}"), R(" 与 "), ...E("1.1×10^{−5}"), R("。这证实闭合的亚临界终态在 "), ...E("ε→0"), R(" 极限下即为 (2.8)，其极点严格位于 "), ...E("ρ(N)=1"), R("。终态本身亦已收敛：在 "), ...E("t_{max}∈[2×10^{4},1.2×10^{5}]"), R(" 上 χ 逐位不变，残余的激活且未传递质量在 "), ...E("10^{−18}"), R(" 以下。")]));

k.push(lead("蒙特卡洛复核。", [R("级联递推 (2.5) 由单群体的直接 Gillespie 精确抽样 [16] 独立复核，该检验既不使用方程组也不使用次代矩阵——遵循计算科学中的代码验证 (code verification) 与解的确认 (solution validation) 的标准流程。在 "), ...E("(m,λ,θ)"), R(" 的八组取值上，"), ...E("2×10^{5}"), R(" 次独立实现给出的估计和递推值一致，偏离多在一个统计标准误以内；其中 "), ...E("m=2、λ=0.35"), R(" 一组初次测得 3.4σ，以 "), ...E("10^{6}"), R(" 次实现、四个独立随机种子复测后分别为 0.64σ、0.74σ、1.12σ 与 1.03σ，确认该偏离纯属统计涨落 (statistical fluctuation) 而非系统偏差。作为整体行为的旁证，"), ...E("ρ(λ)"), R(" 在扫描区间上严格单调递增（由 0.207 增至 2.715），二分因而适定；而以 "), ...E("ε=0.02"), R(" 跨越阈值积分时，"), ...E("R(∞)"), R(" 在亚临界区保持 "), ...E("O(ε)"), R(" 量级、在超临界区抬升为 "), ...E("O(1)"), R(" 的宏观爆发 (macroscopic outbreak)，符合连续相变 (continuous phase transition) 的特征。")]));

k.push(P([R("四项校验合起来：由亚临界终态外推测得的 "), ...E("λ_{c}"), R(" 与 "), ...E("ρ(N)=1"), R(" 在八组构型上一致到 "), ...E("6.5×10^{−4}"), R("、偏离最大 1.51σ，构成对爆发阈值的全非线性独立确认——与 Jacobian 复核互补：后者检验线性化谱，前者检验非线性终态所继承的同一阈值。")]));

module.exports = { children: k };

if (require.main === module) {
  const doc = new Document({
    styles: { default: { document: { run: { font: CJK, size: BODY, color: INK } } } },
    sections: [{ properties: { page: { margin: { top: 1400, bottom: 1400, left: 1440, right: 1440 } } }, children: k }],
  });
  Packer.toBuffer(doc).then(b => { fs.writeFileSync("lambda_c_section_cn.docx", b); console.log("wrote lambda_c_section_cn.docx", b.length); });
}
