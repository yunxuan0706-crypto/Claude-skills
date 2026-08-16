// Section 3, in the same layout as Section 2: numbered subsections,
// section-numbered displayed equations on centre/right tab stops, bold run-in
// labels, NO tables, figure only.
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
const H1 = (t) => new Paragraph({ spacing: { before: 300, after: 150 }, indent: { firstLine: 0 },
  children: [new TextRun({ text: t, font: HEI, bold: true, size: 26, color: INK })] });
const H2 = (t) => new Paragraph({ spacing: { before: 250, after: 120 }, indent: { firstLine: 0 },
  children: [new TextRun({ text: t, font: HEI, bold: true, size: 22, color: INK })] });
const lead = (label, runs) => P([R(label, { bold: true }), ...runs]);

const k = [];

k.push(H1("3. 群体闭合相对度均场的偏差"));

// ---------------- 3.1 ----------------
k.push(H2("3.1. 度均场闭合"));
k.push(P([R("上一节的阈值由群体层面的分支过程给出。为界定该闭合究竟改进了什么，本节把它与同一模型的"), R("度均场闭合", { bold: true }), R("逐项对照。度均场不追踪群体内部的计数，而把群体的激活概率取为系综平均 "), ...E("Θ_{a}=1−(1−φ_{a})^{m_{a}−1}"), R("，其中 "), ...E("φ_{a}=Σ_{k}k_{a}P(k)i_{k}/⟨k^{(a)}⟩"), R(" 为层 a 的成员槽被感染节点占据的比例，于是")], { ind: false }));
k.push(eq("ṡ_{k} = −s_{k}Σ_{a}k_{a}β_{a}Θ_{a},    i̇_{k} = s_{k}Σ_{a}k_{a}β_{a}Θ_{a} − μ i_{k}", "3.1"));
k.push(P([R("在无病定态附近线性化，其增长模式由矩阵")], { ind: false }));
k.push(eq("N^{MF}_{ab} = (m_{b}−1)·λ_{b}·⟨k^{(a)}k^{(b)}⟩/⟨k^{(a)}⟩", "3.2"));
k.push(P([R("的谱半径决定，阈值同样由 "), ...E("ρ(N^{MF})=1"), R(" 给出。(3.2) 与 (3.1) 是同一件事的两种写法：把 (3.1) 在无病定态作数值 Jacobian，其谱横标在 "), ...E("ρ(N^{MF})=1"), R(" 所给的阈值处为零，残差约 "), ...E("10^{−11}"), R("。")], { ind: false }));

// ---------------- 3.2 ----------------
k.push(H2("3.2. 三处偏差及其分解"));
k.push(P([R("把 (3.2) 与上一节的精确次代矩阵逐元相除，差别恰为三处，且"), R("方向并不一致", { bold: true }), R("。其一，均场以 "), ...E("λ_{b}"), R(" 取代传递概率 "), ...E("T_{b}=λ_{b}/(1+λ_{b})"), R("，即把整个感染期的传递概率当作 "), ...E("β/μ"), R("，忽略了"), R("群体一旦向某成员传递即对其用尽", { italics: true }), R("；其二，均场缺少余度相减 "), ...E("−δ_{ab}"), R("，即允许感染沿来路群体折返。这两项都使均场"), R("高估", { bold: true }), R("再生数。其三，均场以 "), ...E("(m_{b}−1)T_{b}"), R(" 取代级联规模 "), ...E("C(m_{b},λ_{b},θ_{b})"), R("：线性化时各成员的贡献被独立相加，看不见群内感染对激活时长的延长，这一项使均场"), R("低估", { bold: true }), R("再生数。")], { ind: false }));
k.push(P([R("三者可以写成三个独立的开关。以 "), ...E("X_{ab}=⟨k^{(a)}k^{(b)}⟩/⟨k^{(a)}⟩"), R(" 记度比，令 "), ...E("s_{T},s_{D},s_{C}∈{0,1}"), R("，则")]));
k.push(eq("N(s_{T},s_{D},s_{C})_{ab} = (m_{b}−1)λ_{b}·f_{T}^{s_{T}}·f_{C}^{s_{C}}·[X_{ab} − s_{D}δ_{ab}]", "3.3"));
k.push(P([R("其中两个标量因子为")], { ind: false }));
k.push(eq("f_{T} = T_{b}/λ_{b} = 1/(1+λ_{b}) < 1,    f_{C} = C(m_{b},λ_{b},θ_{b})/[(m_{b}−1)T_{b}] > 1", "3.4"));
k.push(P([R("于是 "), ...E("(0,0,0)"), R(" 恰为均场 (3.2)，"), ...E("(1,1,1)"), R(" 恰为精确矩阵，二者在数值上分别与各自的独立实现相符至 "), ...E("10^{−9}"), R("。余度相减不是逐层的标量，其作用须在谱半径上定义：记 "), ...E("f_{D}=ρ(X−δ)/ρ(X)<1"), R("，则三项对谱半径的净效应恰为三个因子之积")], { ind: false }));
k.push(eq("ρ/ρ^{MF} = f_{T}·f_{D}·f_{C}", "3.5"));
k.push(P([R("该乘积关系已逐点验证（相对差 "), ...E("<10^{−10}"), R("）。(3.5) 小于 1 时均场高估 ρ，大于 1 时低估。")], { ind: false }));

// ---------------- 3.3 ----------------
k.push(H2("3.3. 数值结果"));
k.push(P([R("以下取 "), ...E("P={(3,3),(5,5)}"), R("（层度 "), ...E("k∈{3,5}"), R(" 且两层完全相关）、"), ...E("m_{1}=m_{2}=m"), R("、"), ...E("θ=1"), R("。图 3(a) 给出三个因子随 λ 的竞争："), ...E("f_{T}"), R(" 单调下降，"), ...E("f_{D}"), R(" 为常数 0.882，"), ...E("f_{C}"), R(" 先升后降；其乘积即净效应。逐项开关给出的阈值最能说明各项的分量：在 "), ...E("m=3"), R(" 处，相对均场只开 T 使阈值升为 1.063 倍、只开 D 升为 1.133 倍、只开 C 反而降为 0.976 倍，三项全开为 1.177 倍。"), R("余度相减是最大的单项，而级联与前两项反向。", { bold: true })], { ind: false }));
k.push(lead("阈值后果。", [R("三项净效应作用于阈值，得 "), ...E("λ_{c}^{MF}/λ_{c}"), R(" 在 "), ...E("m=2,3,4,6,8"), R(" 处依次为 0.765、0.849、0.879、0.903 与 0.913，并随 m 单调升至 "), ...E("m=16"), R(" 处的 0.927（图 3c）。"), R("均场始终低估阈值，差距随 m 缩小而不变号。", { bold: true }), R(" 在 "), ...E("m=2"), R(" 处该比值退化为一条解析恒等式：级联消失（"), ...E("f_{C}=1"), R("），而两侧各自以自身的传递概率解出阈值，"), ...E("f_{T}"), R(" 因而相消，只余")]));
k.push(eq("λ_{c}^{MF}/λ_{c} = (X−1)/X", "3.6"));
k.push(P([R("其中 X 为诸 "), ...E("X_{ab}"), R(" 的公共值。3-正则、本构型与 10-正则分别给出 2/3、13/17 与 9/10，数值验证一致至 "), ...E("10^{−9}"), R("；两层独立时 "), ...E("X_{aa}≠X_{ab}"), R("，(3.6) 相应地不成立。")], { ind: false }));
k.push(lead("净偏差变号。", [R("阈值处的图像并不能推广到整个 "), ...E("(m,λ)"), R(" 平面。图 3(b) 给出净偏差在该平面上的符号与量级：绝大部分区域为暖色，即均场高估 ρ；但在"), R("大群体、中等 λ", { bold: true }), R(" 的一片区域内 (3.5) 越过 1，级联项压过前两项，"), R("均场转为低估 ρ", { bold: true }), R("。对本构型，变号自 "), ...E("m=8"), R(" 起出现，窗口为 "), ...E("λ∈[0.086,0.394]"), R("，并随 m 加宽，"), ...E("m=16"), R(" 处 "), ...E("ρ/ρ^{MF}"), R(" 最大达 1.58。变号的"), R("起点依赖于度分布", { italics: true }), R("——它由 "), ...E("f_{D}=(X−1)/X"), R(" 定标，度越异质则越早变号：在六组分布上起点跨 "), ...E("m=6"), R(" 至 9（10-正则与重尾 "), ...E("k∈{2,20}"), R(" 为 6，3-正则与 "), ...E("k∈{1,3}"), R(" 为 9）——但变号的"), R("存在性并不依赖", { italics: true }), R("：六组分布无一例外。")]));
k.push(lead("变号区恒为超临界。", [R("变号区与阈值之间存在一条结构性的分隔：变号窗口的下沿始终位于 "), ...E("λ_{c}"), R(" 之上。该比值随 m 单调下降——"), ...E("m=8"), R(" 处约 4.7 倍、"), ...E("m=16"), R(" 处约 2.6 倍、至 "), ...E("m=50"), R(" 仍有约 2.1 倍——但在测试范围内从未接近 1，图 3b 中 "), ...E("λ_{c}"), R(" 虚线因而始终位于变号区之外。这解释了上一段两个陈述何以并存："), R("再生数的偏差会变号，阈值的偏差不会", { bold: true }), R("——阈值处 λ 太小，级联因子尚未长大。二者是两个不同的命题，不可互推。")]));

// figure
const img = fs.readFileSync("figure3_meanfield.png");
k.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 160, after: 70 }, indent: { firstLine: 0 },
  children: [new ImageRun({ type: "png", data: img, transformation: { width: 468, height: 352 } })] }));
k.push(new Paragraph({ alignment: AlignmentType.BOTH, spacing: { after: 180, line: 264 }, indent: { firstLine: 0 },
  children: [
    R("图 3", { bold: true, size: SMALL }),
    R("　逐项关闭度均场的三处偏差。(a) 三个因子随 λ 的竞争（", { size: SMALL }), ...E("m=10", SMALL),
    R("）：", { size: SMALL }), ...E("f_{C}", SMALL), R(" 为级联（使均场低估），", { size: SMALL }),
    ...E("f_{D}", SMALL), R(" 为余度相减、", { size: SMALL }), ...E("f_{T}", SMALL),
    R(" 为传递概率（二者使均场高估），黑线为三者之积；阴影标出乘积超过 1 的区间。(b) 净偏差在 ", { size: SMALL }),
    ...E("(m,λ)", SMALL), R(" 平面上的符号与量级，色标为 ", { size: SMALL }), ...E("log_{2}(ρ^{MF}/ρ)", SMALL),
    R("，黄带（黑实线）为变号边界，冷色为均场低估、暖色为高估；黑虚线为阈值 ", { size: SMALL }), ...E("λ_{c}", SMALL),
    R("，始终落在变号区之外。(c) 阈值后果 ", { size: SMALL }), ...E("λ_{c}^{MF}/λ_{c}", SMALL),
    R(" 随 m 的变化，恒小于 1。(d) 精确仿真校验：点为多重超图上的 Gillespie 抽样（", { size: SMALL }),
    ...E("N=4×10^{4}", SMALL), R("），与精确理论一致至 ", { size: SMALL }), ...E("1.3σ", SMALL),
    R("，而均场被排除 3.5σ 至 42.7σ。", { size: SMALL }),
  ] }));

// ---------------- 3.4 ----------------
k.push(H2("3.4. 数值校验"));
k.push(P([R("本节的结论依赖于两件事：均场侧的实现无误，以及精确侧确实是真值。二者分别校验如下。脚本见 verify_mf.py 与 verify_sim.py。")], { ind: false }));
k.push(lead("均场侧的内部一致。", [R("(3.2) 的矩阵与 (3.1) 的 ODE 是两条独立实现，其阈值一致：在 "), ...E("ρ(N^{MF})=1"), R(" 处 ODE 的无病态 Jacobian 谱横标为零，残差约 "), ...E("10^{−11}"), R("（"), ...E("m=2,3,4,6"), R("）。此处曾有一处数值陷阱：以前向差分构造该 Jacobian 会引入随 m 增大的 "), ...E("O(d)"), R(" 截断误差（"), ...E("Θ"), R(" 的曲率含 "), ...E("(m−1)(m−2)"), R("），在 "), ...E("m=6"), R(" 处即达 "), ...E("10^{−7}"), R("，恰与被检验的量同量级；改用中心差分后降至 "), ...E("10^{−11}"), R("。此外还独立积分了 (3.1)：在 "), ...E("0.85λ_{c}^{MF}"), R(" 与 "), ...E("1.15λ_{c}^{MF}"), R(" 两侧分别给出亚临界与超临界行为，闭合了线性化与动力学之间的回路。")]));
k.push(lead("精确侧的真值。", [R("图 3 本身是两个闭合之间的比较，故其结论须有闭合以外的真值支撑。为此在位形模型多重超图上作精确的连续时间 Gillespie 仿真（含真实回路与有限 N，不用任何闭合）。取 "), ...E("N=4×10^{4}"), R("、3000 次播种（分摊于 3 张独立实现），在 "), ...E("λ=0.6λ_{c}"), R(" 与 "), ...E("0.8λ_{c}"), R("、"), ...E("m=3,5,8,12"), R(" 共八个点上，仿真的平均簇规模与精确理论一致至 "), ...E("1.3σ"), R("，而均场的预言被排除 3.5σ 至 42.7σ（图 3d）。这同时检验了 T、余度相减与级联三项。（"), ...E("N=1.2×10^{4}"), R(" 时其中一点偏离 3.1σ，增大 N 后降至 0.4σ，确认为有限尺寸效应。）")]));
k.push(lead("变号区内的级联。", [R("由于变号区恒为超临界，它无法由亚临界簇规模直接测量。作为替代，驱动变号的级联量 "), ...E("C(m,λ)"), R(" 在"), R("变号区自身的参数上", { bold: true }), R("由单群体直接 Gillespie 抽样复核："), ...E("(m,λ)=(8,0.15)、(8,0.30)、(10,0.20)、(12,0.30)、(16,0.15)"), R("（这些点的 (3.5) 均大于 1），"), ...E("2×10^{5}"), R(" 次实现下与递推一致至 "), ...E("1.1σ"), R("。变号所依赖的每一个成分因而都已独立验证；变号本身则是该已验证矩阵在超临界区的预言。")]));

const doc = new Document({
  styles: { default: { document: { run: { font: CJK, size: BODY, color: INK } } } },
  sections: [{ properties: { page: { margin: { top: 1400, bottom: 1400, left: 1440, right: 1440 } } }, children: k }],
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync("meanfield_section_cn.docx", b); console.log("wrote meanfield_section_cn.docx", b.length); });
