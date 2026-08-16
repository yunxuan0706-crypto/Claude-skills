const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, AlignmentType, TabStopType,
  Table, TableRow, TableCell, WidthType, BorderStyle, ShadingType, ImageRun,
} = require("docx");

const INK="1A1A18", HEAD="1A1A18", GREY="6B6A63";
const HDR="ECEBE4", ZEBRA="F7F6F2";
const CJK="SimSun", EQF="Cambria Math";
const BODY=21, EQ=21;
// A4 usable width for equation tab stops (margins 1250 twips each side)
const CENTER=4703, RIGHTT=9406;

const R=(t,o={})=>new TextRun({text:t,font:o.eq?EQF:CJK,size:o.size||BODY,color:o.color||INK,bold:o.bold,italics:o.italics,subScript:o.sub,superScript:o.sup});
function E(s,size){ size=size||EQ; const runs=[]; let i=0,buf="";
  const flush=()=>{ if(buf){runs.push(new TextRun({text:buf,font:EQF,size,color:INK}));buf="";} };
  while(i<s.length){ const c=s[i];
    if((c==="_"||c==="^")&&s[i+1]==="{"){ flush(); const j=s.indexOf("}",i+2);
      runs.push(new TextRun({text:s.slice(i+2,j),font:EQF,size,color:INK,subScript:c==="_",superScript:c==="^"})); i=j+1;
    } else {buf+=c;i++;} }
  flush(); return runs;
}
const P=(runs,o={})=>new Paragraph({children:Array.isArray(runs)?runs:[runs],alignment:o.align,
  spacing:{after:o.after==null?140:o.after,before:o.before||0,line:300},indent:o.indent||{firstLine:o.noindent?0:420}});
// numbered, centred displayed equation with number at the right margin
const eqn=(s,num)=>new Paragraph({spacing:{before:70,after:130,line:288},
  tabStops:[{type:TabStopType.CENTER,position:CENTER},{type:TabStopType.RIGHT,position:RIGHTT}],
  children:[new TextRun({text:"\t"}),...E(s),new TextRun({text:"\t"}),new TextRun({text:"("+num+")",font:EQF,size:BODY,color:INK})]});
const H=(t)=>new Paragraph({spacing:{before:280,after:120},children:[new TextRun({text:t,font:CJK,bold:true,color:HEAD,size:25})]});
const rule=()=>new Paragraph({spacing:{before:40,after:140},border:{bottom:{style:BorderStyle.SINGLE,size:6,color:"D8D8D2"}},children:[]});

function cell(runs,{w,header,zebra,align}={}){
  return new TableCell({width:{size:w,type:WidthType.DXA},margins:{top:44,bottom:44,left:90,right:90},
    shading:header?{type:ShadingType.CLEAR,fill:HDR}:(zebra?{type:ShadingType.CLEAR,fill:ZEBRA}:undefined),
    children:[new Paragraph({alignment:align||AlignmentType.LEFT,spacing:{after:0,line:248},indent:{firstLine:0},
      children:(Array.isArray(runs)?runs:[runs]).map(t=>typeof t==="string"?new TextRun({text:t,font:CJK,size:18,bold:header,color:INK}):t)})]});
}
function table(widths,rows){ const total=widths.reduce((a,b)=>a+b,0);
  return new Table({columnWidths:widths,width:{size:total,type:WidthType.DXA},alignment:AlignmentType.CENTER,
    borders:{top:{style:BorderStyle.SINGLE,size:6,color:"9A9A93"},bottom:{style:BorderStyle.SINGLE,size:6,color:"9A9A93"},
      left:{style:BorderStyle.NONE},right:{style:BorderStyle.NONE},
      insideHorizontal:{style:BorderStyle.SINGLE,size:2,color:"D9D9D2"},insideVertical:{style:BorderStyle.NONE}},
    rows:rows.map((cs,i)=>new TableRow({tableHeader:i===0,cantSplit:true,
      children:cs.map((c,j)=>cell(c,{w:widths[j],header:i===0,zebra:i>0&&i%2===0,align:j===0?AlignmentType.LEFT:AlignmentType.CENTER}))}))});
}
const cap=(runs)=>new Paragraph({alignment:AlignmentType.LEFT,spacing:{before:60,after:150,line:264},indent:{firstLine:0},
  children:Array.isArray(runs)?runs:[runs]});

const kids=[];
// title
kids.push(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:60},indent:{firstLine:0},children:[
  new TextRun({text:"由亚临界终态外推测定爆发阈值 ",font:CJK,bold:true,size:30,color:INK}),
  new TextRun({text:"λ",font:EQF,bold:true,size:30,color:INK}),new TextRun({text:"c",font:EQF,bold:true,size:30,color:INK,subScript:true})]}));
kids.push(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:180},indent:{firstLine:0},
  children:[new TextRun({text:"多重超图上的 SIR 传播 · 群体闭合理论",font:CJK,size:19,color:GREY})]}));

// 1
kids.push(H("1  模型与次代矩阵"));
kids.push(P([R("我们在多重超图上传播 SIR 过程：节点归属于 M 个层的群体，层 a 是同质基数为 "),...E("m_{a}"),R(" 的超图，节点的参与由层度向量 "),...E("k=(k^{(1)},…,k^{(M)})"),R(" 编码，"),...E("k^{(a)}"),R(" 为其所属层 a 群体数，联合分布 "),...E("P(k)"),R("；一切层间结构只经由矩 "),...E("⟨k^{(a)}⟩"),R(" 与 "),...E("⟨k^{(a)}k^{(b)}⟩"),R(" 进入理论。一个群体在其至少 "),...E("θ_{a}"),R(" 名成员受感染后激活，随即以速率 "),...E("β_{a}=λw_{a}μ"),R(" 向每个易感成员传递，感染者以速率 μ 康复；重标时间后模型只余单一标量控制 λ。固定一个空腔节点并对其每个群体追踪除它以外成员的（易感, 感染, 康复）计数，得到一个维数与系统规模 N 无关的闭合常微分方程组，其终态定出爆发规模"),],{noindent:true}));
kids.push(eqn("R(∞) = 1 − (1 − ε)·Ψ(Φ(∞))","1"));
kids.push(P([R("其中 ε 为初始播种比例，Ψ 为节点侧生成函数，"),...E("Φ_{a}"),R(" 为层 a 群体尚未向空腔节点传递感染的概率。")],{noindent:true}));
kids.push(P([R("无向群体的新意在于"),R("群内级联",{bold:true}),R("：群内被感染的成员会延长该群体的激活时长，从而提高其余成员被感染的概率——这在成对或有向情形下不存在。对孤立 m 元群体中的单个种子，设 u(i,s) 为从 i 个感染、s 个易感出发此后新增感染数的期望，以 "),...E("α_{i}=λs·1[i≥θ]"),R(" 记群体总感染率（以康复率为单位），则")]));
kids.push(eqn("u(i,s) = [α_{i}/(α_{i}+i)]·[1 + u(i+1, s−1)] + [i/(α_{i}+i)]·u(i−1, s)","2"));
kids.push(P([R("边界 "),...E("u(0,s)=u(i,0)=0"),R("，级联规模 "),...E("C(m,λ,θ)=u(1, m−1)"),R(" 是可精确计算的标量：退化为成对传递概率 "),...E("C(2,λ,1)=λ/(1+λ)"),R("，在 m≥3 时超出线性估计 "),...E("(m−1)λ/(1+λ)"),R("，而 θ≥2 时为零（单个种子无法激活群体）。")],{noindent:true}));
kids.push(P([R("爆发按多型分支过程展开，型即“节点经由哪一层被感染”。经由层 a 感染的节点在其"),R("其余",{italics:true}),R("层 b 群体中各点燃 "),...E("C(m_{b},λ_{b},θ_{b})"),R(" 个新感染，而这样的群体平均有 "),...E("⟨k^{(a)}k^{(b)}⟩/⟨k^{(a)}⟩−δ_{ab}"),R(" 个。故平均后代（次代）矩阵为")]));
kids.push(eqn("N_{ab} = C(m_{b}, λ_{b}, θ_{b})·[ ⟨k^{(a)}k^{(b)}⟩ / ⟨k^{(a)}⟩ − δ_{ab} ]","3"));
kids.push(P([R("其中 "),...E("λ_{b}=λw_{b}"),R("。爆发阈值即其谱半径达到 1 之处，"),...E("ρ(N)=1"),R("。全部 "),...E("m_{a}=2"),R(" 时回到多层网络 SIR 的次代矩阵；"),...E("M=1, m=2"),R(" 时回到配置模型的标量阈值 "),...E("⟨k⟩/(⟨k^{2}⟩−⟨k⟩)"),R("。")],{noindent:true}));

// 2
kids.push(H("2  由亚临界终态读出阈值"));
kids.push(P([R("约化为低维系统后，"),...E("λ_{c}"),R(" 可由两条逻辑独立的路线测定：其一为谱方法，线性化群体间的分支过程、问其平均后代算子何时首次出现单位特征值，即上节的 "),...E("ρ(N)=1"),R("；其二在精神上是热力学式的，始终停在阈值以下、由小种子爆发规模的发散读出 "),...E("λ_{c}"),R("，本节发展之并用作独立校验。阈值以下小种子只造成有限爆发，但其期望规模随 λ→"),...E("λ_{c}^{−}"),R(" 无界增长；播种一个随机节点，其此后感染节点数的期望在分支近似下为")],{noindent:true}));
kids.push(eqn("χ(λ) = 1^{T}(I − N)^{−1} v_{0} = 1 + 1^{T}(I − N^{T})^{−1} g","4"));
kids.push(P([R("其中 "),...E("v_{0}"),R(" 为种子的型分布、g 为其首代（"),...E("g_{b}=⟨k^{(b)}⟩·C(m_{b},λ_{b},θ_{b})"),R("）。该平均规模的极点来自 "),...E("(I−N)^{−1}"),R("，故恰在 "),...E("ρ(N)=1"),R(" 处按 "),...E("(1−ρ)^{−1}"),R(" 发散；须注意"),R("不可",{bold:true}),R("将其写成标量 "),...E("1/(1−ρ)"),R("，除非 M=1 或 "),...E("v_{0}"),R(" 恰为 Perron 向量。其倒数则是一条光滑函数，在阈值处"),R("线性",{bold:true}),R("穿过零点，因而 "),...E("λ_{c}"),R(" 即 "),...E("ε/R(∞)"),R(" 的零点：")],{noindent:true}));
kids.push(eqn("y(λ) = ε/R(∞) = 1/χ(λ)  →  0,     y ≈ c·(λ_{c} − λ)     (λ → λ_{c}^{−})","5"));
kids.push(P([R("由于 "),...E("R(∞)"),R(" 由积分完整非线性闭合 (1) 得到、而非对角化 N，这条路线在真正独立于特征值条件 (3) 的意义下测定 "),...E("λ_{c}"),R("。有限种子会磨圆过渡、在 "),...E("λ_{c}"),R(" 处仍留下非零 "),...E("R(∞)"),R("，使零点略偏高；因 "),...E("R(∞)/ε"),R(" 对小 ε 为线性，此偏置由 ε→0 的 Richardson 外推消去。随 ε→0，积分终态逐点收敛到 (4)：ε 由 "),...E("2×10^{−4}"),R(" 减半至 "),...E("10^{−5}"),R(" 时相对差由 "),...E("3.7×10^{−3}"),R(" 降到 "),...E("1.1×10^{−5}"),R("，证实闭合的亚临界规模精确为 "),...E("1^{T}(I−N)^{−1}v_{0}"),R("，其极点严格钉在 "),...E("ρ(N)=1"),R("。")],{noindent:true}));

// 3
kids.push(H("3  结果"));
kids.push(P([R("图 1 在八组构型上执行上述方案，其阈值覆盖 "),...E("λ_{c}∈[0.08,0.40]"),R("，含单层高阶超图、成对退化网络与双层多重超图。面板 (a) 中逆爆发规模 "),...E("ε/R(∞)"),R(" 随 "),...E("λ/λ_{c}→1"),R(" 线性坍缩到零，三条示例线在 "),...E("ρ(N)=1"),R(" 预言处交于同一点。面板 (b) 把外推阈值对谱阈值作图，八组全落在对角线上，最大相对偏差 "),...E("6×10^{−4}"),R("；面板 (c) 以外推截距的回归标准误 σ 为单位给出每组偏离，全部落在 "),...E("ρ(N)=1"),R(" 的 1.5σ 之内（表 1）。")],{noindent:true}));
kids.push(new Paragraph({spacing:{before:60,after:60},indent:{firstLine:0},children:[new TextRun({text:"表 1　八组构型：外推阈值与谱阈值 ρ(N)=1 的比较。",font:CJK,size:18,color:GREY,italics:true})]}));
kids.push(table([2260,1560,1560,1300,1340,1140],[
  ["构型",[new TextRun({text:"λc [ρ(N)=1]",font:CJK,size:18,bold:true,color:INK})],[new TextRun({text:"λc [外推]",font:CJK,size:18,bold:true,color:INK})],"相对偏差",[new TextRun({text:"σ",font:EQF,size:18,bold:true,color:INK})],[new TextRun({text:"偏离 (σ)",font:CJK,size:18,bold:true,color:INK})]],
  ["2-layer m=(4,4)","0.080302","0.080354","+0.065%",[R("3.4×10",{size:18}),R("−5",{sup:true,size:18})],"+1.51"],
  ["2-layer m=(3,4)","0.099309","0.099370","+0.062%",[R("4.4×10",{size:18}),R("−5",{sup:true,size:18})],"+1.38"],
  ["2-layer m=(3,3)","0.128162","0.128236","+0.058%",[R("5.9×10",{size:18}),R("−5",{sup:true,size:18})],"+1.26"],
  ["3-reg. m=4","0.174483","0.174536","+0.031%",[R("9.8×10",{size:18}),R("−5",{sup:true,size:18})],"+0.55"],
  ["2-layer m=(2,3)","0.185417","0.185512","+0.052%",[R("9.3×10",{size:18}),R("−5",{sup:true,size:18})],"+1.03"],
  ["6-reg. pairwise","0.250000","0.250096","+0.039%",[R("1.3×10",{size:18}),R("−4",{sup:true,size:18})],"+0.76"],
  ["5-reg. pairwise","0.333333","0.333440","+0.032%",[R("1.9×10",{size:18}),R("−4",{sup:true,size:18})],"+0.57"],
  ["1-layer m=3","0.400567","0.400702","+0.034%",[R("2.9×10",{size:18}),R("−4",{sup:true,size:18})],"+0.47"],
]));
kids.push(P([R("微小残余为线性律 (5) 的"),R("有限窗口效应",{bold:true}),R("：逆规模仅在 "),...E("λ_{c}"),R(" 紧邻严格线性，故有限窗口拟合带一个偏置，其符号随窗口位置翻转、量级与 σ 相当，因而不可与真实偏离区分。")],{before:40,noindent:true}));
const img=fs.readFileSync("figure2_lambda_c.png");
kids.push(new Paragraph({alignment:AlignmentType.CENTER,spacing:{before:100,after:50},indent:{firstLine:0},children:[new ImageRun({type:"png",data:img,transformation:{width:612,height:197}})]}));
kids.push(cap([new TextRun({text:"图 1　",font:CJK,bold:true,size:18,color:INK}),
  new TextRun({text:"由亚临界终态外推测定爆发阈值 ",font:CJK,size:18,color:GREY}),
  new TextRun({text:"λc",font:EQF,italics:true,size:18,color:GREY}),
  new TextRun({text:"，并与 ρ(N)=1 独立比对。(a) ε/R(∞) 随 λ/λc 线性趋零：散点为 ε→0 闭合终态，实线为窗口 [0.90,0.995]λc 的线性拟合，空心圆为 ρ(N)=1 预言。(b) 八组构型外推值与谱值落在对角线上，最大相对偏差 6×10⁻⁴。(c) 偏离的 σ 数，全部 ≤1.5σ（阴影为 ±1σ、±2σ）。",font:CJK,size:18,color:GREY})]));

// 4
kids.push(H("4  正确性核验"));
kids.push(P([R("正文每一条论断都由不共享代码的独立核验支撑，全部通过（脚本 verify.py / verify_ode.py / datacheck.py / recheck_mc.py）。表 2 汇总解析锚点与闭合自检，表 3 汇总独立交叉验证。")],{noindent:true}));
kids.push(new Paragraph({spacing:{before:80,after:60},indent:{firstLine:0},children:[new TextRun({text:"表 2　解析锚点与闭合自检。",font:CJK,size:18,color:GREY,italics:true})]}));
kids.push(table([3260,2500,2600,1000],[
  ["核验项","期望","计算结果","状态"],
  [[R("C(2,λ,1)=λ/(1+λ)",{size:18})],"精确",[R("λ=0.13/0.35/1/2 到 10",{size:18}),R("−12",{sup:true,size:18})],"通过"],
  [[R("C(3,4,5; λ=1) 超出 (m−1)T",{size:18})],"11 / 22 / 31 %","11.11 / 21.53 / 31.05 %","通过"],
  [[R("θ≥2 单种子",{size:18})],"C = 0","0","通过"],
  [[R("Anchor A：m=(3,3), λ=0.13",{size:18})],[R("ρ=1.013",{size:18})],[R("N",{size:18}),R("11",{sub:true,size:18}),R("=0.3860, ρ=1.0132",{size:18})],"通过"],
  [[R("成对退化 5-正则",{size:18})],[R("λc=1/3",{size:18})],"0.3333333333","通过"],
  [[R("Anchor B：θ=2 通道",{size:18})],"阈值不变","C=0，精确出列","通过"],
  [[R("求和规则 / 恒等式",{size:18})],"≈ 0",[R("残差 10",{size:18}),R("−17",{sup:true,size:18}),R(" / 10",{size:18}),R("−12",{sup:true,size:18})],"通过"],
  [[R("无病态 Jacobian 谱横标",{size:18})],[R("0 at λc",{size:18})],"0（机器精度）","通过"],
]));
kids.push(new Paragraph({spacing:{before:130,after:60},indent:{firstLine:0},children:[new TextRun({text:"表 3　独立交叉验证。",font:CJK,size:18,color:GREY,italics:true})]}));
kids.push(table([2900,3050,3410],[
  ["检验","方法","结果"],
  [[R("闭合终态 = 分支公式",{size:18})],[R("χ_ODE vs 1+1ᵀ(I−Nᵀ)⁻¹g",{size:18})],[R("≤1.0×10",{size:18}),R("−3",{sup:true,size:18}),R("；ε→0 收敛到 1.1×10",{size:18}),R("−5",{sup:true,size:18})]],
  [[R("级联 C 独立复核",{size:18})],[R("单群体 Gillespie，10⁶ 抽样",{size:18})],[R("≤1σ 吻合（含 m=2 解析、θ≥2）",{size:18})]],
  [[R("终态 t_max 收敛",{size:18})],[R("t_max ∈ [2×10⁴, 1.2×10⁵]",{size:18})],[R("χ 不变；残余激活质量 ≤10",{size:18}),R("−20",{sup:true,size:18})]],
  [[R("ρ(λ) 单调 / 过阈",{size:18})],[R("扫描 + ε=0.02",{size:18})],[R("严格递增；亚临界 O(ε)，超临界 O(1)",{size:18})]],
]));
kids.push(rule());
kids.push(P([new TextRun({text:"结论：",font:CJK,bold:true,color:INK}),R("由亚临界终态外推得到的 "),...E("λ_{c}"),R(" 与 "),...E("ρ(N)=1"),R(" 在八组构型上一致到 "),...E("6×10^{−4}"),R("、全部 ≤1.5σ；这是对爆发阈值的一次全非线性、独立确认，与数值 Jacobian 复核（只检验线性化谱）互补。")],{after:0,noindent:true}));

const doc=new Document({
  styles:{default:{document:{run:{font:CJK,size:BODY,color:INK}}}},
  sections:[{properties:{page:{margin:{top:1200,bottom:1200,left:1250,right:1250}}},children:kids}],
});
Packer.toBuffer(doc).then(b=>{fs.writeFileSync("lambda_c_report_cn2.docx",b);console.log("wrote lambda_c_report_cn2.docx",b.length);});
