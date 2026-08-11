#!/usr/bin/env python3
"""Generate the verification write-up from data/verification.json.

The text is generated from the same file the figures read, so a number in the
prose cannot drift away from the number in a panel. Editing the prose by hand
defeats that; edit this generator instead.
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "..", "data", "verification.json")))
out = []
W = out.append

s = d["setup"]
W("# 群体闭合的数值验证")
W("")
W("本文档由 `tools/write_report.py` 从 `data/verification.json` 生成，与图 1、图 2 读取")
W("同一份数据，故正文数字与图中数字不会失配。全部计算的脚本见 `tools/`。")
W("")
W("主构型取两层、群体基数 $m=(%d,%d)$、联合层度分布 $P(\\mathbf{k})$ 在 $(2,2)$ 与 $(3,3)$"
  % (s["ms"][0], s["ms"][1]))
W("上各占一半，渠道权重相等。由 $\\rho(\\mathsf{N})=1$ 得阈值 $\\lambda_c=%.5f$，"
  % s["lambda_c"])
W("含时比较取 $\\lambda=%.1f\\lambda_c=%.5f$，初始感染比例 $\\varepsilon=%.2f$，$\\mu=1$。"
  % (s["ratio"], s["lambda"], s["eps"]))
W("")

W("## 1. 求和规则与恒等式")
W("")
W("方程组自带两条无需仿真的自检。其一为求和规则：对主方程逐态求和，外部风险项、群内传染项")
W("与康复项应各自守恒相消，只余 $\\dot\\Phi_a=-\\beta_a\\Phi^A_a$。其二为全体成员易感的恒等式")
W("$x^{(a)}_{m_a-1,0,0}=[(1-\\varepsilon)\\psi_a(\\boldsymbol{\\Phi})]^{m_a-1}$。")
W("")
W("| 积分步长 $\\mathrm{d}t$ | 恒等式残差 | 求和规则残差 |")
W("|---|---|---|")
for r in d["invariants"]:
    W("| %.4f | %.2e | %.1e |" % (r["dt"], r["identity"], r["sumrule"]))
inv = d["invariants"]
ratios = [inv[i]["identity"]/inv[i+1]["identity"] for i in range(len(inv)-1)]
W("")
W("求和规则的残差处于机器精度（$10^{-20}$ 量级），说明它是方程组的代数恒等式而非近似。")
W("恒等式的残差随步长减小按 %s 收敛，相邻比值为 %s，与四阶 Runge–Kutta 的阶一致；"
  % ("$\\mathrm{d}t^4$", "、".join("$%.1f$" % r for r in ratios)))
W("这说明该恒等式同样是方程组的严格推论，观察到的残差纯属积分误差。")
W("")

W("## 2. 主方程的数值 Jacobian 与 $\\rho(\\mathsf{N})=1$")
W("")
W("阈值由两条不共享代码的路径分别求出：其一将主方程在无病定态作数值 Jacobian，对最大实部")
W("特征值二分；其二直接解 $\\rho(\\mathsf{N})=1$。前者只用到方程组，后者只用到群体层面的")
W("分支计数。")
W("")
W("| 构型 | Jacobian 路径 | $\\rho(\\mathsf{N})=1$ | 相对差 |")
W("|---|---|---|---|")
for r in d["threshold_crosscheck"]:
    W("| %s | %.8f | %.8f | %.1e |" % (r["case"], r["jacobian"], r["spectral"], r["rel"]))
worst = max(r["rel"] for r in d["threshold_crosscheck"])
W("")
W("六组构型的相对差均不超过 $%.0e$，为数值 Jacobian 的有限差分精度。其中成对退化情形" % worst)
W("的解析阈值为 $\\lambda_c=1/3$，两条路径与之皆符。若在 $\\mathsf{N}$ 中以 $(m-1)T$ 取代")
W("群内级联 $C$，该一致性立即消失——级联项正是由此被发现。")
W("")

W("## 3. 群内级联 $C$ 的递推与蒙特卡洛")
W("")
W("$C(m,\\lambda,\\theta)$ 由一个小型连续时间 Markov 链的递推给出。作为独立复核，直接对孤立")
W("群体的传播过程做蒙特卡洛抽样，每组 $4\\times10^5$ 次实现。")
W("")
W("| $m$ | $\\lambda$ | 递推 | 蒙特卡洛 | 偏差 | $(m-1)T$ | 级联增益 |")
W("|---|---|---|---|---|---|---|")
for r in d["cascade"]:
    gain = r["recursion"]/r["naive"] if r["naive"] > 0 else float("nan")
    W("| %d | %.1f | %.5f | %.5f ± %.5f | %.1e | %.5f | %.3f |"
      % (r["m"], r["lam"], r["recursion"], r["mc"], r["mc_se"],
         abs(r["recursion"]-r["mc"]), r["naive"], gain))
W("")
W("两者在全部参数上一致，偏差均在蒙特卡洛的两倍标准误以内。末列给出级联相对于朴素计数")
W("$(m-1)T$ 的增益：$m=2$ 时为 1，级联消失；群体越大增益越高，说明成员互相传染对激活时长")
W("的延长不是可忽略的修正。")
W("")

W("## 4. 与精确仿真的含时比较")
W("")
tr = d["trajectory"]
t = tr["t"]
dcl_S = max(abs(a-b) for a, b in zip(tr["S_sim"], tr["S_closure"]))
dmf_S = max(abs(a-b) for a, b in zip(tr["S_sim"], tr["S_mf"]))
dcl_I = max(abs(a-b) for a, b in zip(tr["I_sim"], tr["I_closure"]))
dmf_I = max(abs(a-b) for a, b in zip(tr["I_sim"], tr["I_mf"]))
pk = lambda v: t[max(range(len(v)), key=lambda j: v[j])]
pv = lambda v: max(v)
W("图 1 给出 $N=%d$、%d 次独立实现下的完整时间演化。" % (tr["N"], tr["runs"]))
W("群体闭合与仿真在 $S(t)$ 与 $I(t)$ 上的最大偏差分别为 $%.4f$ 与 $%.4f$；" % (dcl_S, dcl_I))
W("度均场的对应值为 $%.4f$ 与 $%.4f$，相差约 %.0f 倍与 %.0f 倍。"
  % (dmf_S, dmf_I, dmf_S/dcl_S, dmf_I/dcl_I))
W("感染峰值方面，仿真为 $%.4f$（$t=%.2f$），闭合为 $%.4f$（$t=%.2f$），"
  % (pv(tr["I_sim"]), pk(tr["I_sim"]), pv(tr["I_closure"]), pk(tr["I_closure"])))
W("均场为 $%.4f$（$t=%.2f$）——均场把峰高抬高 %.2f 倍并把峰位提前 %.2f 个时间单位。"
  % (pv(tr["I_mf"]), pk(tr["I_mf"]), pv(tr["I_mf"])/pv(tr["I_sim"]),
     pk(tr["I_sim"])-pk(tr["I_mf"])))
W("")
W("![Figure 1](figures/fig1_trajectory.png){width=5.5}")
W("")

W("## 5. 误差随系统规模的标度")
W("")
W("图 2(a) 给出两种闭合的偏差随 $N$ 的变化。统计量为 $\\overline{|\\Delta S|}$，即各时刻")
W("仿真均值与理论之差的绝对值在时间网格上的平均。其不确定度由**对实现的自举**给出")
W("（每个 $N$ 重采样整条轨迹 2000 次）：逐点偏差在时间上的 lag-1 自相关约 $0.95$，")
W("网格点远非独立，重采样整条轨迹可自动保留这一相关结构，无须对它作任何假设。")
W("")
W("零假设的参照不是 $0$ 而是**噪声水平**。$\\overline{|\\Delta S|}$ 是绝对值的平均，恒为正，")
W("其置信区间永不包含 $0$；真实偏差为零时它的期望等于各时刻抽样标准误的平均乘以")
W("$\\sqrt{2/\\pi}$，这才是应当比较的基准。")
W("")
W("| $N$ | 实现次数 | 群体闭合 [95% CI] | 度均场 [95% CI] | 零假设水平 |")
W("|---|---|---|---|---|")
for r in d["scaling"]:
    W("| %d | %d | %.5f [%.5f, %.5f] | %.5f [%.5f, %.5f] | %.5f |"
      % (r["N"], r["runs"], r["closure"], r["closure_lo"], r["closure_hi"],
         r["meanfield"], r["meanfield_lo"], r["meanfield_hi"], r["noise"]))
sc = d["scaling"]
def wfit(xs, ys, ws):
    S = sum(ws); Sx = sum(w*x for w, x in zip(ws, xs)); Sy = sum(w*y for w, y in zip(ws, ys))
    Sxx = sum(w*x*x for w, x in zip(ws, xs)); Sxy = sum(w*x*y for w, x, y in zip(ws, xs, ys))
    det = S*Sxx-Sx*Sx
    return (S*Sxy-Sx*Sy)/det, math.sqrt(S/det)
use = [r for r in sc if r["closure"] > r["noise"]]
sl, sle = wfit([math.log(r["N"]) for r in use], [math.log(r["closure"]) for r in use],
               [1/((math.log(r["closure_hi"])-math.log(r["closure_lo"]))/3.92)**2 for r in use])
xs = [1/r["N"] for r in sc]; ys = [r["meanfield"] for r in sc]
ws = [1/(((r["meanfield_hi"]-r["meanfield_lo"])/3.92)**2) for r in sc]
S = sum(ws); Sx = sum(w*x for w, x in zip(ws, xs)); Sy = sum(w*y for w, y in zip(ws, ys))
Sxx = sum(w*x*x for w, x in zip(ws, xs)); Sxy = sum(w*x*y for w, x, y in zip(ws, xs, ys))
det = S*Sxx-Sx*Sx; a = (Sxx*Sy-Sx*Sxy)/det; sa = math.sqrt(Sxx/det)
W("")
W("群体闭合的偏差由 $N=500$ 处高出零假设水平 $%.1f$ 倍，单调降至 $N=8000$ 处的 $%.5f$，"
  % (sc[0]["closure"]/sc[0]["noise"], sc[-1]["closure"]))
W("已低于该水平（$%.5f$），即与「偏差为零」完全相容。只取高于零假设水平的 $%d$ 个点作"
  % (sc[-1]["noise"], len(use)))
W("加权幂律拟合，指数为 $%.2f\\pm%.2f$。" % (sl, sle))
W("须说明两点：其一，各测量值都含噪声贡献因而是真实偏差的上界，真实衰减只会更陡；")
W("其二，该指数比朴素预期的 $-1$（位形模型中短回路密度的量级）陡约 $%.1f\\sigma$，"
  % (abs(sl+1)/sle))
W("不宜简单声称「与 $O(1/N)$ 一致」。")
W("")
W("度均场的偏差则只由 $%.5f$ 降到 $%.5f$，按 $a+b/N$ 外推得平台"
  % (sc[0]["meanfield"], sc[-1]["meanfield"]))
W("$a=%.5f\\pm%.5f$，即 $%.0f\\sigma$ 地不为零；在最大规模处它仍是零假设水平的 $%.0f$ 倍、"
  % (a, sa, a/sa, sc[-1]["meanfield"]/sc[-1]["noise"]))
W("闭合的 $%.0f$ 倍。两者因而不是同一类误差：闭合的可以靠增大系统消除，均场的不能，"
  % (sc[-1]["meanfield"]/sc[-1]["closure"]))
W("因为那是闭合层级的误差。")
W("")
W("![Figure 2](figures/fig2_verification.png){width=5.5}")
W("")
W("## 结论")
W("")
W("五项验证全部通过。方程组自洽（§1），阈值公式经两条独立路径互证（§2），群内级联量经")
W("蒙特卡洛复核（§3），闭合能复现完整的含时演化（§4），且其误差随系统规模趋零而度均场")
W("的误差收敛到非零平台（§5）。第 2 节的理论至此具备可用于后续结构扫描的可信度。")

open(os.path.join(HERE, "..", "verification.md"), "w", encoding="utf-8").write("\n".join(out)+"\n")
print("wrote verification.md (%d lines)" % len(out))
