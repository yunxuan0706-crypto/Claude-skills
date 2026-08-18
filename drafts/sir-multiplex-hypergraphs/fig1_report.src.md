# 多重超图上的 SIR：群体闭合方程与含时演化

# 一、模型与群体闭合方程

![图 1　模型示意。(a) 多重超图：节点 $u$（深灰）在层 1（$m_1=3$，实线）属于 $e_1,e_2$ 两个群体、在层 2（$m_2=4$，虚线）属于 $e_3$，故 $k^{(1)}(u)=2$、$k^{(2)}(u)=1$。(b) 群内级联：朴素计数下一个种子直接感染 $m-1$ 个成员、群体激活期仅一段，产出 $(m-1)T$；被感染成员延长激活期（红色 extension），级联产出 $C(m,\lambda,\theta)>(m-1)T$。(c) 群体层面分支：经层 $a$ 到达的新感染节点点燃其各层 $b$ 群体，每个产出 $C$ 个新感染，层内（$a=b$）须扣除来路群体的余度相减 $-\delta_{ab}$，得 $N_{ab}=C(X_{ab}-\delta_{ab})$，其中 $X_{ab}=\langle k^{(a)}k^{(b)}\rangle/\langle k^{(a)}\rangle$。](figures/fig_model.png){width=6.5}

## 1.1 多重超图与记号

节点集为 $V$，$N=|V|$。层 $a=1,\dots,M$，层 $a$ 是 $V$ 上的超图
$\mathcal{H}_a=(V,E_a)$，主线取同质基数 $|e|=m_a$。节点的层度
$k^{(a)}(v)=|\{e\in E_a:v\in e\}|$，写成向量 $\mathbf{k}(v)$。引入两个生成函数：

$$
\Psi(\mathbf{x})=\sum_{\mathbf{k}}P(\mathbf{k})\prod_{c}x_c^{k_c}\label{eq:psi}
$$

$$
\psi_a(\mathbf{x})=\frac{1}{\langle k^{(a)}\rangle}\sum_{\mathbf{k}}k_a\,P(\mathbf{k})\prod_c x_c^{\,k_c-\delta_{ac}}\label{eq:psia}
$$

$\Psi$ 描述节点侧，$\psi_a$ 描述沿层 $a$ 的一个成员槽到达某节点时该节点的余度分布，
指数中的 $-\delta_{ac}$ 即余度相减：来路那个群体不再计入。联合层度分布
$P(\mathbf{k})$ 是本模型的结构输入，层间的一切效应最终只经由 $\Psi$ 与 $\psi_a$ 的
交叉项进入动力学。

## 1.2 SIR 动力学

节点状态 $x_v\in\{S,I,R\}$ 跨层共享——个体在任一渠道被感染即在所有渠道中均为感染
态。记层 $a$ 的群体 $e$ 中已感染成员数为 $n_e(t)=|e\cap I(t)|$。当 $n_e\ge\theta_a$
时群体激活，向 $e$ 中每个易感成员各以速率 $\beta_a$ 独立传递，否则静默；感染者以速率
$\mu$ 康复。易感节点 $u$ 的总风险率

$$
\Lambda_u(t)=\sum_{a}\ \sum_{e\in E_a:\,u\in e}\beta_a\,\mathbb{1}[n_e(t)\ge\theta_a]\label{eq:rate}
$$

以 $\mu$ 归一时间，$\lambda_a=\beta_a/\mu$，主线取 $\theta_a=1$。初始时刻每个节点独立地
以概率 $\varepsilon$ 处于感染态、其余易感，$\varepsilon$ 即均匀初始感染比例。与有向超图不同，无向
群体中每个成员既是传染源也是传染靶：一个成员被感染后群体激活，可感染同群其他成员，
后者被感染又延长群体的激活时长，形成群内级联。本模型的新内容均源于此。

## 1.3 空腔构造与群体计数方程组

固定一个检验节点 $u$，令其永久易感（空腔节点）。取一个 $u$ 所属的层 $a$ 群体 $e$，
定义边基变量

$$
\Phi_a(t)=\Pr\bigl[e\ \text{到}\ t\ \text{时刻尚未向}\ u\ \text{传递}\bigr]\label{eq:phi}
$$

在局域树状假设下（位形模型系综 $N\to\infty$ 时局域收敛于超树，回路概率 $O(1/N)$），
$u$ 所属各群体的上游分支互不相交，节点侧因而闭合为

$$
S(t)=(1-\varepsilon)\,\Psi\bigl(\Phi_1(t),\dots,\Phi_M(t)\bigr)\label{eq:S}
$$

以 $e$ 除 $u$ 外的 $m_a-1$ 个成员的计数 $(s,i,r)$（$s+i+r=m_a-1$）为仓室，定义
未传递子概率 $x^{(a)}_{sir}(t)$，则 $\Phi_a=\sum x^{(a)}_{sir}$，激活且未传递的部分为

$$
\Phi^{A}_a(t)=\sum_{i\ge\theta_a}x^{(a)}_{sir}(t)\label{eq:phiA}
$$

一个成员经由其他群体到 $t$ 仍易感的概率为 $(1-\varepsilon)\psi_a(\boldsymbol{\Phi})$，
故其外部风险率

$$
h_a(t)=-\frac{\mathrm{d}}{\mathrm{d}t}\ln\bigl[(1-\varepsilon)\psi_a(\boldsymbol{\Phi})\bigr]
=\sum_{c}\beta_c\,\Phi^{A}_c\,\frac{\partial_c\psi_a(\boldsymbol{\Phi})}{\psi_a(\boldsymbol{\Phi})}\label{eq:h}
$$

一个易感成员的总感染率是外部与群内之和 $h_a+\beta_a\mathbb{1}[i\ge\theta_a]$；群体向
$u$ 的传递以速率 $\beta_a\mathbb{1}[i\ge\theta_a]$ 把概率移出「未传递」类。据此

$$
\begin{aligned}
\dot x^{(a)}_{sir}=\;&\bigl(h_a+\beta_a\mathbb{1}[i-1\ge\theta_a]\bigr)(s+1)\,x^{(a)}_{s+1,i-1,r}
-\bigl(h_a+\beta_a\mathbb{1}[i\ge\theta_a]\bigr)s\,x^{(a)}_{sir}\\
&+\mu\bigl[(i+1)x^{(a)}_{s,i+1,r-1}-i\,x^{(a)}_{sir}\bigr]
-\beta_a\mathbb{1}[i\ge\theta_a]\,x^{(a)}_{sir}
\end{aligned}\label{eq:master}
$$

（越界下标记为零。）对 \eqref{eq:master} 求和，$h_a$ 项、群内项与 $\mu$ 项各自守恒
相消，只余

$$
\dot\Phi_a=-\beta_a\,\Phi^{A}_a\label{eq:phidot}
$$

与 \eqref{eq:h} 自洽。节点侧由 \eqref{eq:S} 与 $\dot R=\mu I$、$I=1-S-R$ 闭合。初始
条件 $x^{(a)}_{s,i,0}(0)=\binom{m_a-1}{i}(1-\varepsilon)^{s}\varepsilon^{i}$，其余为零，
$\Phi_a(0)=1$。方程组共 $\sum_a\binom{m_a+1}{2}+1$ 维（各层的群体计数加上节点侧的单个
独立标量 $R$，其中 $S$ 由 \eqref{eq:S} 决定、$I=1-S-R$ 导出），与系统规模 $N$ 无关；
$M=2$、$m=(3,4)$ 时为 $16+1=17$ 个方程，积分一次即给出完整的 $S(t)$、$I(t)$、$R(t)$
与终态规模。

## 1.4 爆发阈值

以「感染经由哪一层传来」为型，经由层 $a$ 感染的节点其层 $b$ 群体数的期望为
$\langle k^{(a)}k^{(b)}\rangle/\langle k^{(a)}\rangle-\delta_{ab}$，每个被点燃的群体产出
$C(m_b,\lambda_b,\theta_b)$ 个新感染成员（$C$ 为孤立群体内的群内级联期望，由一个小型
连续时间 Markov 链的递推给出，$m=2$ 时退化为成对传递概率）。故多型分支矩阵

$$
\mathsf{N}_{ab}=C(m_b,\lambda_b,\theta_b)\Bigl[\frac{\langle k^{(a)}k^{(b)}\rangle}{\langle k^{(a)}\rangle}-\delta_{ab}\Bigr]\label{eq:N}
$$

爆发阈值由谱半径条件 $\rho(\mathsf{N})=1$ 确定。图 2 的控制参数 $\lambda$ 即以此阈值
$\lambda_c$ 定标。

# 二、含时演化：闭合、均场与精确仿真

作为对照，同一模型的度均场取群体激活概率为系综平均
$\Theta_a=1-(1-\phi_a)^{m_a-1}$，其中 $\phi_a$ 为层 $a$ 成员槽被感染节点占据的比例，

$$
\dot s_\mathbf{k}=-s_\mathbf{k}\sum_a k_a\beta_a\Theta_a,\qquad
\dot i_\mathbf{k}=s_\mathbf{k}\sum_a k_a\beta_a\Theta_a-\mu\, i_\mathbf{k}\label{eq:mf}
$$

均场丢弃了群内相关、邻域相关与群内级联，是评价群体闭合的基线。

图 2 给出 $M=2$、$m=(3,4)$、$N=6000$、$\lambda=1.6\lambda_c$ 下的完整时间演化，联合
层度分布 $P(\mathbf{k})$ 在 $(2,2)$ 与 $(3,3)$ 上各占一半，$\lambda_c=0.0993$，
$\varepsilon=0.02$。群体闭合 \eqref{eq:master} 与精确仿真在 $S(t)$ 与 $I(t)$ 上的最大
偏差分别为 $0.0039$ 与 $0.0011$，均小于图中线宽；感染峰高相差 $0.9\%$。度均场
\eqref{eq:mf} 则把感染峰值由仿真的 $0.086$ 抬到 $0.139$（$1.6$ 倍），并把峰位由
$t\approx4.4$ 提前到 $t\approx3.6$。群体闭合复现了完整的含时演化，度均场在同一参数下
显著失真。

![图 2　多重超图 SIR 的含时演化，$N=6000$、$M=2$、$m=(3,4)$、$\lambda=1.6\lambda_c$、初始感染比例 $\varepsilon=0.02$。黑线为 $400$ 次独立 Gillespie 实现的均值，灰带为 $95\%$ 置信区间，最宽处约与线宽相当；蓝虚线为群体闭合 \eqref{eq:master}，橙点线为度均场 \eqref{eq:mf}。(a) 易感比例 $S(t)$；(b) 感染比例 $I(t)$。](figures/fig1_trajectory.png){width=6.2}

# 三、数值校验

仿真为 \eqref{eq:rate} 所定义连续时间 Markov 过程的精确 Gillespie 抽样，取 $\mu=1$。

**内部自洽。** 方程组自带两条无需仿真的自检。其一为求和规则：对 \eqref{eq:master} 逐态
求和应严格给出 \eqref{eq:phidot}，残差处于机器精度（$10^{-20}$ 量级），说明它是方程组
的代数恒等式。其二为全体成员易感的恒等式
$x^{(a)}_{m_a-1,0,0}=[(1-\varepsilon)\psi_a(\boldsymbol{\Phi})]^{m_a-1}$，其残差随积分
步长按 $\mathrm{d}t^4$ 收敛（$\mathrm{d}t$ 由 $0.02$ 逐次减半时残差为
$1.3\times10^{-10}$、$8.0\times10^{-12}$、$5.0\times10^{-13}$、$3.1\times10^{-14}$，
相邻比约 $16$），与四阶 Runge–Kutta 的阶一致，说明该恒等式是方程组的严格推论、残差
纯属积分误差。

**阈值自洽。** 阈值公式 \eqref{eq:N} 与方程组自身在无病定态的线性失稳独立比对：将
\eqref{eq:master} 作数值 Jacobian、对最大实部特征值二分求阈值，与 $\rho(\mathsf{N})=1$
在六组构型（含反相关度分布、非对称速率与群体规模、三层、以及成对退化情形
$\lambda_c=1/3$）上的相对差均不超过 $1\times10^{-8}$，即数值 Jacobian 的有限差分精度。
两条路径不共享代码；若以 $(m-1)T$ 取代群内级联 $C$，该一致性立即消失。

**含时复现。** 如图 2 所示，群体闭合与精确仿真在整个时间演化上逐点吻合，$S(t)$ 与
$I(t)$ 的最大偏差均小于线宽（$0.0039$ 与 $0.0011$），感染峰高相差 $0.9\%$；对峰邻域作
抛物线插值以消除网格离散后，闭合与仿真的峰位分别为 $t=4.38$ 与 $4.37$，相差约
$0.01$；度均场则系统性偏离。偏差随系统规模的标度进一步表明，群体闭合的误差随 $N$ 增大单调趋零（至
$N=8000$ 已落入蒙特卡洛噪声），而度均场的误差收敛到非零平台——两者不是同一类误差：
前者可靠增大系统消除，后者是闭合层级的误差，不能。

**独立复算。** 上述结果均由与原推导不共享路径的方式复核：$m=2$（成对）退化情形下
方程组可解析积出终态不动点，与 RK4 积分的 $R(\infty)$ 相差 $3\times10^{-14}$；群内级联
$C$ 由孤立群体的直接蒙特卡洛复核，全部十五组参数（$m=2,\dots,6$ 与 $\lambda=0.5,1,2$）
的标准化偏差均不超过 $2.5\sigma$、与无系统偏差相容；阈值本身亦由亚临界终态外推独立
测得，与 \eqref{eq:N} 的解析值相容。
