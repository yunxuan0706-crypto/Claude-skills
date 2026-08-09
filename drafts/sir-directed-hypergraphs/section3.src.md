# III. 并发型边基房室方程

## A. 空腔构造与闭合假设

成对网络上的边基房室模型（EBCM）以"沿一条边尚未传来感染"的概率为核心变量，把 SIR 的含时演化压缩到少数常微分方程 \cite{volz2008,miller2012ebcm}。移到有向超图上有两重并发需要处理：一个头节点同时暴露于 $k_{\rm in}$ 条超边，而每条超边内又有 $\tau$ 个尾成员同时向它施压。前一重可以因子化，后一重不能——本小节说明界线在哪里。

固定一个**检验节点** $u$，人为令其永久保持易感（空腔节点），并取一条以 $u$ 为头成员的超边 $e$，即 $u\in H(e)$。定义边基变量

$$
\Phi(t)=\Pr\bigl[e\ \text{到}\ t\ \text{时刻尚未向}\ u\ \text{传递感染}\bigr]\label{eq:phi}
$$

由 II.B 的速率规则，$e$ 在时刻 $s$ 以速率 $\beta g(n_e(s))$ 向 $u$ 施压，故

$$
\Phi(t)=\Bigl\langle\exp\Bigl[-\beta\int_0^{t}g\bigl(n_e(s)\bigr)\,\mathrm{d}s\Bigr]\Bigr\rangle\label{eq:phiexp}
$$

主线取 $\theta=1$，此时 $g(n_e)=\mathbb{1}[n_e\ge1]$，指数上的积分即 $e$ 的**激活时长**，也就是其各尾成员感染期之**并**的测度。并的测度不等于各成员贡献之和，因而

$$
\Phi(t)\neq\prod_{v\in T(e)}\Phi_v(t)\label{eq:noprod}
$$

朴素的逐成员乘积闭合在此失效。这与 II.D 中破坏命题 2 假设 (ii) 的是同一件事——超边级激活：一条超边只有一个激活时钟，它由全体尾成员共同驱动，又被全体头成员共享。

闭合建立在两条假设上。

**(H1) 局域树状。** 位形模型系综在 $N\to\infty$ 时局域收敛于树，任一有限邻域内出现回路的概率为 $O(1/N)$。据此 $u$ 的 $k_{\rm in}$ 条入超边的上游分支互不相交，第一重并发因子化：

$$
S(t)=(1-\varepsilon)\,\psi_{\rm in}\bigl(\Phi(t)\bigr),\qquad
\psi_{\rm in}(x)=\sum_{k}P(k_{\rm in}=k)\,x^{k}\label{eq:S}
$$

其中 $\varepsilon$ 为初始均匀感染比例。

**(H2) 尾内无相互作用。** 由 $T(e)\cap H(e)=\varnothing$，超边 $e$ 不在其尾成员之间传递感染；结合 (H1)，$T(e)$ 的 $\tau$ 个成员处在互不相交的上游分支，其感染时刻因而相互独立。

(H2) 是有向结构特有的。无向超图中同一条超边内的节点相互施加影响，"尾成员相互独立"根本无从谈起；正是尾与头的分离，使第二重并发退化为"$\tau$ 个独立个体驱动同一个共享时钟"，从而可由一个低维状态刻画。

还需注意一处**度偏置**。沿 $e$ 的一个尾槽到达的节点 $v$ 是按 $k_{\rm out}$ 加权抽取的，而决定 $v$ 何时被感染的是它的 $k_{\rm in}$ 条入超边。由 $T\cap H=\varnothing$，$e$ 本身不属于 $v$ 的入超边，因而**不出现无向 EBCM 中的余度相减**——偏置变量与计数变量是两个不同的量。据此定义尾侧生成函数

$$
\psi_{\rm tail}(x)=\frac{1}{\langle k_{\rm out}\rangle}\sum_{k_{\rm in},k_{\rm out}}k_{\rm out}\,P(k_{\rm in},k_{\rm out})\,x^{k_{\rm in}}\label{eq:psitail}
$$

$\psi_{\rm tail}$ 显式依赖于**联合**分布 $P(k_{\rm in},k_{\rm out})$，而 \eqref{eq:S} 的 $\psi_{\rm in}$ 只依赖于其边际。这是 II.A 的入出度相关 $r_{io}$ 进入动力学的解析通道，也是下文阈值公式中唯一承载 $r_{io}$ 的位置。

## B. 方程组

由 (H2)，$e$ 的 $\tau$ 个尾成员相互独立，其联合状态可由计数 $(a,b,c)$ 概括：$a$ 个易感、$b$ 个感染、$c$ 个康复，$a+b+c=\tau$。定义**未传递子概率**

$$
x_{abc}(t)=\Pr\bigl[e\ \text{尚未向}\ u\ \text{传递}\ \ \wedge\ \ T(e)\ \text{处于状态}\ (a,b,c)\bigr]\label{eq:xabc}
$$

于是 $\Phi=\sum_{a+b+c=\tau}x_{abc}$，而

$$
\Phi_A(t)=\sum_{b\ge\theta}x_{abc}(t)\label{eq:phiA}
$$

为"$e$ 处于激活态且尚未传递"的概率。

尾成员的感染是一个时间非齐次的 Markov 过程。由 \eqref{eq:S} 与 \eqref{eq:psitail}，一个尾成员到 $t$ 时刻仍易感的概率为 $(1-\varepsilon)\psi_{\rm tail}(\Phi(t))$，其瞬时风险率因而是

$$
h(t)=-\frac{\mathrm{d}}{\mathrm{d}t}\ln\bigl[(1-\varepsilon)\psi_{\rm tail}(\Phi)\bigr]
=\beta\,\Phi_A\,\frac{\psi'_{\rm tail}(\Phi)}{\psi_{\rm tail}(\Phi)}\label{eq:h}
$$

末一步用了下面的 \eqref{eq:phidot}。向 $u$ 的传递以速率 $\beta\,\mathbb{1}[b\ge\theta]$ 把概率移出"未传递"类，即一个**带湮灭的 Markov 链**。据此

$$
\dot x_{abc}=h\bigl[(a+1)x_{a+1,b-1,c}-a\,x_{abc}\bigr]
+\mu\bigl[(b+1)x_{a,b+1,c-1}-b\,x_{abc}\bigr]
-\beta\,\mathbb{1}[b\ge\theta]\,x_{abc}\label{eq:master}
$$

（越界下标记为零）。对 \eqref{eq:master} 求和，$h$ 项与 $\mu$ 项各自守恒而相消，只余

$$
\dot\Phi=-\beta\,\Phi_A\label{eq:phidot}
$$

与 \eqref{eq:h} 自洽。终态由 $\dot R=\mu I$、$I=1-S-R$ 与 \eqref{eq:S} 闭合。初始条件为 $x_{a,\tau-a,0}(0)=\binom{\tau}{a}(1-\varepsilon)^{a}\varepsilon^{\tau-a}$，其余为零。

方程 \eqref{eq:S} 与 \eqref{eq:master}–\eqref{eq:phidot} 构成一个 $\binom{\tau+2}{2}+1$ 维闭系统：$\tau=2$ 时 7 个方程，$\tau=3$ 时 11 个，维数与系统规模 $N$ 无关。

**$\tau=1$ 的退化。** 此时只有三个变量，\eqref{eq:master} 化为

$$
\dot x_{100}=-h\,x_{100},\qquad
\dot x_{010}=h\,x_{100}-(\mu+\beta)\,x_{010},\qquad
\dot x_{001}=\mu\,x_{010}\label{eq:tau1}
$$

这正是成对 EBCM 的 $(\theta_S,\theta_I,\theta_R)$ 三元组 \cite{miller2012ebcm}，唯一的替换是把余度生成函数换成 \eqref{eq:psitail} 的 $\psi_{\rm tail}$。

**一条精确恒等式。** 全部尾成员易感蕴含 $e$ 从未激活，因而必未传递，故

$$
x_{\tau00}(t)=\bigl[(1-\varepsilon)\,\psi_{\rm tail}(\Phi(t))\bigr]^{\tau}\label{eq:identity}
$$

对一切 $t$ 成立。\eqref{eq:identity} 不是额外假设，而是 \eqref{eq:master} 的推论：对右端求导得 $-\tau h\,x_{\tau00}$，与 \eqref{eq:master} 在 $(a,b,c)=(\tau,0,0)$ 处逐项相同。它为数值积分提供了一条无需仿真的自检。

**$\theta\ge2$。** \eqref{eq:master} 中只有湮灭项的指示函数依赖于 $\theta$，故方程组对任意 $\theta$ 形式不变。但下一小节的阈值分析不再适用：$\theta\ge2$ 时单个感染尾成员不足以激活超边，\eqref{eq:master} 在无病定态附近的线性化没有增长模式，这与 II.B 所述 bootstrap 图像一致，须按第 IV 节的 $k$-core 途径处理。

## C. 爆发阈值

令

$$
\kappa\equiv\psi'_{\rm tail}(1)=\frac{\langle k_{\rm in}k_{\rm out}\rangle}{\langle k_{\rm out}\rangle}\label{eq:kappa}
$$

$\varepsilon\to0$ 时无病定态为 $x_{\tau00}=1$、其余为零、$\Phi=1$。一阶量只有 $p\equiv x_{\tau-1,1,0}$ 与 $q\equiv x_{\tau-1,0,1}$，$b\ge2$ 的态为二阶。由 \eqref{eq:h} 有 $h\simeq\beta\kappa p$，代入 \eqref{eq:master} 得

$$
\dot p=\bigl[\tau\beta\kappa-(\mu+\beta)\bigr]\,p\label{eq:linear}
$$

$q$ 不反馈到 $p$，故阈值由 $\tau\beta\kappa=\mu+\beta$ 单独给出。以 $\lambda=\beta/\mu$ 表示，

$$
\lambda_c=\frac{1}{\tau\kappa-1},\qquad \tau\kappa>1\label{eq:lc}
$$

$\tau\kappa\le1$ 时不存在有限阈值。

**分支解释自洽。** 一个感染的尾成员使其超边激活，超边在该成员的感染期内以速率 $\beta$ 向每个头成员施压，故对给定头成员的传递概率为 $T=\beta/(\beta+\mu)$。从尾侧计数得 $R_0=\tau\kappa T$；从头侧计数，新感染节点的出度按 $k_{\rm in}$ 偏置，得 $R_0=\eta T\langle k_{\rm in}k_{\rm out}\rangle/\langle k_{\rm in}\rangle$。由握手关系 \eqref{eq:handshake} 有 $\langle k_{\rm out}\rangle/\langle k_{\rm in}\rangle=\tau/\eta$，两式因而恒等：

$$
\tau\,\frac{\langle k_{\rm in}k_{\rm out}\rangle}{\langle k_{\rm out}\rangle}
=\eta\,\frac{\langle k_{\rm in}k_{\rm out}\rangle}{\langle k_{\rm in}\rangle}\label{eq:twocounts}
$$

两种数法给出同一个 $R_0$，且 $R_0=1$ 与 \eqref{eq:lc} 等价。

**$\tau=\eta=1$ 的退化。** \eqref{eq:lc} 化为 $\lambda_c=\bigl[\langle k_{\rm in}k_{\rm out}\rangle/\langle k_{\rm out}\rangle-1\bigr]^{-1}$，即有向随机图上 SIR 的已知阈值 \cite{boguna2005directed,meyers2006directed}。需要强调，令 $k_{\rm in}=k_{\rm out}$ 并**不**回到无向网络的结果：无向情形的余度相减源于"来路那条边不可再用"，而有向情形下来路超边根本不在 $v$ 的入超边之列。\eqref{eq:lc} 分母中的 $-1$ 来自 $T=\beta/(\beta+\mu)$ 里的 $\beta$，与余度无关。两者是不同的结构，不应混为一谈。

**阈值对 $r_{io}$ 的依赖。** 由 $\langle k_{\rm in}k_{\rm out}\rangle=\langle k_{\rm in}\rangle\langle k_{\rm out}\rangle+r_{io}\sigma_{k_{\rm in}}\sigma_{k_{\rm out}}$ 与握手关系 \eqref{eq:handshake}，

$$
\tau\kappa=\tau\langle k_{\rm in}\rangle+\frac{N}{M}\,r_{io}\,\sigma_{k_{\rm in}}\sigma_{k_{\rm out}}\label{eq:taukappa}
$$

于是

$$
\lambda_c=\Bigl[\tau\langle k_{\rm in}\rangle+\tfrac{N}{M}\,r_{io}\,\sigma_{k_{\rm in}}\sigma_{k_{\rm out}}-1\Bigr]^{-1}\label{eq:lcrio}
$$

在两条边际度序列固定时 $\lambda_c$ 是 $r_{io}$ 的严格减函数：正的入出度相关促进传播。这是第 VI 节 $r_{io}$ 扫描的定量预言。

**一个否定性推论。** \eqref{eq:lc} 只依赖于阶数 $\tau$ 与双度序列（经由 $\kappa$）。II.C 的重叠 $\alpha^{ab}$ 在双度序列固定时仍可自由调节，因而**在树状闭合内 $\lambda_c$ 与 $\alpha$ 无关**。这不是疏漏而是可检验的预言：重叠恰好度量超边之间共享节点的程度，也就是 (H1) 所排除的短回路。第 VI 节若测得阈值随 $\alpha$ 移动，其幅度即为树状闭合失效的定量刻度；若测不到，则 $\alpha$ 的作用只体现在阈值以外的量上。无论哪一种结果，都是有内容的。

## D. 与均场闭合的比较

同一模型的度均场闭合是：设节点 $v$ 的感染概率为 $i_v$，把其入超边的激活概率取为系综平均 $\Theta=1-(1-\phi)^{\tau}$（$\theta=1$），其中 $\phi=\sum_v k_{\rm out}(v)i_v\big/\sum_v k_{\rm out}(v)$ 为尾槽被感染节点占据的比例，则

$$
\dot s_v=-\beta\,k_{\rm in}(v)\,\Theta\,s_v,\qquad
\dot i_v=\beta\,k_{\rm in}(v)\,\Theta\,s_v-\mu\,i_v\label{eq:mf}
$$

线性化给出 $\dot\phi=(\beta\tau\kappa-\mu)\phi$，即

$$
\lambda_c^{\rm MF}=\frac{1}{\tau\kappa}\label{eq:lcmf}
$$

与 \eqref{eq:lc} 相比，均场丢掉的恰是分母中的 $-1$，其来源可以精确指出：均场把一个尾成员在整个感染期内向同一头成员的传递概率记作 $\beta/\mu=\lambda$，而实际为 $\beta/(\beta+\mu)=\lambda/(1+\lambda)$——超边一旦向某头成员传递即对它用尽，均场没有这一条件。相对误差

$$
\frac{\lambda_c}{\lambda_c^{\rm MF}}=\frac{\tau\kappa}{\tau\kappa-1}\label{eq:gain}
$$

在 $\tau\kappa\to1^{+}$ 时发散：**均场的误差在阈值邻域最大**，而这正是最需要理论的区域。这与成对网络上 EBCM 相对均场的增益同源 \cite{volz2008,miller2012ebcm}，本文把它推广到有向超图并给出闭式。

需要说明，\eqref{eq:mf} 是我们为本模型构造的均场对照，而非 \cite{li2024directed} 中模型的复述——后者是社会传播而非 SIR，其方向性由单一标量强度参数调节。此处比较的是**闭合层级**，不是两项工作的模型。

## E. 数值校验

校验分三层：闭合内部的自洽、闭合对精确仿真的复现、以及两条结构预言。全部脚本见 `tools/ebcm_directed.py`；仿真为 II.B 所定义连续时间 Markov 过程的精确 Gillespie 抽样 \cite{gillespie1977}，取 $\mu=1$。

**内部自洽。** 恒等式 \eqref{eq:identity} 的残差以 Euler 步长线性收敛：$\mathrm{d}t$ 由 $0.008$ 逐次减半至 $0.001$ 时，最大残差由 $7.13\times10^{-4}$ 降至 $8.89\times10^{-5}$，相邻比值三次均为 $2.00$。残差因而是积分误差而非闭合的破缺。两种分支计数 \eqref{eq:twocounts} 在 $(\tau,\eta)=(2,2),(3,1),(1,3),(2,3)$ 上的相对差均不超过 $2\times10^{-16}$。$\kappa$ 在 $4000$ 次保度双边交换下的改变量恰为零，与 III.C 的否定性推论一致。

阈值公式 \eqref{eq:lc} 与方程组自身的失稳点独立比对：对 \eqref{eq:master} 由 $\varepsilon=10^{-12}$ 出发积分，以终态是否宏观为判据对 $\lambda$ 二分，全过程不使用 $\kappa$。所得阈值随积分窗口 $t_{\max}$ 单调趋近预言值，

| $t_{\max}$ | 200 | 400 | 1600 | 6400 |
|---|---|---|---|---|
| 相对偏差 | $2.7\times10^{-2}$ | $1.1\times10^{-2}$ | $1.6\times10^{-3}$ | $4.8\times10^{-5}$ |

偏差来自阈值邻域增长率趋零、有限窗口内尚未展开，而非公式本身。$\tau=\eta=1$ 时 \eqref{eq:lc} 与有向随机图的已知阈值逐位相同。

**对精确仿真的复现。** 取 $\tau=\eta=2$、$\lambda=1.6\lambda_c$、初始感染比例 $\varepsilon=0.01$，比较 $S(t)$ 在 $t\le12$ 上的平均偏差：

| $N$ | EBCM | 均场 | 终态 $R(\infty)$（仿真 / EBCM / 均场） |
|---|---|---|---|
| 1000 | 0.0355 | 0.199 | 0.362 / 0.407 / 0.586 |
| 4000 | 0.0073 | 0.165 | 0.410 / 0.413 / 0.585 |

两个尺度是独立抽取的结构。EBCM 的偏差随 $N$ 增大约缩小五倍，终态误差由 $12\%$ 降到 $0.8\%$——树状闭合的误差应当且确实随 $N$ 消失；均场的偏差则几乎不随 $N$ 变化，说明那是闭合层级的误差而非有限尺寸效应。

**阈值邻域。** \eqref{eq:gain} 预言均场的误差在 $\tau\kappa\to1^{+}$ 时发散。在 $N=8000$ 上沿 $\lambda/\lambda_c$ 扫描终态规模：

| $\lambda/\lambda_c$ | 0.7 | 0.9 | 1.0 | 1.2 | 1.6 | 2.4 |
|---|---|---|---|---|---|---|
| 仿真 | 0.0196 | 0.0449 | 0.0791 | 0.199 | 0.401 | 0.565 |
| EBCM | 0.0197 | 0.0469 | 0.0862 | 0.213 | 0.401 | 0.566 |
| 均场 | 0.0491 | 0.229 | 0.315 | 0.439 | 0.576 | 0.689 |

EBCM 在整条扫描线上跟随仿真，包括阈值处本身；均场在 $\lambda=0.9\lambda_c$ 处高估五倍、在 $\lambda_c$ 处高估四倍，深超临界区才收敛到 $1.2$ 倍。这正是 \eqref{eq:gain} 所预言的形状：均场把阈值压低到 $\lambda_c^{\rm MF}=0.75\lambda_c$，因而在真实的亚临界区已经给出宏观爆发。第 I 节所称"边基闭合相对均场的增益"由此定量化——增益不是均匀的百分数，而是集中在阈值邻域。

**结构预言。** 保持两条边际度序列不变、只改变 $k_{\rm in}$ 与 $k_{\rm out}$ 在节点上的配对，$\tau=\eta=2$：

| $r_{io}$ | $-0.795$ | $-0.394$ | $-0.027$ | $+0.387$ | $+0.795$ |
|---|---|---|---|---|---|
| $\kappa$ | 1.197 | 1.596 | 2.013 | 2.372 | 2.790 |
| $\lambda_c$ | 0.717 | 0.456 | 0.330 | 0.267 | 0.218 |

阈值在 $r_{io}$ 的可及范围内变动三倍以上，与 \eqref{eq:lcrio} 一致。这是第 VI 节 $r_{io}$ 扫描的解析基线；同时须记住 II.D 的宇称结论——$r_{io}$ 移动阈值，却不产生方向对称性破缺，二者并不矛盾。
