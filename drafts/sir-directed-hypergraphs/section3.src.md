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

末一步用了下面的 \eqref{eq:phidot}。向 $u$ 的传递以速率 $\beta\,\mathbb{1}[b\ge\theta]$ 把概率移出"未传递"类。等价地，在 $(a,b,c)$ 之外添一个"已传递"的吸收态，则整体是一条标准的连续时间 Markov 链，而 $x_{abc}$ 是其未被吸收部分的**子概率**，$\sum_{abc}x_{abc}=\Phi\le1$；概率论中这称为以速率 $\beta\,\mathbb{1}[b\ge\theta]$ **消灭**（killing）的 Markov 链，其半群为次随机的，而 \eqref{eq:phiexp} 正是它的 Feynman–Kac 表示。据此

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

这正是成对 EBCM 的三元组 $(\Phi_S,\Phi_I,\Phi_R)$ \cite{miller2012ebcm}——依次为伙伴处于易感、已感染但尚未传递、已康复且未曾传递的概率——唯一的实质替换是把余度生成函数换成 \eqref{eq:psitail} 的 $\psi_{\rm tail}$。需注意记号差异：本文的 $\Phi$ 对应该文中三者之和 $\Theta$，而本文的 $x_{abc}$ 对应其 $\Phi$ 分量。

**$\tau=2$ 的显式形式。** 六个状态按 $(a,b,c)$ 记为 $x_{200},x_{110},x_{101},x_{020},x_{011},x_{002}$，\eqref{eq:master} 展开为

$$
\dot x_{200}=-2h\,x_{200},\qquad
\dot x_{110}=2h\,x_{200}-(h+\mu+\beta)\,x_{110},\qquad
\dot x_{101}=\mu\,x_{110}-h\,x_{101}\label{eq:tau2a}
$$

$$
\dot x_{020}=h\,x_{110}-(2\mu+\beta)\,x_{020},\qquad
\dot x_{011}=h\,x_{101}+2\mu\,x_{020}-(\mu+\beta)\,x_{011},\qquad
\dot x_{002}=\mu\,x_{011}\label{eq:tau2b}
$$

前三式对应尾集中至多一个成员曾被感染的情形，后三式对应两个成员；$\Phi_A=x_{110}+x_{020}+x_{011}$。逐式相加，$h$ 项与 $\mu$ 项成对相消，只余 $-\beta\Phi_A$，即 \eqref{eq:phidot}。手写的这六个方程与按 \eqref{eq:master} 自动生成的方程组，其积分在 $t=14$ 处相差 $6\times10^{-17}$。

**一条精确恒等式。** 全部尾成员易感蕴含 $e$ 从未激活，因而必未传递，故

$$
x_{\tau00}(t)=\bigl[(1-\varepsilon)\,\psi_{\rm tail}(\Phi(t))\bigr]^{\tau}\label{eq:identity}
$$

对一切 $t$ 成立。\eqref{eq:identity} 不是额外假设，而是 \eqref{eq:master} 的推论：对右端求导得 $-\tau h\,x_{\tau00}$，与 \eqref{eq:master} 在 $(a,b,c)=(\tau,0,0)$ 处逐项相同。它为数值积分提供了一条无需仿真的自检。

**求解到什么程度。** 有必要把"求解"的含义说清楚。$\theta=1$ 时，\eqref{eq:S} 与 \eqref{eq:master}–\eqref{eq:phidot} 构成一个封闭的低维系统，其积分给出完整的 $S(t)$、$I(t)$、$R(t)$ 与终态规模，精度由 III.E 逐项校验；爆发阈值 \eqref{eq:lc} 是闭式，$\tau=1$ 时终态亦有闭式 \eqref{eq:finalsize}。但 $\tau\ge2$ 时终态**没有**闭式，原因是结构性的而非技术性的：$\Phi_\infty$ 取决于各尾成员感染期之**并**的长度；$\tau=1$ 时这个并只有一段，其长度服从 $\mathrm{Exp}(\mu)$ 且与起始时刻无关，故可只用终态量闭合，而 $\tau\ge2$ 时并的长度依赖于各感染时刻的相对先后，终态量不足以确定它。这与 \eqref{eq:noprod} 是同一个障碍。取线性核 $g(n)=n$ 时传播事件逐（尾成员，头成员）对独立，闭式重新出现。$\theta\ge2$ 的爆发条件不属本节范围。

**$\theta\ge2$。** \eqref{eq:master} 中只有消灭项的指示函数依赖于 $\theta$，故方程组对任意 $\theta$ 形式不变。但下一小节的阈值分析不再适用：$\theta\ge2$ 时单个感染尾成员不足以激活超边，\eqref{eq:master} 在无病定态附近的线性化没有增长模式，这与 II.B 所述 bootstrap 图像一致，须按第 IV 节的 $k$-core 途径处理。

## C. 爆发阈值

令

$$
\kappa\equiv\psi'_{\rm tail}(1)=\frac{\langle k_{\rm in}k_{\rm out}\rangle}{\langle k_{\rm out}\rangle}\label{eq:kappa}
$$

$\varepsilon\to0$ 时无病定态为 $x_{\tau00}=1$、其余为零、$\Phi=1$。一阶量是恰有一个尾成员曾被感染的两个态，$p\equiv x_{\tau-1,1,0}$ 与 $q\equiv x_{\tau-1,0,1}$；判据是曾被感染的成员数 $b+c$，$b+c\ge2$ 的态为二阶——注意 $x_{\tau-2,1,1}$ 虽然 $b=1$ 也在其中。由两个种子规模积分读出各态随 $\varepsilon$ 的幂次，测得阶数逐态等于 $b+c$（$\tau=3$，见 `tools/check_order_counting.py`）。$\Phi_A$ 的一阶部分因而只有 $p$。由 \eqref{eq:h} 有 $h\simeq\beta\kappa p$，代入 \eqref{eq:master} 得

$$
\dot p=\bigl[\tau\beta\kappa-(\mu+\beta)\bigr]\,p\label{eq:linear}
$$

$q$ 不反馈到 $p$，故阈值由 $\tau\beta\kappa=\mu+\beta$ 单独给出。以 $\lambda=\beta/\mu$ 表示，

$$
\lambda_c=\frac{1}{\tau\kappa-1},\qquad \tau\kappa>1\label{eq:lc}
$$

$\tau\kappa\le1$ 时不存在有限阈值。

**分支解释自洽。** 一个感染的尾成员使其超边激活，超边在该成员的感染期内以速率 $\beta$ 向每个头成员施压，故对给定头成员的传递概率为 $T=\beta/(\beta+\mu)$。从尾侧计数得 $R_0=\tau\kappa T$；从头侧计数，新感染节点沿入边到达、其出度因而按 $k_{\rm in}$ 偏置，得 $R_0=\eta T\langle k_{\rm in}k_{\rm out}\rangle/\langle k_{\rm in}\rangle$——成对情形下 $\langle k_{\rm in}k_{\rm out}\rangle/\langle k_{\rm in}\rangle$ 正是有向图中"后继节点的平均出度"，亦即其再生数 \cite{allard2023review}。由握手关系 \eqref{eq:handshake} 有 $\langle k_{\rm out}\rangle/\langle k_{\rm in}\rangle=\tau/\eta$，两式因而恒等：

$$
\tau\,\frac{\langle k_{\rm in}k_{\rm out}\rangle}{\langle k_{\rm out}\rangle}
=\eta\,\frac{\langle k_{\rm in}k_{\rm out}\rangle}{\langle k_{\rm in}\rangle}\label{eq:twocounts}
$$

两种数法给出同一个 $R_0$，且 $R_0=1$ 与 \eqref{eq:lc} 等价。

**$\tau=\eta=1$ 的退化。** \eqref{eq:lc} 化为 $\lambda_c=\bigl[\langle k_{\rm in}k_{\rm out}\rangle/\langle k_{\rm out}\rangle-1\bigr]^{-1}$，即有向随机图上 SIR 的已知阈值 \cite{meyers2006directed,boguna2005directed,allard2023review}。\cite{meyers2006directed} 的半有向阈值条件在纯有向极限下化为临界传播率 $T_c=\langle k_{\rm in}\rangle/\langle k_{\rm in}k_{\rm out}\rangle$，\cite{allard2023review} 亦以再生数 $R_0=\langle k_{\rm in}k_{\rm out}\rangle/\langle k_{\rm in}\rangle$ 给出同一结果；由 $T=\lambda/(1+\lambda)$ 代入 \eqref{eq:lc} 得 $T_c=1/\kappa=\langle k_{\rm out}\rangle/\langle k_{\rm in}k_{\rm out}\rangle$，而 $\tau=\eta$ 时握手关系 \eqref{eq:handshake} 给出 $\langle k_{\rm in}\rangle=\langle k_{\rm out}\rangle$，两式因而恒等。需要强调，令 $k_{\rm in}=k_{\rm out}$ 并**不**回到无向网络的结果：无向情形的余度相减源于"来路那条边不可再用"，而有向情形下来路超边根本不在 $v$ 的入超边之列。\eqref{eq:lc} 分母中的 $-1$ 来自 $T=\beta/(\beta+\mu)$ 里的 $\beta$，与余度无关。两者是不同的结构，不应混为一谈。

**阈值对 $r_{io}$ 的依赖。** 由 $\langle k_{\rm in}k_{\rm out}\rangle=\langle k_{\rm in}\rangle\langle k_{\rm out}\rangle+r_{io}\sigma_{k_{\rm in}}\sigma_{k_{\rm out}}$ 与握手关系 \eqref{eq:handshake}，

$$
\tau\kappa=\tau\langle k_{\rm in}\rangle+\frac{N}{M}\,r_{io}\,\sigma_{k_{\rm in}}\sigma_{k_{\rm out}}\label{eq:taukappa}
$$

于是

$$
\lambda_c=\Bigl[\tau\langle k_{\rm in}\rangle+\tfrac{N}{M}\,r_{io}\,\sigma_{k_{\rm in}}\sigma_{k_{\rm out}}-1\Bigr]^{-1}\label{eq:lcrio}
$$

在两条边际度序列固定时 $\lambda_c$ 是 $r_{io}$ 的严格减函数：正的入出度相关促进传播。这一方向与成对有向图上的 SIR 一致——那里 $R_0$ 同样经 $\langle k_{\rm in}k_{\rm out}\rangle$ 依赖于入出度的协方差 \cite{allard2023review}；\eqref{eq:lcrio} 把该依赖推广到有向超图，并显式给出 $\tau$ 与 $N/M$ 的系数。需要提醒的是这一方向并非普适于一切传播机制：有向网络的阈值型级联中，正的入出度相关反而**提高**系统稳健性 \cite{xu2018directionality}。二者并不冲突——阈值模型按活跃邻居的比例判定激活，高入度抬高激活门槛；SIR 中高入度只增加暴露，故 $k_{\rm in}$ 与 $k_{\rm out}$ 同时偏大的节点在 SIR 中放大 $R_0$，在阈值级联中却既难被激活又本可广播。这是第 VI 节 $r_{io}$ 扫描的定量预言。

**一个否定性推论。** \eqref{eq:lc} 只依赖于阶数 $\tau$ 与双度序列（经由 $\kappa$）。更强地说，\eqref{eq:S} 与 \eqref{eq:psitail} 也只读取双度序列，故**整个闭合而不止阈值都看不见重叠 $\alpha$**。由 II.C 的求和规则，双度序列锁死的只是四个通道的重叠总量 $C_{ab}$，并不锁死 $\alpha^{ab}=C_{ab}/m_{ab}$ 本身，$\alpha$ 因而在固定度序列下仍是可动的自由度。这不是疏漏而是可检验的预言：重叠恰好度量超边之间共享节点的程度，也就是 (H1) 所排除的短回路。第 VI 节若测得阈值随 $\alpha$ 移动，其幅度即为树状闭合失效的定量刻度；若测不到，则 $\alpha$ 的作用只体现在阈值以外的量上。但须说明，可动不等于易动：单纯的保度双边交换只能把 $\alpha^{\parallel}$ 挪动 $0.005$–$0.006$（$N=150$ 与 $N=500$ 两处测量），在如此窄的区间上仿真的变化不可分辨，因而检验这一预言必须先有第 V 节的定向重连工具把 $\alpha$ 推到有意义的范围。极性失衡 $\Delta\alpha$ 是例外，它可被驱动到 $0.108$，故第 VI 节对它的上界测量不受此限。

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

校验分三层：闭合内部的自洽、闭合对精确仿真的复现、以及结构预言。脚本见 `tools/ebcm_directed.py` 与 `tools/section3_data.py`；仿真为 II.B 所定义连续时间 Markov 过程的精确 Gillespie 抽样 \cite{gillespie1977}，取 $\mu=1$。

**内部自洽。** 恒等式 \eqref{eq:identity} 的残差随 Euler 步长线性收敛：$\mathrm{d}t$ 由 $0.008$ 逐次减半至 $0.001$ 时，最大残差由 $7.13\times10^{-4}$ 降至 $8.89\times10^{-5}$，相邻比值三次均为 $2.00$。残差因而是积分误差而非闭合的破缺。两种分支计数 \eqref{eq:twocounts} 在 $(\tau,\eta)=(2,2),(3,1),(1,3),(2,3)$ 上的相对差不超过 $2\times10^{-16}$；\eqref{eq:lcrio} 与 \eqref{eq:lc} 在全部 $45$ 个位形模型实现上的相对差不超过 $4\times10^{-16}$，即经由 $r_{io}$ 的改写是恒等的。$\kappa$ 在 $4000$ 次保度双边交换下的改变量恰为零，与本节的否定性推论一致。

阈值公式 \eqref{eq:lc} 与方程组自身的失稳点独立比对：对 \eqref{eq:master} 由 $\varepsilon=10^{-12}$ 出发积分，以终态是否宏观为判据对 $\lambda$ 二分，全过程不使用 $\kappa$。所得阈值随积分窗口 $t_{\max}$ 单调趋近预言值——相对偏差在 $t_{\max}=200,400,1600,6400$ 处依次为 $2.7\times10^{-2}$、$1.1\times10^{-2}$、$1.6\times10^{-3}$、$4.8\times10^{-5}$。偏差来自阈值邻域增长率趋零、有限窗口内尚未展开，而非公式本身。

**含时复现。** 图 1 给出 $\tau=\eta=2$、$N=6000$、$\lambda=1.6\lambda_c$ 下的完整时间演化。并发型 EBCM 与仿真在 $S(t)$ 与 $I(t)$ 上的最大偏差分别为 $0.0031$ 与 $0.0008$，均小于图中线宽；峰位完全重合（$t=4.50$），峰高相差 $1.3\%$。均场则把峰值由 $0.059$ 抬到 $0.117$（约一倍），并把峰位由 $t=4.50$ 提前到 $t=3.75$。

![图 1　有向超图 SIR 的含时演化，$N=6000$、$\tau=\eta=2$、$\lambda=1.6\lambda_c$、初始感染比例 $\varepsilon=0.01$。黑线为 $400$ 次独立 Gillespie 实现的均值，灰带为 $95\%$ 置信区间，最宽处约与线宽相当；蓝虚线为并发型 EBCM \eqref{eq:master}，橙点线为均场闭合 \eqref{eq:mf}。(a) 易感比例；(b) 感染比例。](figures/fig1_trajectory_zh.png){width=6.5}

**定量校验。** 图 2 汇总三项定量比较。

![图 2　闭合的定量校验。(a) $t\le12$ 上 $S(t)$ 的平均绝对偏差随系统规模的变化；初始感染比例 $\varepsilon=0.01$，除 $N$ 外生成参数、仿真次数与 $\lambda/\lambda_c$ 全部固定，每点为 $4$ 个独立结构的均值，误差棒为结构间标准误，灰点划线为加权幂律拟合，指数 $-1.2$（子集拟合散布约 $\pm0.2$，与 $-1$ 相容）。(b) 终态规模随 $\lambda/\lambda_c$ 的变化，$N=6000$，初始感染比例 $\varepsilon=0.005$，每点 $400$ 次实现，误差棒为 $95\%$ 置信区间。(c) 爆发阈值随入出度相关的变化，$N=3000$，每点 $5$ 个独立双度序列，两个方向的误差棒均为标准误，实线为 \eqref{eq:lcrio}。](figures/fig2_validation_zh.png){width=6.5}

图 2(a) 分离出两种误差的不同性质。EBCM 的偏差随 $N$ 幂律衰减，四点加权拟合给出指数 $-1.22\pm0.05$；但这一误差棒只是拟合的统计误差，改用子集重拟合得 $-1.32\pm0.11$、$-1.21\pm0.10$、$-1.03\pm0.31$，因而实际不确定度约 $\pm0.2$，数据尚不能把该指数与 $-1$ 区分开。$-1$ 正是树状闭合所预期的：位形模型中短回路的密度为 $O(1/N)$。可以确定的是偏差趋于零。均场的偏差同样随 $N$ 下降（$0.220\to0.156$），却按 $a+b/N$ 外推到非零平台 $a=0.146\pm0.001$（四种子集拟合的散布）。两者因而不是同一类误差：EBCM 的可以靠增大系统消除，均场的不能，因为那是闭合层级的误差。$N=4000$ 时二者相差 $28$ 倍。

图 2(b) 检验 \eqref{eq:gain} 所预言的形状——均场的误差集中在阈值邻域。均场把阈值压低到 $\lambda_c^{\rm MF}=0.75\lambda_c$，因而在真实的亚临界区已给出宏观爆发：$\lambda=0.9\lambda_c$ 处高估 $5.1$ 倍、$\lambda=\lambda_c$ 处高估 $4.3$ 倍，深超临界区才收敛到 $1.2$ 倍。第 I 节所称"边基闭合相对均场的增益"由此定量化：增益不是均匀的百分数，而集中在最需要理论的区域。

同一张图也暴露出 EBCM 自身的限度，须如实报告。偏差在深超临界区最小（$\lambda\ge1.6\lambda_c$ 处 $\le0.7\%$），亚临界侧为 $+3.3\%$（$0.6\lambda_c$）与 $+1.7\%$（$0.8\lambda_c$），而在 $\lambda=\lambda_c$ 与 $1.1\lambda_c$ 处 EBCM 系统性高估 $18\%$ 与 $19\%$（分别为 $6.9\sigma$ 与 $8.4\sigma$）。偏差在亚临界区与阈值邻域一律为正，并在阈值邻域出现尖峰，深超临界区则收敛到零。这一高估是有限尺寸效应而非闭合误差，已直接验证：在 $\lambda=\lambda_c$ 处固定 $\varepsilon=0.005$、每点取 $3$ 个独立结构逐一增大系统，高估由 $N=1000$ 的 $+54\%$ 单调降至 $N=16000$ 的 $+4.4\%$——$N$ 增大 $16$ 倍，高估缩小约 $12$ 倍。这与 \eqref{eq:S}–\eqref{eq:master} 是 $N\to\infty$ 理论一致（阈值邻域关联长度发散，有限系统的偏离在此最大），也与均场那个不随 $N$ 消失的平台形成对照。即便如此，生产运行仍须在阈值邻域做有限尺寸标度外推，不能直接引用有限 $N$ 的 EBCM 值。

图 2(c) 给出阈值对结构的依赖。保持两条边际度序列不变、只重新配置 $k_{\rm in}$ 与 $k_{\rm out}$ 在节点上的配对，$\lambda_c$ 在 $r_{io}\in[-0.89,+0.89]$ 上从 $0.803$ 降到 $0.209$，变动近四倍，与 \eqref{eq:lcrio} 逐点吻合。需要说明，图 2(c) 的纵轴由 \eqref{eq:lc} 经 $\kappa$ 算得，因而它检验的是"$r_{io}$ 能把阈值推到多远"以及改写 \eqref{eq:lcrio} 的恒等性，而非对阈值的独立测量；对 $\lambda_c(r_{io})$ 的仿真检验属于第 VI 节。同时须记住 II.D 的宇称结论——$r_{io}$ 移动阈值，却不产生方向对称性破缺，二者并不矛盾。

**独立复算。** 上述数字均由与原推导不共享路径的方式复核一遍（脚本 `tools/section3_audit.py`）。$\tau=1$ 时方程组可解析积出终态：由 $\dot\Phi=-\beta x_{010}$ 与 $\dot x_{001}=\mu x_{010}$ 得 $x_{001}=\lambda^{-1}(1-\Phi)$，令 $t\to\infty$ 即

$$
\Phi_\infty=\frac{1+\lambda\,(1-\varepsilon)\,\psi_{\rm tail}(\Phi_\infty)}{1+\lambda}\label{eq:finalsize}
$$

这是一个标量不动点，求解不涉及积分器；三组参数下它与 RK4 积分的 $R(\infty)$ 相差不超过 $8\times10^{-14}$。$R_0=\tau\kappa\,\beta/(\beta+\mu)$ 另由直接抽样复核：按头槽抽取节点（即按 $k_{\rm in}$ 偏置，与新感染节点的到达方式一致），赋以 $\mathrm{Exp}(\mu)$ 感染期，对其每条出超边的每个头成员各跑一个 $\mathrm{Exp}(\beta)$ 时钟并计数子代；$\lambda=0.3,0.8,2.0$ 时测得值与预言的相对差分别为 $4.0\times10^{-3}$、$1.3\times10^{-3}$、$5.4\times10^{-4}$，均在统计误差以内。该检验既不使用方程组也不使用 \eqref{eq:finalsize}。均场阈值 \eqref{eq:lcmf} 同样以其自身终态二分复核，$t_{\max}=6400$ 时相对偏差 $1.3\times10^{-4}$。

阈值本身也已由仿真独立测定，此前它只是理论值。亚临界区单个种子引发的簇平均总规模为 $1/(1-R_0)$，故 $\varepsilon/R(\infty)\to1-R_0(\lambda)$，该量在 $\lambda_c$ 处归零。在 $N=8000$ 上取 $\lambda/\lambda_c\in[0.4,0.8]$ 五点、初始感染比例 $\varepsilon=0.002$、每点 $1500$ 次实现，加权外推得 $\lambda_c=0.3410\pm0.0058$，与 \eqref{eq:lc} 给出的 $0.3327$ 相差 $1.4\sigma$。拟合斜率给出 $\tau\kappa-1=2.85\pm0.05$，较公式值 $3.01$ 低约 $5\%$；这一偏低方向是预期的——纯分支形式忽略易感耗尽与不同种子所生簇的重叠，二者都压低 $R(\infty)$，从而抬高 $\varepsilon/R(\infty)$。

另有四项针对数值实现而非理论的核查，结果一并记录：位形模型的修复步在全部五种配置下未丢弃任何超边，实现的双度序列与请求值逐点相同；RK4 在 $\mathrm{d}t$ 由 $0.02$ 减至 $0.0025$ 时 $R(14)$ 的变化小于 $10^{-12}$，故积分步长不是误差来源；仿真按 Bernoulli$(\varepsilon)$ 播种与按固定种子数播种的终态相差 $0.5\sigma$，故 ODE 假设的"恰好 $\varepsilon$"不引入系统偏置；增量式 Gillespie 采样器与逐事件重算的朴素实现相差 $1.1\sigma$。
