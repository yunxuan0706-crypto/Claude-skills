# II. 模型与结构量

## A. 有向超图与记号

有向超图 $\mathcal{H}=(V,E)$ 由节点集 $V$（$N=|V|$）与超边集 $E$（$M=|E|$）构成。每条超边是一个有序对

$$
e=\bigl(T(e),\,H(e)\bigr),\qquad T(e),H(e)\subseteq V,\qquad T(e)\cap H(e)=\varnothing
$$

其中尾集 $T(e)$ 为施加影响的一方，头集 $H(e)$ 为接受影响的一方。尾集与头集不相交，这既排除了自作用，也是后文重叠退化关系 \eqref{eq:reduction} 成立的前提。主线取同质阶数，即对所有 $e$ 有 $|T(e)|=\tau$、$|H(e)|=\eta$；异质阶数的推广见附录。$\tau=\eta=1$ 时 $\mathcal{H}$ 退化为普通有向图；忽略尾/头之分则退化为无向超图。

节点的出度与入度分别计其作为尾成员与头成员的次数，

$$
k_{\rm out}(v)=\bigl|\{e\in E: v\in T(e)\}\bigr|,\qquad
k_{\rm in}(v)=\bigl|\{e\in E: v\in H(e)\}\bigr|
$$

二者的联合分布记为 $P(k_{\rm in},k_{\rm out})$。逐边计数给出两条握手关系

$$
\sum_{v\in V}k_{\rm out}(v)=\tau M,\qquad \sum_{v\in V}k_{\rm in}(v)=\eta M
$$

因此 $\langle k_{\rm out}\rangle/\langle k_{\rm in}\rangle=\tau/\eta$：尾集与头集基数不等时，两个平均度必然不等。本文的第一个结构参数是入出度相关

$$
r_{io}=\frac{\langle k_{\rm in}k_{\rm out}\rangle-\langle k_{\rm in}\rangle\langle k_{\rm out}\rangle}{\sigma_{k_{\rm in}}\sigma_{k_{\rm out}}}\label{eq:rio}
$$

即 $(k_{\rm in},k_{\rm out})$ 在节点上的 Pearson 相关系数。它刻画的是同一节点两个方向的度值关系，与相连节点之间的度—度同配性 [1] 是不同的量。

## B. SIR 动力学

每个节点在任一时刻处于易感（S）、感染（I）、康复（R）三态之一。感染只沿尾→头方向传递：一条超边的尾集成员被感染后，对其头集中的易感节点施加压力，反向不发生传递。

记 $I(t)\subseteq V$ 为 $t$ 时刻的感染节点集，并令

$$
n_e(t)=\bigl|T(e)\cap I(t)\bigr|\label{eq:ne}
$$

为超边 $e$ 尾集中已感染的成员数。超边 $e$ 对其头集中每个易感节点 $u\in H(e)\cap S(t)$ 施加的感染率取

$$
\lambda_e(t)=\beta\,g\bigl(n_e(t)\bigr),\qquad g(0)=0
$$

其中 $\beta>0$ 为传播强度，$g$ 为感染核。本文主线取群体阈值核

$$
g(n)=\Theta(n-\theta)=\begin{cases}1,& n\ge\theta\\ 0,& n<\theta\end{cases}
$$

即尾集中感染成员数达到 $\theta$ 时超边被激活，此后以恒定率向头集传递。$\theta=1$ 对应简单传播：任一尾成员被感染即足以激活超边。$\theta\ge2$ 则要求多个尾成员的联合状态，属于阈值型复杂传播 [2, 3]，两者的渗流对应物不同，将在第 IV 节分别处理。取线性核 $g(n)=n$ 可回收"每个感染尾成员独立传递"的情形；含时间累积的剂量核见附录。

易感节点 $u$ 的总感染率是它所在全部头集的贡献之和，

$$
\Lambda_u(t)=\sum_{e:\,u\in H(e)}\lambda_e(t)\label{eq:rate}
$$

感染节点以恒定率 $\mu$ 独立康复，$I\to R$，康复后不再参与传递。$\{S,I,R\}^{V}$ 上的状态演化因而是一个连续时间 Markov 过程，可用 Gillespie 算法精确抽样 [4]，这是第 III 节数值校验的基础。非指数康复期分布的推广见附录。

式 \eqref{eq:ne}–\eqref{eq:rate} 中方向性的作用是不对称的：节点的 $k_{\rm out}$ 决定它能把感染推给多少个头集，$k_{\rm in}$ 决定它暴露于多少个尾集。这一不对称正是 $r_{io}$ 得以进入动力学的通道。

## C. 有向超边重叠

第二个结构参数刻画超边之间如何共享节点。在无向超图中，两条超边的重叠只有"共享几个节点"一个自由度 [5, 6]；引入尾/头之分后，共享节点在两条超边中各自扮演的角色也成为信息。

对任意有序超边对 $(e,f)$，$e\neq f$，定义重叠列联矩阵的四个分量

$$
N_{ab}(e,f)=\bigl|A(e)\cap B(f)\bigr|,\qquad a,b\in\{T,H\}
$$

其中 $A=T$ 当 $a=T$、$A=H$ 当 $a=H$，$B$ 同理。四个通道各有独立的动力学含义：$N_{TT}$ 为**共发**，共享节点在两条超边中同为尾成员，同一感染源向两个不同头集扇出；$N_{HH}$ 为**共收**，共享节点同为头成员，来自两个不同尾集的压力在同一节点上汇聚，群体阈值 $\theta$ 正是在此累积跨超边的剂量；$N_{HT}$ 为**串联**，$e$ 的头集落入 $f$ 的尾集，构成 $e\to f$ 的传播链；$N_{TH}$ 为反向串联。由 $T(e)\cap H(e)=\varnothing$，令 $\bar e=T(e)\cup H(e)$ 即得

$$
\bigl|\bar e\cap\bar f\bigr|=N_{TT}+N_{TH}+N_{HT}+N_{HH}\label{eq:reduction}
$$

也就是说，对任意一对超边，四通道之和精确等于其无向重叠。本文的重叠量因而是现行定义的方向细化，而非另一套互不相容的定义，这使得后文结果可与无向文献直接对照。

归一化到 $[0,1]$ 取

$$
n_{ab}(e,f)=\frac{N_{ab}(e,f)}{\min\bigl(|A(e)|,|B(f)|\bigr)}
$$

聚合方式需要谨慎。对固定的 $e$ 逐节点计数可得四条求和规则，例如

$$
\sum_{f\neq e}N_{TT}(e,f)=\sum_{v\in T(e)}\bigl(k_{\rm out}(v)-1\bigr),\qquad
\sum_{f\neq e}N_{HT}(e,f)=\sum_{v\in H(e)}k_{\rm out}(v)
$$

其余两式同理。右端只依赖度序列，因此**在保度约束下四个通道的重叠总量是常数**：无条件平均 $\langle n_{ab}\rangle$ 正比于该常数除以超边对数，恒定不变，无法作为可调参数。可变的是重叠的堆积方式——是让少数超边对深度共享，还是让共享摊薄到大量超边对上。据此定义条件平均

$$
\alpha^{ab}_{\sigma\sigma'}=\Bigl\langle\, n_{ab}(e,f)\ \Bigm|\ N_{ab}(e,f)\ge1,\ \sigma(e)=\sigma,\ \sigma(f)=\sigma' \,\Bigr\rangle\label{eq:alpha}
$$

即在确实发生重叠的超边对上取平均，其中 $\sigma(e)=(|T(e)|,|H(e)|)$ 为超边的阶。记 $m_{ab}$ 为该通道中发生重叠的超边对数，则 $\alpha^{ab}\propto\text{常数}/m_{ab}$，调节 $\alpha$ 等价于调节 $m_{ab}$；这决定了第 V 节重连算法的实现方式。

式 \eqref{eq:alpha} 中阶指标取对角块（$\sigma=\sigma'$）即阶内重叠、取非对角块即阶间重叠，与无向情形下超边重叠矩阵的对角/非对角划分一致 [6]；本文只是在其每个矩阵元上再展开出四个方向通道。同质阶数下阶指标退化，以下略去并简记为 $\alpha^{ab}$。

互惠方向不由上述张量刻画。两条超边是否互相指回，取决于一组超边能否共同完成反向覆盖，这是超边级而非超边对级的判据。本文直接采用已有的强互惠定义 [7]：若存在一组超边，其尾集之并覆盖 $H(e)$ 且头集之并覆盖 $T(e)$，则称 $e$ 被强互惠；取

$$
\alpha^{\rightleftarrows}=\frac{\bigl|\{e\in E:\ e\ \text{被强互惠}\}\bigr|}{M}
$$

采用这一定义是为了与已有测量口径一致，便于第 VII 节与真实数据的报告值并列。需要注意它与 $\alpha^{ab}$ 不同类：退化关系 \eqref{eq:reduction} 对它不成立，它也不进入第 IV 节转移算子的权重。

## D. 方向翻转与宇称

定义全局方向翻转 $\mathcal{R}$ 为对所有超边交换尾集与头集，

$$
\mathcal{R}:\ e=\bigl(T(e),H(e)\bigr)\ \longmapsto\ \bigl(H(e),T(e)\bigr)
$$

$\mathcal{R}$ 是有向超图空间上的对合（$\mathcal{R}^2=\mathrm{id}$）。它把每个节点的入度与出度互换，把 $\tau$ 与 $\eta$ 互换，并在四个重叠通道上作用为

$$
\mathcal{R}:\quad \alpha^{TT}\leftrightarrow\alpha^{HH},\qquad \alpha^{HT}\leftrightarrow\alpha^{TH}
$$

据此把四通道重组为宇称本征量：串联重叠 $\alpha^{\rightrightarrows}=\tfrac12(\alpha^{HT}+\alpha^{TH})$ 与并联重叠 $\alpha^{\parallel}=\tfrac12(\alpha^{TT}+\alpha^{HH})$ 在 $\mathcal{R}$ 下为偶，而极性失衡

$$
\Delta\alpha=\tfrac12\bigl(\alpha^{TT}-\alpha^{HH}\bigr)
$$

为奇。$\alpha^{\parallel}$ 必须在扫描 $\alpha^{\rightrightarrows}$ 或 $\Delta\alpha$ 时固定，否则三种效应混在同一条扫描线上。

入出度相关 $r_{io}$ 在 $\mathcal{R}$ 下为**偶**：翻转互换每个节点的 $(k_{\rm in},k_{\rm out})$，而 \eqref{eq:rio} 的 Pearson 相关系数对其两个变量对称，故 $r_{io}(\mathcal{R}\mathcal{H})=r_{io}(\mathcal{H})$。据此可把结构量按宇称分类，$\mathcal{R}$-奇者共三类：尾头基数之差 $\tau-\eta$、联合度分布关于对角线的反对称部分 $P(k_{\rm in},k_{\rm out})-P(k_{\rm out},k_{\rm in})$，以及极性失衡 $\Delta\alpha$。

这一分类对终态有直接约束。设结构系综 $\mathbb{P}[\mathcal{H}]$ 生成随机有向超图，$S[\mathcal{H}]$ 为在 $\mathcal{H}$ 上按 II.B 演化所得的终态爆发规模。

> **命题 1.** 若 $\mathbb{P}$ 在 $\mathcal{R}$ 下不变，即 $\mathbb{P}[\mathcal{R}\mathcal{H}]=\mathbb{P}[\mathcal{H}]$ 对所有 $\mathcal{H}$ 成立，则 $\langle S\rangle_{\mathbb{P}}=\langle S\circ\mathcal{R}\rangle_{\mathbb{P}}$。

证明只需一次换元：$\mathcal{R}$ 是结构空间上的可测对合，故 $\sum_{\mathcal{H}}\mathbb{P}[\mathcal{H}]\,S[\mathcal{R}\mathcal{H}]=\sum_{\mathcal{H}}\mathbb{P}[\mathcal{R}\mathcal{H}]\,S[\mathcal{H}]=\sum_{\mathcal{H}}\mathbb{P}[\mathcal{H}]\,S[\mathcal{H}]$。

命题 1 的逆否形式给出方向对称性破缺的必要条件：**系综层面的破缺要求上述三类 $\mathcal{R}$-奇量至少有一个非零。** 由此立即可知，在 $\tau=\eta$ 且 $P$ 关于对角线对称的设定下单独扫描 $r_{io}$ 不可能产生破缺——尽管 $r_{io}$ 确实移动阈值。这一推论决定了第 VI 节破缺实验的设计：固定前两类 $\mathcal{R}$-奇量为零，只扫 $\Delta\alpha$，并以 $r_{io}$ 扫描作为无破缺的阴性对照。需要强调命题 1 是关于系综平均的陈述；对单个实现 $\mathcal{H}$，$S[\mathcal{H}]$ 与 $S[\mathcal{R}\mathcal{H}]$ 一般不等，其差为有限尺寸涨落。

---

## 参考文献

1. M. E. J. Newman, *Physical Review Letters*, **89**, 208701 (2002). doi:10.1103/PhysRevLett.89.208701
2. D. J. Watts, *Proceedings of the National Academy of Sciences*, **99**, 5766–5771 (2002). doi:10.1073/pnas.082090499
3. D. Centola, M. Macy, *American Journal of Sociology*, **113**, 702–734 (2007). doi:10.1086/521848
4. D. T. Gillespie, *The Journal of Physical Chemistry*, **81**, 2340–2361 (1977). doi:10.1021/j100540a008
5. F. Malizia, S. Lamata-Otín, M. Frasca, V. Latora, J. Gómez-Gardeñes, *Nature Communications*, **16** (2025). doi:10.1038/s41467-024-55506-1
6. S. Lamata-Otín, F. Malizia, V. Latora, M. Frasca, J. Gómez-Gardeñes, *Physical Review E*, **111**, 034302 (2025). doi:10.1103/PhysRevE.111.034302
7. Q. F. Lotito, A. Vendramini, A. Montresor, F. Battiston, *Communications Physics*, **9** (2026). doi:10.1038/s42005-025-02472-9
