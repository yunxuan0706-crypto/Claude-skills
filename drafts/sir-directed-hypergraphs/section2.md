# II. 模型与结构量

## A. 有向超图与记号

有向超图 $\mathcal{H}=(V,E)$ 由节点集 $V$（$N=|V|$）与超边集 $E$（$M=|E|$）构成。$E$ 为集合，即不含重复超边；第 V 节的生成器在配对后显式修复重复超边与 $T\cap H$ 相交的违规。每条超边是一个有序对

$$
e=\bigl(T(e),\,H(e)\bigr),\qquad T(e),H(e)\subseteq V,\qquad T(e)\cap H(e)=\varnothing
$$

其中尾集 $T(e)$ 为施加影响的一方，头集 $H(e)$ 为接受影响的一方。二者不相交，这既排除了自作用，也是后文重叠退化关系 \eqref{eq:reduction} 成立的前提。

主线取同质阶数，即对所有 $e$ 有 $|T(e)|=\tau$、$|H(e)|=\eta$；异质阶数的推广见附录。$\tau=\eta=1$ 时 $\mathcal{H}$ 退化为普通有向图；取 $T(e)\cup H(e)$ 为超边的支撑集则退化为无向超图。

节点的出度与入度分别计其作为尾成员与头成员的次数，

$$
k_{\rm out}(v)=\bigl|\{e\in E: v\in T(e)\}\bigr|,\qquad
k_{\rm in}(v)=\bigl|\{e\in E: v\in H(e)\}\bigr|
$$

二者的联合分布记为 $P(k_{\rm in},k_{\rm out})$。逐边计数给出两条握手关系

$$
\sum_{v\in V}k_{\rm out}(v)=\tau M,\qquad \sum_{v\in V}k_{\rm in}(v)=\eta M\label{eq:handshake}
$$

因此 $\langle k_{\rm out}\rangle/\langle k_{\rm in}\rangle=\tau/\eta$：尾集与头集基数不等时，两个平均度必然不等。

本文所用的随机系综为有向超图位形模型：给定双度序列 $d=\{(k_{\rm in}(v),k_{\rm out}(v))\}_{v\in V}$ 与阶数 $(\tau,\eta)$，对全部尾槽与头槽做均匀随机配对，再修复违规超边；具体实现见第 V 节。

本文的第一个结构参数是入出度相关

$$
r_{io}=\frac{\langle k_{\rm in}k_{\rm out}\rangle-\langle k_{\rm in}\rangle\langle k_{\rm out}\rangle}{\sigma_{k_{\rm in}}\sigma_{k_{\rm out}}}\label{eq:rio}
$$

即 $(k_{\rm in},k_{\rm out})$ 在节点上的 Pearson 相关系数——刻画同一节点两个方向的度值关系，而非相连节点之间的度值相似性 [1]。

必须强调 $r_{io}$ 是**双度序列 $d$ 的函数**，与超边如何连接无关。因此它由系综的输入决定，一经给定即不再随超边的重排而改变：调节 $r_{io}$ 只能在构造 $d$ 的阶段进行（在两条边际度序列固定的前提下重新配置 $k_{\rm in}$ 与 $k_{\rm out}$ 在节点上的配对），而下一节的重叠 $\alpha$ 则相反，它在 $d$ 固定后由超边重排调节。两个结构参数因而作用在生成流程的不同环节，互不干扰——这一分离正是第 V 节两套工具的分工依据。

## B. SIR 动力学

每个节点在任一时刻处于易感（S）、感染（I）、康复（R）三态之一。感染只沿尾→头方向传递：超边尾集中的感染成员驱动其头集中的易感节点，反向不发生传递。

记 $I(t)\subseteq V$ 为 $t$ 时刻的感染节点集，并令

$$
n_e(t)=\bigl|T(e)\cap I(t)\bigr|\label{eq:ne}
$$

为超边 $e$ 尾集中已感染的成员数。超边 $e$ 对其头集中每个尚处易感态的节点 $u\in H(e)$ 施加的感染率取

$$
\lambda_e(t)=\beta\,g\bigl(n_e(t)\bigr),\qquad g(0)=0
$$

其中 $\beta>0$ 为传播强度，$g$ 为感染核。本文主线取群体阈值核

$$
g(n)=\Theta(n-\theta)=\begin{cases}1,& n\ge\theta\\ 0,& n<\theta\end{cases}
$$

尾集感染成员数 $n_e\ge\theta$ 时超边处于激活态，以恒定率向头集传递；成员康复使 $n_e$ 回落到 $\theta$ 以下时超边失活。$\theta=1$ 对应简单传播，任一尾成员感染即足以激活超边；$\theta\ge2$ 则要求多个尾成员的联合状态，属于阈值型复杂传播 [2, 3]。两者的渗流对应物不同，将在第 IV 节分别处理。

阈值放在尾侧、且逐超边判定，是本文主线的机制选择。取线性核 $g(n)=n$ 可回收"每个感染尾成员独立传递"的情形；把阈值改放头侧、令头节点在跨超边累积的剂量超过阈值时被感染，则得到另一类机制，由附录的剂量核统一处理。

易感节点 $u$ 的总感染率是它所在全部头集的贡献之和，

$$
\Lambda_u(t)=\sum_{e:\,u\in H(e)}\lambda_e(t)\label{eq:rate}
$$

感染节点以恒定率 $\mu$ 独立康复，$I\to R$，康复后不再参与传递。以 $\mu$ 归一时间尺度后，过程只依赖比值

$$
\lambda\equiv\beta/\mu
$$

故给定感染核后，$\lambda$ 是唯一的连续动力学控制参数，爆发阈值指其临界值 $\lambda_c$；后文所称"阈值随结构参数变化"即指 $\lambda_c$ 对 $r_{io}$、$\alpha$ 等量的依赖。$\{S,I,R\}^{V}$ 上的演化因而是连续时间 Markov 过程，可用 Gillespie 算法精确抽样 [4]——这是第 III 节数值校验的基础。非指数康复期分布的推广见附录。

初始条件需随 $\theta$ 而变。$\theta=1$ 时取单个均匀随机选取的种子节点感染即可。$\theta\ge2$ 时单个种子无法激活任何超边——任一超边的感染尾成员数至多为 $1<\theta$——过程立即停在初态，爆发规模恒为 $1/N$；因此必须取一个规模不小于 $\theta$ 的种子集，我们取均匀随机的初始感染比例 $f$。这与 bootstrap 渗流以初始激活比例为控制参数的做法一致 [5]，$f$ 因而是 $\theta\ge2$ 分支的第二个控制参数。

有限系统几乎必然在有限时间内到达无感染者的吸收态。记该吸收态中康复节点的占比为爆发规模 $S$，记 $S[\mathcal{H}]$ 为给定 $\mathcal{H}$ 时 $S$ 对动力学随机性与种子选取的条件期望。两种初始条件下种子分布都均匀，因而与超边方向无关，这一点 II.D 将用到。

式 \eqref{eq:ne}–\eqref{eq:rate} 中方向性的作用是不对称的：节点的 $k_{\rm out}$ 决定它能把感染推给多少个头集，$k_{\rm in}$ 决定它暴露于多少个尾集。这一不对称正是 $r_{io}$ 得以进入动力学的通道。

## C. 有向超边重叠

第二个结构参数刻画超边之间如何共享节点。在无向超图中，两条超边的重叠只有"共享几个节点"一个自由度 [6, 7]；引入尾/头之分后，共享节点在两条超边中各自扮演的角色也成为信息。

对任意有序超边对 $(e,f)$，$e\neq f$，定义重叠列联矩阵的四个分量

$$
N_{ab}(e,f)=\bigl|A(e)\cap B(f)\bigr|,\qquad a,b\in\{T,H\}
$$

其中 $A=T$ 当 $a=T$、$A=H$ 当 $a=H$，$B$ 同理。四个通道各有独立的动力学含义：$N_{TT}$ 为**共发**，共享节点在两条超边中同为尾成员，同一感染源向两个不同头集扇出；$N_{HH}$ 为**共收**，共享节点同为头成员，两条超边的感染率按 \eqref{eq:rate} 在该节点上叠加（附录的头侧剂量核下，跨超边剂量亦在此累积）；$N_{HT}$ 为**串联**，$e$ 的头集落入 $f$ 的尾集，构成 $e\to f$ 的传播链；$N_{TH}$ 为反向串联。由 $T(e)\cap H(e)=\varnothing$，令 $\bar e=T(e)\cup H(e)$ 即得

$$
\bigl|\bar e\cap\bar f\bigr|=N_{TT}+N_{TH}+N_{HT}+N_{HH}\label{eq:reduction}
$$

也就是说，对任意一对超边，四通道之和精确等于其无向重叠：本节的重叠量是现行定义的方向细化，因而后文结果可与无向文献直接对照。

归一化到 $[0,1]$ 取（同质阶数下分母逐通道为常数：$TT$ 取 $\tau$、$HH$ 取 $\eta$、两个串联通道取 $\min(\tau,\eta)$）

$$
n_{ab}(e,f)=\frac{N_{ab}(e,f)}{\min\bigl(|A(e)|,|B(f)|\bigr)}\label{eq:norm}
$$

聚合方式受一条硬约束限制。对固定的 $e$ 逐节点计数可得四条求和规则，例如

$$
\sum_{f\neq e}N_{TT}(e,f)=\sum_{v\in T(e)}\bigl(k_{\rm out}(v)-1\bigr),\qquad
\sum_{f\neq e}N_{HT}(e,f)=\sum_{v\in H(e)}k_{\rm out}(v)
$$

其余两式同理。右端只依赖度序列，因此**在保度约束下四个通道的重叠总量是常数**；同质阶数下 \eqref{eq:norm} 的分母逐通道为常数，故无条件平均 $\langle n_{ab}\rangle$ 同样被度序列锁死，无法作为可调参数。可变的是重叠的堆积方式——是让少数超边对深度共享，还是让共享摊薄到大量超边对上。据此定义条件平均

$$
\alpha^{ab}_{\sigma\sigma'}=\Bigl\langle\, n_{ab}(e,f)\ \Bigm|\ N_{ab}(e,f)\ge1,\ \sigma(e)=\sigma,\ \sigma(f)=\sigma' \,\Bigr\rangle\label{eq:alpha}
$$

即在确实发生重叠的超边对上取平均，其中 $\sigma(e)=(|T(e)|,|H(e)|)$ 为超边的阶；该通道无重叠时 $\alpha^{ab}$ 不定义。记 $m_{ab}$ 为该通道中发生重叠的超边对数，则 $\alpha^{ab}=C_{ab}/m_{ab}$，常数 $C_{ab}$ 仅由度序列决定。因此调节 $\alpha$ 等价于调节 $m_{ab}$，这决定了第 V 节重连算法的实现方式。

式 \eqref{eq:alpha} 中阶指标取对角块（$\sigma=\sigma'$）即阶内重叠、取非对角块即阶间重叠，与无向情形下超边重叠矩阵的对角/非对角划分一致 [7]；本文只是在其每个矩阵元上再展开出四个方向通道。同质阶数下阶指标退化，以下略去并简记为 $\alpha^{ab}$；第 VII 节的真实数据阶数异质，届时阶指标重新生效。

互惠方向不由上述张量刻画。两条超边是否互相指回，取决于一组超边能否共同完成反向覆盖，这是超边级而非超边对级的判据。本文直接采用已有的强互惠定义 [8]：若存在一组超边，其尾集之并覆盖 $H(e)$ 且头集之并覆盖 $T(e)$，则称 $e$ 被强互惠；取

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

两个串联通道并不独立。对有序超边对恒有 $N_{TH}(e,f)=N_{HT}(f,e)$，故同质阶数下 $\alpha^{HT}$ 与 $\alpha^{TH}$ 逐实现相等（异质阶数下互为阶指标的转置），串联方向只携带一个分量，记 $\alpha^{\rightrightarrows}=\alpha^{HT}$。四个通道在聚合层面因而只有三个独立分量。

按宇称重组：$\alpha^{\rightrightarrows}$ 与并联重叠 $\alpha^{\parallel}=\tfrac12(\alpha^{TT}+\alpha^{HH})$ 在 $\mathcal{R}$ 下为偶，而极性失衡

$$
\Delta\alpha=\tfrac12\bigl(\alpha^{TT}-\alpha^{HH}\bigr)
$$

为奇。$\alpha^{\parallel}$ 必须在扫描 $\alpha^{\rightrightarrows}$ 或 $\Delta\alpha$ 时固定，否则三种效应混在同一条扫描线上。

入出度相关 $r_{io}$ 在 $\mathcal{R}$ 下为**偶**：翻转互换每个节点的 $(k_{\rm in},k_{\rm out})$，而 \eqref{eq:rio} 的 Pearson 相关系数对其两个变量对称，故 $r_{io}(\mathcal{R}\mathcal{H})=r_{io}(\mathcal{H})$。

据此可把结构量按宇称分类。本文所用结构量中，$\mathcal{R}$-奇的来源有两个且相互独立：双度序列 $d$ 关于分量互换的不对称性，以及极性失衡 $\Delta\alpha$。尾头基数之差 $\tau-\eta$ 不是第三个独立来源——由 \eqref{eq:handshake}，$d$ 互换对称蕴含 $\tau=\eta$，故 $\tau\neq\eta$ 只是序列不对称的一个充分标志；反之序列不对称时 $\tau=\eta$ 仍可成立。

这一分类对终态有直接约束。设结构系综 $\mathbb{P}[\mathcal{H}]$ 生成随机有向超图，$S[\mathcal{H}]$ 为 II.B 定义的条件期望爆发规模；由于种子均匀选取，$S$ 的定义本身不引入方向偏好。

> **命题 1.** 若 $\mathbb{P}$ 在 $\mathcal{R}$ 下不变，即 $\mathbb{P}[\mathcal{R}\mathcal{H}]=\mathbb{P}[\mathcal{H}]$ 对所有 $\mathcal{H}$ 成立，则 $\langle S\rangle_{\mathbb{P}}=\langle S\circ\mathcal{R}\rangle_{\mathbb{P}}$。

证明只需一次换元：$\mathcal{R}$ 是结构空间上的可测对合，故 $\sum_{\mathcal{H}}\mathbb{P}[\mathcal{H}]\,S[\mathcal{R}\mathcal{H}]=\sum_{\mathcal{H}}\mathbb{P}[\mathcal{R}\mathcal{H}]\,S[\mathcal{H}]=\sum_{\mathcal{H}}\mathbb{P}[\mathcal{H}]\,S[\mathcal{H}]$。

命题 1 的逆否形式给出方向对称性破缺的必要条件：**系综层面的破缺要求生成系综本身破坏 $\mathcal{R}$-不变性**。在本文所用结构量中，能承载这一破坏的只有上述两个 $\mathcal{R}$-奇的来源。

但必要条件被满足并不意味着破缺可以在 $\langle S\rangle$ 上被看到。事实上存在一条更强的不变性，它与系综是否 $\mathcal{R}$-不变无关：

> **命题 2.** 设终态可表示为某个随机结构上的可达性——即从种子 $v$ 出发的最终感染集为 $v$ 的出分支 $\mathrm{out}(v)$，如 [9] 的半有向流行渗流网络所给出。则对**任意**有向超图 $\mathcal{H}$，均匀单种子下的平均爆发规模在方向翻转下不变。

证明同样只需一次计数：$\sum_v|\mathrm{out}(v)|=\#\{(v,u):u\ \text{自}\ v\ \text{可达}\}=\sum_v|\mathrm{in}(v)|$，而 $\mathcal{R}$ 把可达关系转置，可达有序对的**总数**因而逐实现不变；除以 $N^2$ 即得 $\langle S\rangle$ 不变。

命题 2 的后果是：$\langle S\rangle$ **不能作为方向对称性破缺的序参量**。它在任何结构上都对称，与 $\Delta\alpha$ 是否为零无关。第 VI 节的数值确认了这一点：即便在 $\tau=3,\eta=1$ 这种最不对称的构型上，$\langle S\rangle$ 在翻转前后的差异也始终落在统计涨落内。

破缺出现在可达对**如何分布**，而非其总数。按 [9]，大爆发的**概率**由巨入分支的相对规模决定，而爆发的**终态规模**由巨出分支决定；$\mathcal{R}$ 恰好交换这两者。因此本文取

$$
\Pi=\Pr\bigl(S>S_{\rm cut}\bigr),\qquad \bar S=\mathbb{E}\bigl[S\mid S>S_{\rm cut}\bigr]
$$

为破缺的序参量，其中 $S_{\rm cut}$ 为区分自限爆发与大爆发的阈值。命题 2 同时给出一条可证伪的预言：$\langle S\rangle\simeq\Pi\bar S$ 在翻转下守恒，故 $\Pi$ 与 $\bar S$ **若破缺必反向、相互补偿**。第 VI 节将检验这一补偿关系。

在此基础上，$r_{io}$ 的角色是一个可直接检验的推论：它不能单独产生破缺。设双度序列 $d$ 关于分量互换对称，即 $\mathcal{R}d$ 与 $d$ 至多相差一个节点重标号。$\mathcal{R}$ 把位形模型系综 $\mathrm{CM}(d)$ 映为 $\mathrm{CM}(\mathcal{R}d)$，后者是前者在该重标号下的像；而 $S$ 与节点标号无关，故命题 1 的结论成立，系综平均的爆发规模在方向翻转下相同。关键在于：这类对称序列可以取遍 $r_{io}\in[-1,1]$ 的整个范围（例如 $d$ 由若干对 $(a,b)$ 与 $(b,a)$ 组成时序列对称而 $r_{io}$ 可为任意值），因此**沿 $r_{io}$ 扫描的整条线上都不应出现破缺，尽管 $r_{io}$ 确实移动阈值**。

这决定了第 VI 节破缺实验的设计：全部测量以 $(\Pi,\bar S)$ 为序参量而非 $\langle S\rangle$；在互换对称的双度序列上做 $r_{io}$ 扫描作为阴性对照；破缺实验保持 $d$ 对称、只用保度重排调节 $\Delta\alpha$，并同时报告 $\langle S\rangle$ 以验证它按命题 2 保持不变——后者是整套测量的内部一致性检查。

命题 1 是关于系综平均的陈述。对单个实现 $\mathcal{H}$，$S[\mathcal{H}]$ 与 $S[\mathcal{R}\mathcal{H}]$ 一般不等，其差为有限尺寸涨落，随 $N$ 增大而趋零。

---

## 参考文献

1. M. E. J. Newman, *Physical Review Letters*, **89**, 208701 (2002). doi:10.1103/PhysRevLett.89.208701
2. D. J. Watts, *Proceedings of the National Academy of Sciences*, **99**, 5766–5771 (2002). doi:10.1073/pnas.082090499
3. D. Centola, M. Macy, *American Journal of Sociology*, **113**, 702–734 (2007). doi:10.1086/521848
4. D. T. Gillespie, *The Journal of Physical Chemistry*, **81**, 2340–2361 (1977). doi:10.1021/j100540a008
5. G. J. Baxter, S. N. Dorogovtsev, A. V. Goltsev, J. F. F. Mendes, *Physical Review E*, **82**, 011103 (2010). doi:10.1103/PhysRevE.82.011103
6. F. Malizia, S. Lamata-Otín, M. Frasca, V. Latora, J. Gómez-Gardeñes, *Nature Communications*, **16** (2025). doi:10.1038/s41467-024-55506-1
7. S. Lamata-Otín, F. Malizia, V. Latora, M. Frasca, J. Gómez-Gardeñes, *Physical Review E*, **111**, 034302 (2025). doi:10.1103/PhysRevE.111.034302
8. Q. F. Lotito, A. Vendramini, A. Montresor, F. Battiston, *Communications Physics*, **9** (2026). doi:10.1038/s42005-025-02472-9
9. E. Kenah, J. M. Robins, *Physical Review E*, **76**, 036113 (2007). doi:10.1103/PhysRevE.76.036113
