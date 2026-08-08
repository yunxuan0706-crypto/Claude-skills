# I. 引言

在成对接触网络上，SIR 传播的理论已沿两条彼此独立的路径发展成熟。终态一侧由渗流承担：任意度分布随机图的连通结构可由生成函数方法解析给出 [1]，而键渗流映射把最终感染规模化归为占据概率 $\phi$ 下的巨簇规模 [2]；这一思路可上溯至一般流行过程与动力学渗流的等价性 [3]。由此立即可见度分布异质性的作用：临界占据概率为 $\phi_c=\langle k\rangle/(\langle k^2\rangle-\langle k\rangle)$，二阶矩发散时阈值趋于零——这一结论最早在无标度网络上的 SIS 模型中给出 [4]。Kenah 与 Robins 随后指出，SIR 与键渗流之间并非严格同构：与之精确同构的是一个半有向的"流行渗流网络"，其中爆发规模分布对应出分支、阈值对应巨强连通分支的出现、终态规模对应巨出分支 [5]。值得注意的是，方向性由此已经内在地进入了 SIR 的渗流表述，尽管接触网络本身是无向的。

含时一侧由方程法承担。Volz 用成对近似给出随机网络上 SIR 的低维闭合 [6]，Miller 等将其整理为边基房室模型（EBCM），以"沿一条边尚未传来感染"的概率为核心变量，用少数常微分方程同时容纳异质接触率与伙伴关系的有限持续时间 [7]。传染与康复过程的分布形式并非细节：偏离马尔可夫假设会显著改变流行阈值 [8]。两条路径分工明确——渗流给出阈值与终态，方程给出峰值与时间演化——并可相互印证 [9]。

方向性作为接触网络本身的属性也早已被处理，且入出度相关从一开始就是其中的关键量。有向随机图各巨分支的规模由入出度联合分布 $P(k_{\rm in},k_{\rm out})$ 决定；当该分布不可分解——即入度与出度相关时——巨强连通分支的规模会偏离巨入分支与巨出分支规模之积 [10]。在任意两点度相关与双向边同时存在的一般情形下，有向渗流的阈值与各巨分支规模已可解析求出 [11]，这套机制并已直接用于有向接触网络上的流行预测 [12]。入出度相关对动力学的影响亦有直接证据：在有向网络上的阈值型传播中，正的入出度相关在较宽的平均入度范围内提高系统稳健性，负相关则相反 [13]。它与无向网络中的度—度同配性 [14] 并非同一件事——后者刻画相连节点之间度值的相似性，前者刻画同一节点两个方向的度值关系。

与此同时，传播研究的重心转向高阶交互——许多传染并非沿边逐个发生，而是在一个群体内集体发生 [15, 16]，而这类高阶交互在真实数据中普遍存在 [17]。单纯复形上的高阶影响与强化机制可诱导不连续相变，并出现健康态与流行态共存的双稳区 [18]；超图上的社会传播已有适用于任意超图的解析框架，展现一阶与二阶相变、双稳与滞后 [19]；而把异质暴露与最小感染剂量结合，会导出一个普适的非线性感染核，并随之出现不连续相变与超指数增长 [20]——这与经典的阈值型复杂传播 [21, 22] 一脉相承。

结构对高阶传播的影响同样明确。链接度分布的异质性可抑制爆炸式相变 [23]；超核分解给出一种中心性，传播过程局域于中心超核 [24]；由超边之间共享节点数定义的高阶连通分支，其巨分支的存在是单一种子实现全局爆发的必要条件 [25]；哪些群体最适合作为种子也已被刻画 [26]。尤为关键的是，超边之间的重叠决定了转变是否爆炸式，且其作用方向值得注意：只有在阶内重叠较低时才出现爆炸性与双稳 [27]。重叠的组织方式可由一个超边重叠矩阵刻画，其对角元与非对角元分别对应阶内与阶间重叠 [28]，而阶间相关驱动系统沿不同路径走向爆炸式传播 [29]。超图上的渗流理论亦已建立 [30]，相关进展见综述 [31]。这些工作共享同一前提：超边是无向的，一条超边内的所有节点对称地相互施加影响。

然而现实中的高阶交互往往是有向的——施加影响的一方与接受影响的一方并不对称。有向超图正是为刻画这种不对称而设的形式化：每条超边分为施加影响的尾集与接受影响的头集（二者不相交），这一结构在算法与运筹领域已使用三十余年 [32]。代谢与化学反应把一组底物转化为一组产物，是其典型实例 [33]；比特币交易与引用数据亦已按有向超图分析 [34]。针对有向超图的结构刻画工具近年迅速成型：真实数据的微观组织、高阶互惠性与 motif 已被系统刻画 [34, 35]；无向超图中超边重叠的测度与生成器已给出 [36]，位形模型与保度随机化亦有成熟实现 [37, 38]；而保入出度序列的均匀采样与非退化零模型在有向超图上也已可用 [39, 40]。

动力学一侧则由 Li 等开启。他们在有向超图上建立社会传播模型，发现不连续相变的双稳区间随有向强度减弱而收缩，表明方向性确实改变高阶传播的相变结构 [41]。该工作采用均场闭合，不涉及渗流，其方向性由单一标量强度参数调节。均场忽略邻域相关，误差在阈值邻域最大；在成对网络上，正是 EBCM 补上了这一层 [6, 7]。因此，有向超图上的传播阈值究竟如何依赖于结构相关，目前既没有超出均场的含时理论，也没有可用于解析扫描的渗流理论。

本文沿方程与渗流两条路径推进这一问题，并证明二者在同一极限下相合。方程一侧，我们提出有向超图 SIR 的并发型 EBCM 方程；其技术困难在于一个头节点同时承受多条超边、多个尾成员的并发压力，闭合变量不能简单相乘。我们给出闭合方案，在尾集基数 $\tau=1,2$ 的可控情形下手推验证，并以精确 Gillespie 仿真 [42] 在阈值邻域校验，从而定量给出边基闭合相对均场的增益。渗流一侧，我们把成对网络上的 SIR–渗流映射及其半有向修正 [2, 5] 与有向渗流的生成函数机制 [10, 11] 推广到有向超图，得到爆发阈值随结构参数变化的半解析依赖。两条路径互不依赖，因而可互为验证；我们进一步证明二者的自洽方程在 $\mu\to0$ 极限下同构，给出同一阈值。

为使结构依赖可被受控地检验，我们给出有向超边重叠的可计算定义。对每一对有序超边构造 $2\times2$ 重叠列联矩阵，其四个通道分别对应共享节点在两条超边中同为尾成员、同为头成员、或一为头一为尾（两个方向各一个通道）。按阶聚合后得到一个张量，其阶对角块与非对角块分别对应已有的阶内与阶间重叠 [27, 28]；对任意一对超边，四个通道之和精确等于无向情形下的重叠，因此本文结果可与无向文献直接对照。互惠通道直接采用文献中的强互惠判据 [34]。在此定义下，我们在严格保持逐节点入出度的约束下，运行以入出度相关 $r_{io}$ 与重叠 $\alpha$ 为目标的定向重连，从而沿单一结构坐标做受控扫描；零模型基线取自已有的均匀采样器 [39, 40]。

这里先指出一个与直觉相反的结果。全局方向翻转 $\mathcal{R}$ 是结构空间上的对合；若生成系综在 $\mathcal{R}$ 下不变，则爆发规模的系综平均在方向翻转下必然不变。因此方向对称性破缺要求某个 $\mathcal{R}$-奇不变量非零。而入出度相关 $r_{io}$ 在 $\mathcal{R}$ 下是偶的——方向翻转互换每个节点的入度与出度，而相关系数对其两个变量对称。这意味着单独调节 $r_{io}$ 不可能产生系综层面的方向不对称，尽管它确实移动阈值。本文所用结构量中，$\mathcal{R}$-奇者有三个：尾集与头集基数之差、联合度分布关于对角线的反对称部分，以及共发重叠与共收重叠之差 $\Delta\alpha$。把前两个置零后，$\Delta\alpha$ 就是破缺的最低阶来源；我们据此设计破缺实验，并给出 $r_{io}$ 扫描下无破缺的阴性对照。

还需说明一点：当群体阈值 $\theta\ge2$ 时，激活规则不再是单条超边的独立事件，模型转而落入 bootstrap 与 $k$-core 型渗流，其相变结构不同于普通渗流——$k$-core 的出现是混合型相变 [43]，bootstrap 渗流则可能同时存在一个连续阈值与一个更高的混合型阈值，且后者在二阶矩发散时消失 [44]；高阶渗流与 $K$-core 渗流之间的联系在多重超图上已被建立 [45]。本文对 $\theta=1$ 与 $\theta\ge2$ 分别处理。

本文余下部分安排如下。第 II 节给出有向超图与 SIR 动力学的定义、有向超边重叠的可计算形式，以及结构量在方向翻转下的宇称分类。第 III 节推导并发型 EBCM 方程并给出仿真校验。第 IV 节建立 SIR 到有向超图渗流的映射与生成函数阈值理论。第 V 节介绍保度生成与定向重连工具。第 VI 节给出阈值对结构参数的依赖、方向对称性破缺，以及两方法的交叉验证。第 VII 节在真实有向超图数据上检验模型。第 VIII 节为讨论与结论。

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

即 $(k_{\rm in},k_{\rm out})$ 在节点上的 Pearson 相关系数。它刻画的是同一节点两个方向的度值关系，与相连节点之间的度—度同配性 [14] 是不同的量。

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

即当尾集感染成员数 $n_e\ge\theta$ 时超边处于激活态，以恒定率向头集传递；成员康复使 $n_e$ 回落到 $\theta$ 以下时超边失活。$\theta=1$ 对应简单传播：任一尾成员处于感染态即足以激活超边。$\theta\ge2$ 则要求多个尾成员的联合状态，属于阈值型复杂传播 [21, 22]，两者的渗流对应物不同，将在第 IV 节分别处理。取线性核 $g(n)=n$ 可回收"每个感染尾成员独立传递"的情形。另一类机制把阈值放在头侧——头节点在跨超边累积的剂量超过阈值时被感染；该变体由附录的剂量核统一处理，主线不采用。

易感节点 $u$ 的总感染率是它所在全部头集的贡献之和，

$$
\Lambda_u(t)=\sum_{e:\,u\in H(e)}\lambda_e(t)\label{eq:rate}
$$

感染节点以恒定率 $\mu$ 独立康复，$I\to R$，康复后不再参与传递。$\{S,I,R\}^{V}$ 上的状态演化因而是一个连续时间 Markov 过程，可用 Gillespie 算法精确抽样 [42]，这是第 III 节数值校验的基础。非指数康复期分布的推广见附录。初始条件取单个均匀随机选取的种子节点感染、其余易感；有限系统几乎必然在有限时间内到达无感染者的吸收态。记吸收态中康复节点的占比为爆发规模 $S$，并记 $S[\mathcal{H}]$ 为给定 $\mathcal{H}$ 时 $S$ 对动力学随机性与种子选取的条件期望；均匀种子的分布与超边方向无关，这一点将在 II.D 用到。

式 \eqref{eq:ne}–\eqref{eq:rate} 中方向性的作用是不对称的：节点的 $k_{\rm out}$ 决定它能把感染推给多少个头集，$k_{\rm in}$ 决定它暴露于多少个尾集。这一不对称正是 $r_{io}$ 得以进入动力学的通道。

## C. 有向超边重叠

第二个结构参数刻画超边之间如何共享节点。在无向超图中，两条超边的重叠只有"共享几个节点"一个自由度 [27, 28]；引入尾/头之分后，共享节点在两条超边中各自扮演的角色也成为信息。

对任意有序超边对 $(e,f)$，$e\neq f$，定义重叠列联矩阵的四个分量

$$
N_{ab}(e,f)=\bigl|A(e)\cap B(f)\bigr|,\qquad a,b\in\{T,H\}
$$

其中 $A=T$ 当 $a=T$、$A=H$ 当 $a=H$，$B$ 同理。四个通道各有独立的动力学含义：$N_{TT}$ 为**共发**，共享节点在两条超边中同为尾成员，同一感染源向两个不同头集扇出；$N_{HH}$ 为**共收**，共享节点同为头成员，两条超边的感染率按 \eqref{eq:rate} 在该节点上叠加（附录的头侧剂量核下，跨超边剂量亦在此累积）；$N_{HT}$ 为**串联**，$e$ 的头集落入 $f$ 的尾集，构成 $e\to f$ 的传播链；$N_{TH}$ 为反向串联。由 $T(e)\cap H(e)=\varnothing$，令 $\bar e=T(e)\cup H(e)$ 即得

$$
\bigl|\bar e\cap\bar f\bigr|=N_{TT}+N_{TH}+N_{HT}+N_{HH}\label{eq:reduction}
$$

也就是说，对任意一对超边，四通道之和精确等于其无向重叠。本文的重叠量因而是现行定义的方向细化，而非另一套互不相容的定义，这使得后文结果可与无向文献直接对照。

归一化到 $[0,1]$ 取

$$
n_{ab}(e,f)=\frac{N_{ab}(e,f)}{\min\bigl(|A(e)|,|B(f)|\bigr)}\label{eq:norm}
$$

聚合方式需要谨慎。对固定的 $e$ 逐节点计数可得四条求和规则，例如

$$
\sum_{f\neq e}N_{TT}(e,f)=\sum_{v\in T(e)}\bigl(k_{\rm out}(v)-1\bigr),\qquad
\sum_{f\neq e}N_{HT}(e,f)=\sum_{v\in H(e)}k_{\rm out}(v)
$$

其余两式同理。右端只依赖度序列，因此**在保度约束下四个通道的重叠总量是常数**；同质阶数下 \eqref{eq:norm} 的分母逐通道为常数，故无条件平均 $\langle n_{ab}\rangle$ 同样被度序列锁死，无法作为可调参数。可变的是重叠的堆积方式——是让少数超边对深度共享，还是让共享摊薄到大量超边对上。据此定义条件平均

$$
\alpha^{ab}_{\sigma\sigma'}=\Bigl\langle\, n_{ab}(e,f)\ \Bigm|\ N_{ab}(e,f)\ge1,\ \sigma(e)=\sigma,\ \sigma(f)=\sigma' \,\Bigr\rangle\label{eq:alpha}
$$

即在确实发生重叠的超边对上取平均，其中 $\sigma(e)=(|T(e)|,|H(e)|)$ 为超边的阶。记 $m_{ab}$ 为该通道中发生重叠的超边对数，则 $\alpha^{ab}=C_{ab}/m_{ab}$，其中常数 $C_{ab}$ 仅由度序列决定；调节 $\alpha$ 等价于调节 $m_{ab}$；这决定了第 V 节重连算法的实现方式。

式 \eqref{eq:alpha} 中阶指标取对角块（$\sigma=\sigma'$）即阶内重叠、取非对角块即阶间重叠，与无向情形下超边重叠矩阵的对角/非对角划分一致 [28]；本文只是在其每个矩阵元上再展开出四个方向通道。同质阶数下阶指标退化，以下略去并简记为 $\alpha^{ab}$。

互惠方向不由上述张量刻画。两条超边是否互相指回，取决于一组超边能否共同完成反向覆盖，这是超边级而非超边对级的判据。本文直接采用已有的强互惠定义 [34]：若存在一组超边，其尾集之并覆盖 $H(e)$ 且头集之并覆盖 $T(e)$，则称 $e$ 被强互惠；取

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

注意对有序超边对恒有 $N_{TH}(e,f)=N_{HT}(f,e)$：同质阶数下聚合量 $\alpha^{HT}$ 与 $\alpha^{TH}$ 逐实现相等（异质阶数下互为阶指标的转置），串联方向只携带一个独立分量，记 $\alpha^{\rightrightarrows}=\alpha^{HT}$。四个通道在聚合层面因而携带三个独立分量，按宇称重组：$\alpha^{\rightrightarrows}$ 与并联重叠 $\alpha^{\parallel}=\tfrac12(\alpha^{TT}+\alpha^{HH})$ 在 $\mathcal{R}$ 下为偶，而极性失衡

$$
\Delta\alpha=\tfrac12\bigl(\alpha^{TT}-\alpha^{HH}\bigr)
$$

为奇。$\alpha^{\parallel}$ 必须在扫描 $\alpha^{\rightrightarrows}$ 或 $\Delta\alpha$ 时固定，否则三种效应混在同一条扫描线上。

入出度相关 $r_{io}$ 在 $\mathcal{R}$ 下为**偶**：翻转互换每个节点的 $(k_{\rm in},k_{\rm out})$，而 \eqref{eq:rio} 的 Pearson 相关系数对其两个变量对称，故 $r_{io}(\mathcal{R}\mathcal{H})=r_{io}(\mathcal{H})$。据此可把结构量按宇称分类；本文所用结构量中，$\mathcal{R}$-奇者有三个：尾头基数之差 $\tau-\eta$、联合度分布关于对角线的反对称部分 $P(k_{\rm in},k_{\rm out})-P(k_{\rm out},k_{\rm in})$，以及极性失衡 $\Delta\alpha$。

这一分类对终态有直接约束。设结构系综 $\mathbb{P}[\mathcal{H}]$ 生成随机有向超图，$S[\mathcal{H}]$ 为 II.B 定义的条件期望爆发规模；由于种子均匀选取，$S$ 的定义本身不引入方向偏好。

> **命题 1.** 若 $\mathbb{P}$ 在 $\mathcal{R}$ 下不变，即 $\mathbb{P}[\mathcal{R}\mathcal{H}]=\mathbb{P}[\mathcal{H}]$ 对所有 $\mathcal{H}$ 成立，则 $\langle S\rangle_{\mathbb{P}}=\langle S\circ\mathcal{R}\rangle_{\mathbb{P}}$。

证明只需一次换元：$\mathcal{R}$ 是结构空间上的可测对合，故 $\sum_{\mathcal{H}}\mathbb{P}[\mathcal{H}]\,S[\mathcal{R}\mathcal{H}]=\sum_{\mathcal{H}}\mathbb{P}[\mathcal{R}\mathcal{H}]\,S[\mathcal{H}]=\sum_{\mathcal{H}}\mathbb{P}[\mathcal{H}]\,S[\mathcal{H}]$。

命题 1 的逆否形式给出方向对称性破缺的必要条件：**系综层面的破缺要求生成系综本身破坏 $\mathcal{R}$-不变性**；在本文所用结构量中，能承载这一破坏的只有上述三个 $\mathcal{R}$-奇量。特别地，$\tau=\eta$ 且 $P$ 关于对角线对称的位形模型系综是 $\mathcal{R}$-不变的，而以 $\mathcal{R}$-偶量为目标的定向重连（如 $r_{io}$ 扫描）保持该不变性，故单独扫描 $r_{io}$ 不可能产生系综层面的破缺——尽管 $r_{io}$ 确实移动阈值。这一推论决定了第 VI 节破缺实验的设计：固定前两个 $\mathcal{R}$-奇量为零，只扫 $\Delta\alpha$，并以 $r_{io}$ 扫描作为无破缺的阴性对照。需要强调命题 1 是关于系综平均的陈述；对单个实现 $\mathcal{H}$，$S[\mathcal{H}]$ 与 $S[\mathcal{R}\mathcal{H}]$ 一般不等，其差为有限尺寸涨落。

---

## 参考文献

1. M. E. J. Newman, S. H. Strogatz, D. J. Watts, *Physical Review E*, **64**, 026118 (2001). doi:10.1103/PhysRevE.64.026118
2. M. E. J. Newman, *Physical Review E*, **66**, 016128 (2002). doi:10.1103/PhysRevE.66.016128
3. P. Grassberger, *Mathematical Biosciences*, **63**, 157–172 (1983). doi:10.1016/0025-5564(82)90036-0
4. R. Pastor-Satorras, A. Vespignani, *Physical Review Letters*, **86**, 3200–3203 (2001). doi:10.1103/PhysRevLett.86.3200
5. E. Kenah, J. M. Robins, *Physical Review E*, **76**, 036113 (2007). doi:10.1103/PhysRevE.76.036113
6. E. Volz, *Journal of Mathematical Biology*, **56**, 293–310 (2008). doi:10.1007/s00285-007-0116-4
7. J. C. Miller, A. C. Slim, E. M. Volz, *Journal of the Royal Society Interface*, **9**, 890–906 (2012). doi:10.1098/rsif.2011.0403
8. P. Van Mieghem, R. van de Bovenkamp, *Physical Review Letters*, **110**, 108701 (2013). doi:10.1103/PhysRevLett.110.108701
9. R. Pastor-Satorras, C. Castellano, P. Van Mieghem, A. Vespignani, *Reviews of Modern Physics*, **87**, 925–979 (2015). doi:10.1103/RevModPhys.87.925
10. S. N. Dorogovtsev, J. F. F. Mendes, A. N. Samukhin, *Physical Review E*, **64**, 025101 (2001). doi:10.1103/PhysRevE.64.025101
11. M. Boguñá, M. Á. Serrano, *Physical Review E*, **72**, 016106 (2005). doi:10.1103/PhysRevE.72.016106
12. L. A. Meyers, M. E. J. Newman, B. Pourbohloul, *Journal of Theoretical Biology*, **240**, 400–418 (2006). doi:10.1016/j.jtbi.2005.10.004
13. X. J. Xu, J. Y. Li, X. Fu, L. J. Zhang, *Scientific Reports*, **8** (2018). doi:10.1038/s41598-018-22508-1
14. M. E. J. Newman, *Physical Review Letters*, **89**, 208701 (2002). doi:10.1103/PhysRevLett.89.208701
15. F. Battiston, G. Cencetti, I. Iacopini, V. Latora, M. Lucas, A. Patania, J. G. Young, G. Petri, *Physics Reports*, **874**, 1–92 (2020). doi:10.1016/j.physrep.2020.05.004
16. F. Battiston, E. Amico, A. Barrat, G. Bianconi, G. Ferraz de Arruda, B. Franceschiello, I. Iacopini, S. Kéfi et al., *Nature Physics*, **17**, 1093–1098 (2021). doi:10.1038/s41567-021-01371-4
17. A. R. Benson, R. Abebe, M. T. Schaub, A. Jadbabaie, J. Kleinberg, *Proceedings of the National Academy of Sciences*, **115**, E11221–E11230 (2018). doi:10.1073/pnas.1800683115
18. I. Iacopini, G. Petri, A. Barrat, V. Latora, *Nature Communications*, **10**, 2485 (2019). doi:10.1038/s41467-019-10431-6
19. G. Ferraz de Arruda, G. Petri, Y. Moreno, *Physical Review Research*, **2**, 023032 (2020). doi:10.1103/PhysRevResearch.2.023032
20. G. St-Onge, H. Sun, A. Allard, L. Hébert-Dufresne, G. Bianconi, *Physical Review Letters*, **127**, 158301 (2021). doi:10.1103/PhysRevLett.127.158301
21. D. J. Watts, *Proceedings of the National Academy of Sciences*, **99**, 5766–5771 (2002). doi:10.1073/pnas.082090499
22. D. Centola, M. Macy, *American Journal of Sociology*, **113**, 702–734 (2007). doi:10.1086/521848
23. N. W. Landry, J. G. Restrepo, *Chaos*, **30**, 103117 (2020). doi:10.1063/5.0020034
24. M. Mancastroppa, I. Iacopini, G. Petri, A. Barrat, *Nature Communications*, **14**, 6223 (2023). doi:10.1038/s41467-023-41887-2
25. J. H. Kim, K. I. Goh, *Physical Review Letters*, **132**, 087401 (2024). doi:10.1103/PhysRevLett.132.087401
26. G. St-Onge, I. Iacopini, V. Latora, A. Barrat, G. Petri, A. Allard, L. Hébert-Dufresne, *Communications Physics*, **5**, 25 (2022). doi:10.1038/s42005-021-00788-w
27. F. Malizia, S. Lamata-Otín, M. Frasca, V. Latora, J. Gómez-Gardeñes, *Nature Communications*, **16** (2025). doi:10.1038/s41467-024-55506-1
28. S. Lamata-Otín, F. Malizia, V. Latora, M. Frasca, J. Gómez-Gardeñes, *Physical Review E*, **111**, 034302 (2025). doi:10.1103/PhysRevE.111.034302
29. F. Malizia, A. Guzmán, I. Iacopini, I. Z. Kiss, *Physical Review Letters* (2025). doi:10.1103/z3d5-94zb
30. G. Bianconi, S. N. Dorogovtsev, *Physical Review E*, **109**, 014306 (2024). doi:10.1103/PhysRevE.109.014306
31. G. Ferraz de Arruda, A. Aleta, Y. Moreno, *Nature Reviews Physics*, **6**, 468–482 (2024). doi:10.1038/s42254-024-00733-0
32. G. Gallo, G. Longo, S. Pallottino, S. Nguyen, *Discrete Applied Mathematics*, **42**, 177–201 (1993). doi:10.1016/0166-218X(93)90045-P
33. S. Klamt, U. U. Haus, F. Theis, *PLoS Computational Biology*, **5**, e1000385 (2009). doi:10.1371/journal.pcbi.1000385
34. Q. F. Lotito, A. Vendramini, A. Montresor, F. Battiston, *Communications Physics*, **9** (2026). doi:10.1038/s42005-025-02472-9
35. S. Kim, M. Choe, J. Yoo, K. Shin, *Data Mining and Knowledge Discovery* (2023). doi:10.1007/s10618-023-00955-3
36. G. Lee, M. Choe, K. Shin, *Proceedings of the Web Conference 2021* (2021). doi:10.1145/3442381.3450010
37. P. S. Chodrow, *Journal of Complex Networks*, **8**, cnaa018 (2020). doi:10.1093/comnet/cnaa018
38. K. Nakajima, K. Shudo, N. Masuda, *IEEE Transactions on Network Science and Engineering*, **9**, 1139–1153 (2022). doi:10.1109/TNSE.2021.3133380
39. Y. J. Kraakman, C. Stegehuis, *Discrete Mathematics*, **349**, 114961 (2026). doi:10.1016/j.disc.2025.114961
40. M. Abuissa, M. Riondato, E. Upfal, *Data Mining and Knowledge Discovery*, **40** (2026). doi:10.1007/s10618-026-01209-8
41. J. Li, X. Wu, J. Lü, L. Lei, *Communications Physics*, **7** (2024). doi:10.1038/s42005-024-01614-9
42. D. T. Gillespie, *The Journal of Physical Chemistry*, **81**, 2340–2361 (1977). doi:10.1021/j100540a008
43. S. N. Dorogovtsev, A. V. Goltsev, J. F. F. Mendes, *Physical Review Letters*, **96**, 040601 (2006). doi:10.1103/PhysRevLett.96.040601
44. G. J. Baxter, S. N. Dorogovtsev, A. V. Goltsev, J. F. F. Mendes, *Physical Review E*, **82**, 011103 (2010). doi:10.1103/PhysRevE.82.011103
45. H. Sun, G. Bianconi, *Physical Review E*, **104**, 034306 (2021). doi:10.1103/PhysRevE.104.034306
