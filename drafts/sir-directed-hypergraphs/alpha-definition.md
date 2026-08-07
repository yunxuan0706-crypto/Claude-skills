# 有向超边重叠 α 的可计算定义

> 解决大纲 II.2 / 风险 2.2（α 的操作化定义与解析入口）
> 供 II.2 正文、V.4 扫描工具、以及 I.2 边界声明引用

---

## 0. 为什么原定义不能用

大纲原文把 α 描述为"尾集内部两两节点间的有向关系与超边整体方向一致的比例"。这个描述有三处使它无法进入理论：

1. **它不是超边之间的量。** 它描述的是**单条超边内部**尾集的成对投影，而"重叠"在文献中一律指**两条超边之间**共享节点的程度 [malizia2025overlap, lamataotin2025matrix]。用同一个词指不同的对象，审稿人第一轮就会要求澄清。
2. **它需要一个不存在的成对底图。** "尾集内部两两节点间的有向关系"预设了节点之间已有成对有向边；但在纯超图模型里没有这样的底图，必须先选一个投影规则（clique 展开？star 展开？），而不同投影会给出不同的 α，结论随之改变。
3. **它没有解析入口。** 正因为它不是超边间的量，它不会出现在任何自洽方程的分支结构里，因此无法进入 PGF 或房室闭合——这正是风险 2.2 记录的困难。

下面给出的定义把 α 重新安置成**超边之间的、有向分解的重叠**。它同时满足三个要求：无向极限下**精确退化**为顶刊现行定义、在保度约束下**确实可调**、并且有**明确的解析入口**。

---

## 1. 定义

### 1.1 重叠列联矩阵

设有向超图 $\mathcal{H}=(V,E)$，超边 $e=(T(e),H(e))$，$T(e)\cap H(e)=\varnothing$。

对任意**有序**超边对 $(e,f)$，$e\neq f$，定义 $2\times2$ **重叠列联矩阵**

$$
N(e,f)=\begin{pmatrix} N_{TT} & N_{TH} \\ N_{HT} & N_{HH}\end{pmatrix},
\qquad N_{ab}(e,f)=\bigl|\,A(e)\cap B(f)\,\bigr|,\quad a,b\in\{T,H\}
$$

其中 $A=T$ 当 $a=T$、$A=H$ 当 $a=H$，$B$ 同理。四个分量各有独立的动力学含义：

| 分量 | 结构含义 | 动力学含义 |
|---|---|---|
| $N_{TT}$ | **共发**（co-sending）：一个节点同时是 $e$ 与 $f$ 的施压方 | 源冗余；同一感染源向两个不同头集扇出 |
| $N_{HH}$ | **共收**（co-receiving）：一个节点同时是两条超边的受压方 | 剂量汇聚；群体阈值 $\theta$ 在此累积**跨超边**的压力 |
| $N_{HT}$ | **串联正向**：$e$ 的头集落入 $f$ 的尾集 | 传播链本身：$e$ 激活后经该节点驱动 $f$ |
| $N_{TH}$ | **串联反向**：$f$ 的头集落入 $e$ 的尾集 | 反向链；与 $N_{HT}$ 同时非零即为高阶互惠 [lotito2026directed] |

"共发/共收"这组词不是我们自造的：Lotito 等在有向超图的微观组织框架里正是用 **co-sending / co-receiving** 来刻画源集与靶集的重叠 [lotito2026directed]。本定义与该词汇表一致。

### 1.2 归一化与阶分辨聚合

归一化到 $[0,1]$：

$$
n_{ab}(e,f)=\frac{N_{ab}(e,f)}{\min\bigl(|A(e)|,|B(f)|\bigr)}
$$

记超边的**阶** $\sigma(e)=(|T(e)|,|H(e)|)$。定义**有向重叠张量**

$$
\boxed{\ \alpha^{ab}_{\sigma\sigma'}=\Bigl\langle\, n_{ab}(e,f)\ \Bigm|\ N_{ab}(e,f)\ge 1,\ \sigma(e)=\sigma,\ \sigma(f)=\sigma' \,\Bigr\rangle\ }
$$

即：**在确实发生重叠的超边对上，重叠的平均深度**。条件化（$N_{ab}\ge1$）不是修饰，是必需的——理由见 §3。

阶对角块 $\sigma=\sigma'$ 是**阶内重叠**（intra-order），非对角块是**阶间重叠**（inter-order）。这一对角/非对角的划分与 Lamata-Otín 等的超边重叠矩阵完全一致 [lamataotin2025matrix]；本定义在其每个矩阵元上再展开出 $2\times2$ 的方向通道，因此是**张量**而非矩阵。

> 术语提醒：大纲把 α 叫作"阶间重叠"。在上式里，"阶间"专指 $\sigma\neq\sigma'$。若主线采用同质阶数（$|T(e)|=\tau$、$|H(e)|=\eta$ 对所有 $e$），则**只剩阶内块**，此时应改称**有向超边重叠**，不要再写"阶间"，否则与 [lamataotin2025matrix, malizia2025gbcm] 的术语冲突。异质阶数留到附录时"阶间"才重新有意义。

### 1.3 方向对称性分解

全局方向翻转 $\mathcal{R}$：对所有 $e$ 交换 $T(e)\leftrightarrow H(e)$。它是空间上的对合（$\mathcal{R}^2=\mathrm{id}$），四个通道在其下的变换为

$$
\mathcal{R}:\quad \alpha^{TT}\leftrightarrow\alpha^{HH},\qquad \alpha^{HT}\leftrightarrow\alpha^{TH}
$$

据此把四通道重组为**宇称本征量**（略去阶指标）：

$$
\alpha^{\rightrightarrows}=\tfrac12\bigl(\alpha^{HT}+\alpha^{TH}\bigr)
\qquad\text{（串联／同向，}\mathcal{R}\text{-偶）}
$$
$$
\alpha^{\rightleftarrows}=\bigl\langle \min\bigl(n_{HT},n_{TH}\bigr)\bigr\rangle
\qquad\text{（互惠／反向，}\mathcal{R}\text{-偶）}
$$
$$
\alpha^{\parallel}=\tfrac12\bigl(\alpha^{TT}+\alpha^{HH}\bigr)
\qquad\text{（并联，}\mathcal{R}\text{-偶）}
$$
$$
\Delta\alpha=\tfrac12\bigl(\alpha^{TT}-\alpha^{HH}\bigr)
\qquad\text{（极性失衡，}\boldsymbol{\mathcal{R}\text{-奇}}\text{）}
$$

这就把大纲的两个记号 $\alpha^{\rightrightarrows}/\alpha^{\rightleftarrows}$ 落到了确定的可计算式上，同时暴露出**大纲遗漏的两个自由度** $\alpha^{\parallel}$ 与 $\Delta\alpha$。$\alpha^{\parallel}$ 必须在扫描 $\alpha^{\rightrightarrows}$ 时**固定住**，否则扫描把三种效应混在一起；$\Delta\alpha$ 则是下一节的主角。

---

## 2. 三条必须写进正文的性质

### 性质 1（退化）：无向极限下精确回收顶刊定义

令 $\bar e=T(e)\cup H(e)$ 为超边的支撑集。由 $T(e)\cap H(e)=\varnothing$，

$$
|\bar e\cap\bar f|=N_{TT}+N_{TH}+N_{HT}+N_{HH}
$$

因此通道求和后的张量**逐超边对地**等于 [malizia2025overlap, lamataotin2025matrix] 的无向超边重叠。

> **这条性质是 I.2 中"α 是 refinement 而非另起炉灶"这一增量声明的全部依据。** 没有它，审稿人有权认为我们换了个定义使自己的结论无法与已有文献对照。

### 性质 2（宇称）：方向对称性破缺的必要条件

设结构系综 $P[\mathcal{H}]$ 生成随机有向超图，$S[\mathcal{H}]$ 为 SIR 终态爆发规模。

**命题.** 若 $P$ 在 $\mathcal{R}$ 下不变，即 $P[\mathcal{R}\mathcal{H}]=P[\mathcal{H}]$，则 $\langle S\rangle_{P}$ 在方向翻转下不变；因而**系综层面的方向对称性破缺要求至少一个 $\mathcal{R}$-奇不变量非零**。

*证明.* $\mathcal{R}$ 是有向超图空间上的可测对合。$\langle S\rangle_{\mathcal{R}P}=\sum_{\mathcal{H}}P[\mathcal{R}\mathcal{H}]S[\mathcal{H}]=\sum_{\mathcal{H}}P[\mathcal{H}]S[\mathcal{H}]=\langle S\rangle_P$，第二个等号用 $P$ 的 $\mathcal{R}$-不变性并对 $\mathcal{H}\mapsto\mathcal{R}\mathcal{H}$ 换元。∎

把系综的 $\mathcal{R}$-奇不变量列全，就得到一张**破缺来源清单**：

| $\mathcal{R}$-奇不变量 | 说明 |
|---|---|
| $\tau-\eta$ | 尾集与头集基数不等。注意 $\langle k_{\rm out}\rangle N=\tau\lvert E\rvert$、$\langle k_{\rm in}\rangle N=\eta\lvert E\rvert$，故 $\tau\neq\eta$ 自动使两个平均度不等 |
| $P(k_{\rm in},k_{\rm out})-P(k_{\rm out},k_{\rm in})$ | 联合度分布关于对角线的不对称部分（如 $\langle k_{\rm in}^2k_{\rm out}\rangle-\langle k_{\rm in}k_{\rm out}^2\rangle$） |
| $\Delta\alpha=\tfrac12(\alpha^{TT}-\alpha^{HH})$ | 共发重叠与共收重叠之差 |

**一个立即可用的推论：$r_{io}$ 是 $\mathcal{R}$-偶的，因此它本身不能产生方向对称性破缺。** 理由：$\mathcal{R}$ 把每个节点的 $(k_{\rm in},k_{\rm out})$ 互换，而 Pearson 相关系数对两个变量对称，故 $r_{io}(\mathcal{R}\mathcal{H})=r_{io}(\mathcal{H})$。

> 这个推论直接改写大纲 V.5 与 C4 的做法：**在 $\tau=\eta$、$P$ 关于对角线对称的设定下扫描 $r_{io}$，系综平均的爆发规模必然对称**，观测到的任何差异只是有限尺寸涨落。破缺实验必须显式打开一个 $\mathcal{R}$-奇量——最干净的做法是固定 $\tau=\eta$ 与 $P$ 的对称性，只扫 $\Delta\alpha$，从而把破缺**单独归因**到 $\Delta\alpha$。这把大纲里"严格证明为拉伸目标"的部分变成了一个可证的必要条件加一个干净的数值设计。

### 性质 3（可调性）：为什么必须条件化

**求和规则.** 对固定的 $e$，把 $N_{ab}(e,f)$ 对所有 $f\neq e$ 求和，按共享节点计数得

$$
\sum_{f\neq e}N_{TT}(e,f)=\sum_{v\in T(e)}\bigl(k_{\rm out}(v)-1\bigr),
\qquad
\sum_{f\neq e}N_{HH}(e,f)=\sum_{v\in H(e)}\bigl(k_{\rm in}(v)-1\bigr)
$$
$$
\sum_{f\neq e}N_{HT}(e,f)=\sum_{v\in H(e)}k_{\rm out}(v),
\qquad
\sum_{f\neq e}N_{TH}(e,f)=\sum_{v\in T(e)}k_{\rm in}(v)
$$

右端**只依赖度序列**。所以：

> **在严格保度的约束下，四个通道的重叠总量是常数，不可调。**

这是一条硬约束，而大纲 V.4（"在保度前提下调节 α"）默认了 α 可调却从未验证。若把 α 定义为无条件均值 $\langle n_{ab}\rangle$，则它正比于上述固定总量除以超边对数，**恒为常数，扫描无从谈起**——工具 2 会在实现阶段撞上这堵墙。

条件化解决了这个问题：固定总量下，仍可自由改变重叠的**堆积方式**——是让少数超边对深度共享，还是让共享摊薄到大量超边对上。形式化地，记 $m_{ab}=\#\{(e,f):N_{ab}\ge1\}$ 为发生重叠的对数，则

$$
\alpha^{ab}\ \propto\ \frac{\text{固定总量}}{m_{ab}}
$$

**调 α 等价于调 $m_{ab}$**，而 $m_{ab}$ 在保度纤维内确实可变。这同时给出扫描的物理图像与实现路径：*降低* α = 把重叠摊开到更多超边对；*提高* α = 把重叠集中到更少的对上。这与 Malizia 等"低阶内重叠 → 爆炸式转变"的结论在方向上可直接对照 [malizia2025overlap]。

---

## 3. 解析入口（风险 2.2）

重叠破坏局部树状假设，节点级 PGF 因而不闭合。可行路径是把自洽方程从**节点级抬到超边级**：

1. **群体级房室闭合。** Malizia 等的 group-based compartmental modeling（GBCM）正是为"跨阶相关"设计的均场框架，能解析地分离各阶贡献 [malizia2025gbcm]。把它有向化：状态变量取"超边 $e$ 尚未被激活"的概率，转移项按 §1.1 的四个通道分列。
2. **四通道转移算子。** 超边级自洽方程的线性化给出分块转移算子

$$
\mathcal{M}=\begin{pmatrix}\mathcal{M}^{TT} & \mathcal{M}^{TH}\\ \mathcal{M}^{HT} & \mathcal{M}^{HH}\end{pmatrix},
\qquad \mathcal{M}^{ab}\ \text{的权重由}\ \alpha^{ab}\ \text{与}\ \phi\ \text{给出}
$$

爆发阈值即 $\rho(\mathcal{M})=1$（谱半径）。$\alpha$ 于是**显式**进入阈值条件，风险 2.2 关闭。

3. **宇称的算子表述.** $\mathcal{R}$ 在 $\mathcal{M}$ 上的作用是反对角共轭 $\mathcal{M}\mapsto J\mathcal{M}^{\!\top}J$（$J$ 为通道交换）。$\Delta\alpha=0$ 时 $\mathcal{M}$ 在该共轭下不变，谱不变，阈值对称——与性质 2 一致。这也正是附录 D"张量后向算子"所需的算子语言。

**诚实标注：** 第 2 步中 $\mathcal{M}^{ab}$ 的权重需要 $\alpha^{ab}$ 之外的信息（重叠对的度分布），因此严格来说 $\alpha$ 是**闭合的一阶控制量**而非完备描述；GBCM 在无向情形也是这样的近似 [malizia2025gbcm]。正文应写成"至一阶闭合"，不要写成精确。

---

## 4. 对下游各节的影响

| 位置 | 需要的改动 |
|---|---|
| **II.2** | 用 §1 替换现有描述；同质阶数下术语从"阶间重叠"改为"有向超边重叠" |
| **V.4** | 扫描目标从 $\alpha$ 的均值改为条件均值；实现上以 $m_{ab}$ 为调节手柄；必须同时固定 $\alpha^{\parallel}$ |
| **V.5 / C4** | 按性质 2 重新设计：固定 $\tau=\eta$ 与 $P$ 的对角对称性，只扫 $\Delta\alpha$，把破缺单独归因 |
| **IV.2–IV.3** | 阈值条件写成 $\rho(\mathcal{M})=1$，$\alpha^{ab}$ 显式出现在 $\mathcal{M}$ 中 |
| **附录 D** | 后向算子即 §3 的 $\mathcal{M}$；宇称由反对角共轭给出 |
| **I.2** | 增量声明依据性质 1（退化）与性质 2（宇称），二者都可证 |

## 5. 仍需拍板

- **归一化的分母。** 上文取 $\min(|A(e)|,|B(f)|)$。同质阶数下这就是 $\tau$、$\eta$ 或 $\min(\tau,\eta)$，无歧义；异质阶数下另有取 $\sqrt{|A||B|}$（余弦型）的选项，两者在异质性强时结论可能不同。建议主线用 $\min$，附录做敏感性检查。
- **$\alpha^{\rightleftarrows}$ 的形式。** 上文取 $\langle\min(n_{HT},n_{TH})\rangle$。另一可选形式是 Lotito 等互惠性定义下的 exact/strong/weak 三档 [lotito2026directed]。若第 VII 部分要与真实数据对照，建议直接采用其 strong 定义以便复用其测量口径。
