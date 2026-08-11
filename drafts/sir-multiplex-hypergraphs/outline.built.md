# 多重超图上的 SIR 传播：多渠道群体交互的仓室理论

**工作标题**（备选）：Compartmental theory of SIR epidemics on multiplex hypergraphs

---

## 0–1. 定位与引言（略写）

同一种病经由多条**群体渠道**传播，个体在任一渠道被感染即在所有渠道中都是感染态。
成对多层已解决 [3, 4, 11]；
高阶单层已解决 [8]；高阶多层只做到渗流
[9] 或均场且各层角色不同
[5, 10]。本文补上**超出均场的含时
仓室理论**，并把爆发阈值写成一个层间矩阵的谱半径。

---

## 2. 模型

### 2.1 多重超图与记号

节点集 $V$，$N=|V|$。层 $a=1,\dots,M$，层 $a$ 是 $V$ 上的超图
$\mathcal{H}_a=(V,E_a)$，主线取同质基数 $|e|=m_a$（异质推广见附录 A）。

- 层度 $k^{(a)}(v)=|\{e\in E_a:v\in e\}|$，写成向量 $\mathbf{k}(v)=(k^{(1)},\dots,k^{(M)})$
- **联合层度分布** $P(\mathbf{k})$ —— 本文的第一个结构输入
- 握手关系 $\sum_v k^{(a)}(v)=m_a|E_a|$，即 $\langle k^{(a)}\rangle N=m_a|E_a|$
- 多元生成函数
  $$
  \Psi(\mathbf{x})=\sum_{\mathbf{k}}P(\mathbf{k})\prod_{c}x_c^{k_c}\label{eq:psi}
  $$
- **层 $a$ 的余度生成函数**（沿层 $a$ 的一个成员槽到达节点，来路那个群体不再计入）
  $$
  \psi_a(\mathbf{x})=\frac{1}{\langle k^{(a)}\rangle}\sum_{\mathbf{k}}k_a\,P(\mathbf{k})\prod_c x_c^{\,k_c-\delta_{ac}}\label{eq:psia}
  $$

$\psi_a$ 的下标 $-\delta_{ac}$ 就是余度相减。它与有向工作相反：那里尾侧偏置与入度
计数是两个不同的量，故无相减；这里成员槽与计数是同一个量，相减回来了。

### 2.2 SIR 动力学

节点状态 $x_v\in\{S,I,R\}$ **跨层共享**（"一处感染、处处感染"的形式化）。

层 $a$ 的群体 $e$ 中已感染成员数 $n_e(t)=|e\cap I(t)|$。$n_e\ge\theta_a$ 时群体激活，
向 $e$ 中**每个**易感成员各以速率 $\beta_a$ 独立传递；$n_e<\theta_a$ 时静默。
感染者以速率 $\mu$ 康复。易感节点 $u$ 的总风险率

$$
\Lambda_u(t)=\sum_{a}\ \sum_{e\in E_a:\,u\in e}\beta_a\,\mathbb{1}[n_e(t)\ge\theta_a]\label{eq:rate}
$$

以 $\mu$ 归一时间，$\lambda_a=\beta_a/\mu$；主线取 $\lambda_a=\lambda w_a$，
$w_a$ 为渠道权重，$\lambda$ 为唯一的标量控制参数。主线取 $\theta_a=1$，
$\theta_a\ge2$ 见 3.7。

**关键结构差异（决定后文一切）。** 无向群体中每个成员既是源也是靶。一个成员被感染
后群体激活，它可以感染同群其他成员；后者被感染又**延长**群体的激活时长，从而进一步
提高剩余成员被感染的概率。这条"群内级联"在有向情形不存在（尾集与头集不相交，头成员
被感染不延长激活），是本模型全部新内容的来源。

### 2.3 层间结构量

两条互相独立、作用在生成流程不同环节的量：

- **层间参与相关** $\rho_{ab}=\mathrm{corr}\bigl(k^{(a)},k^{(b)}\bigr)$。由联合度序列
  决定，与群体如何连接无关，只能在构造度序列的阶段调节。
- **层间群体重叠**
  $$
  o_{ab}=\frac{1}{|E_a||E_b|}\sum_{e\in E_a}\sum_{f\in E_b}\frac{|e\cap f|}{\min(m_a,m_b)}\label{eq:overlap}
  $$
  在度序列固定后由群体重排调节。

> **开放项（第一件要做的事）**：$o_{ab}$ 在保联合度约束下的可动范围。有向工作里同类
> 量只能挪动 $0.005$–$0.006$，仿真分辨不出，最后必须专门造重连工具。**先测再决定**
> 要不要把它当主线结构参数。

### 2.4 可约性

> **命题 2.1（可约性）.** 设节点状态跨层共享，且各层群体按 \eqref{eq:rate} 独立传递。
> 令 $\widetilde{\mathcal{H}}$ 为把全部 $E_a$ 并成一个**多重集**（两层可以含有节点集
> 相同的群体，合并后必须保留重数）、每条超边保留其层标签 $a$（从而保留
> $m_a,\theta_a,\beta_a$）所得的**单层带型超图**。则多重超图上的过程与
> $\widetilde{\mathcal{H}}$ 上的 SIR 过程同分布。

**证明.** 两个过程的状态空间同为 $\{S,I,R\}^V$。生成元是各超边贡献之和，而
\eqref{eq:rate} 中每条超边的贡献只依赖于该超边的成员集合、其 $(\theta_a,\beta_a)$
与当前构型，不依赖于它属于哪一层；康复项逐节点相同。故两个生成元逐项相等。$\square$

命题 2.1 必须放在正文显著位置。它不是缺陷，而是把贡献逼清楚的工具：

- 多层买到的**不是新动力学**，而是新的**结构自由度**——联合层度分布（含 $\rho_{ab}$）
  与层间群体重叠。单层带型超图的位形模型默认各型度独立；真实多层不独立。
- 它也给出一条实用推论：**任何单层带型超图的理论都自动求解本模型**，因此第 3 节的
  闭合可以在单层语言下推导、在多层语言下解读。

**失效边界**（引向第 7 节）：若节点把跨层暴露**合并后**再判阈值，风险率取
$\beta f\bigl(\sum_a n_{e_a}\bigr)$ 而非 $\sum_a\beta_a g_a(n_{e_a})$，
则生成元不再可分解，命题 2.1 失效。

---

## 3. 仓室理论

### 3.1 三个层级

全文的定量比较固定在这三者之间。

| 层级 | 仓室 | 忽略了什么 | 维数 |
|---|---|---|---|
| 度均场（基线） | 按层度类 $\mathbf{k}$ 的节点 | 群内相关、邻域相关、群内级联 | 层度类数 |
| **群体闭合（主结果）** | 群体的成员计数 $(s,i,r)$，逐层 | 群体间相关（$O(1/N)$） | $\sum_a\binom{m_a+1}{2}+2$ |
| 精确 | 全构型 | —— | $3^N$ |

### 3.2 空腔构造与两条假设

固定一个**检验节点** $u$，人为令其永久易感（空腔节点）。取一个 $u$ 所属的层 $a$
群体 $e$，定义边基变量

$$
\Phi_a(t)=\Pr\bigl[e\ \text{到}\ t\ \text{时刻尚未向}\ u\ \text{传递}\bigr]\label{eq:phi}
$$

**(H1) 局域树状。** 位形模型系综在 $N\to\infty$ 时局域收敛于超树，任一有限邻域内出现
回路的概率为 $O(1/N)$。据此 $u$ 所属的各群体的上游分支互不相交，因而

$$
S(t)=(1-\varepsilon)\,\Psi\bigl(\Phi_1(t),\dots,\Phi_M(t)\bigr)\label{eq:S}
$$

**多元生成函数在此进入**：层间相关只经由 $\Psi$（节点侧）与 $\psi_a$（成员侧，
见 \eqref{eq:psia}）的交叉项传入动力学，$\langle k^{(a)}k^{(b)}\rangle$ 一族因而
承载全部层间效应。除此之外没有别的通道——这正是 6.5 那条"闭合看不见层间群体重叠"
的推论的来源。

**(H2′) 成员的外部暴露相互独立且同分布。** 去掉 $e$ 之后，$e$ 的各成员落在互不相交的
分支中（由 (H1)），故它们**经由自身其他群体**被感染的过程相互独立；又由位形模型的
交换性，它们的层度是从 $k^{(a)}$ 偏置分布中独立同分布抽取的。

注意 (H2′) 与有向工作的 (H2) 不同：那里断言的是成员之间**完全**不相互作用；这里成员
之间通过 $e$ 本身强烈相互作用（群内级联），只有**外部**暴露独立。群内耦合不被忽略，
而是进入下面的计数状态。

### 3.3 群体计数仓室与方程组

**状态。** 对 $u$ 所属的层 $a$ 群体 $e$，以其**除 $u$ 外**的 $m_a-1$ 个成员的计数
$(s,i,r)$（易感/感染/康复，$s+i+r=m_a-1$）为仓室，定义未传递子概率

$$
x^{(a)}_{sir}(t)=\Pr\bigl[e\ \text{尚未向}\ u\ \text{传递}\ \wedge\ e\setminus\{u\}\ \text{处于}\ (s,i,r)\bigr]\label{eq:xsir}
$$

于是 $\Phi_a=\sum_{s+i+r=m_a-1}x^{(a)}_{sir}$，并记激活且未传递的部分

$$
\Phi^{A}_a(t)=\sum_{i\ge\theta_a}x^{(a)}_{sir}(t)\label{eq:phiA}
$$

（$u$ 易感，故不计入 $n_e$，群体的激活只看 $i$。）

**外部风险率。** 由 (H2′)，$e$ 的一个成员**经由其他群体**到 $t$ 仍易感的概率为
$(1-\varepsilon)\psi_a(\boldsymbol{\Phi})$，其中 $\psi_a$ 见 \eqref{eq:psia}。
故其瞬时外部风险率

$$
h_a(t)=-\frac{\mathrm{d}}{\mathrm{d}t}\ln\bigl[(1-\varepsilon)\psi_a(\boldsymbol{\Phi})\bigr]
=\sum_{c}\beta_c\,\Phi^{A}_c\,\frac{\partial_c\psi_a(\boldsymbol{\Phi})}{\psi_a(\boldsymbol{\Phi})}\label{eq:h}
$$

末一步用了下面的 \eqref{eq:phidot}。这是有向工作中标量 $h=\beta\Phi_A\psi'/\psi$ 的
**多层推广**：单个导数变成对所有层求和的偏导数。

**主方程。** 一个易感成员的总感染率是外部与群内之和 $h_a+\beta_a\mathbb{1}[i\ge\theta_a]$；
群体向 $u$ 的传递以速率 $\beta_a\mathbb{1}[i\ge\theta_a]$ 把概率移出"未传递"类。据此

$$
\begin{aligned}
\dot x^{(a)}_{sir}=\;&\bigl(h_a+\beta_a\mathbb{1}[i-1\ge\theta_a]\bigr)(s+1)\,x^{(a)}_{s+1,i-1,r}
-\bigl(h_a+\beta_a\mathbb{1}[i\ge\theta_a]\bigr)s\,x^{(a)}_{sir}\\
&+\mu\bigl[(i+1)x^{(a)}_{s,i+1,r-1}-i\,x^{(a)}_{sir}\bigr]
-\beta_a\mathbb{1}[i\ge\theta_a]\,x^{(a)}_{sir}
\end{aligned}\label{eq:master}
$$

（越界下标记为零。）等价地，添一个"已传递"的吸收态即得一条标准的连续时间 Markov 链，
$x^{(a)}_{sir}$ 是其未被吸收部分的子概率；这是以速率
$\beta_a\mathbb{1}[i\ge\theta_a]$ 消灭（killing）的 Markov 链。

对 \eqref{eq:master} 求和，$h_a$ 项、$\beta_a$ 群内项与 $\mu$ 项各自守恒相消，只余

$$
\dot\Phi_a=-\beta_a\,\Phi^{A}_a\label{eq:phidot}
$$

与 \eqref{eq:h} 自洽。节点侧由 \eqref{eq:S} 与 $\dot R=\mu I$、$I=1-S-R$ 闭合。

**初始条件**（均匀初始感染比例 $\varepsilon$）：

$$
x^{(a)}_{s,i,0}(0)=\binom{m_a-1}{i}(1-\varepsilon)^{s}\varepsilon^{i},\qquad s+i=m_a-1
$$

其余为零，$\Phi_a(0)=1$。

**维数。** 层 $a$ 贡献 $\binom{m_a+1}{2}$ 个状态（$s+i+r=m_a-1$ 的非负解数），加上
$R$ 与 $I$ 两个标量，共 $\sum_a\binom{m_a+1}{2}+2$，**与 $N$ 无关**。
$M=2$、$m_1=3,m_2=4$ 时是 $6+10+2=18$ 个方程。

### 3.4 闭合的地位

> **命题 3.1（树极限精确性）.** 在 (H1) 与 (H2′) 下，\eqref{eq:S}、\eqref{eq:h}–\eqref{eq:phidot}
> 给出的 $S(t)$ 在 $N\to\infty$ 时与精确过程一致。

**证明要点（目前只是论证，尚非严格证明；见下）.**（i）由 (H1)，$u$ 所属各群体的上游分支互不相交，\eqref{eq:S} 的乘积形式
（即 $\Psi$）成立。（ii）由 (H2′)，$e$ 的 $m_a-1$ 个成员的外部感染时刻相互独立且共享
同一生存函数 $(1-\varepsilon)\psi_a(\boldsymbol{\Phi})$；对独立同分布个体，其"尚存活
人数"是以**总体风险率** $h_a=-\mathrm{d}\ln(\cdot)/\mathrm{d}t$ 为速率的 Markov 计数
过程，故 \eqref{eq:master} 中的 $s\,h_a$ 项精确。（iii）群内传染只依赖于 $i$，
且对全体易感成员对称，故计数 $(s,i,r)$ 是充分统计量。$\square$

> **这一条更正了本大纲上一版的说法。** 上一版称"无向情形下 (H2) 失效、闭合只是近似"。
> 把方程写出来之后可以看到：失效的是**成员之间无相互作用**这一条，但群内相互作用被
> 计数状态完整地吸收了，真正需要的只有**外部暴露独立**(H2′)，而它在树极限下成立。
> 结论：闭合的地位与有向工作相同，仍须由 4.2 的规模扫描实测确认。
>
> **数值证据（已做，见 `tools/scaling_run.py`）。** $\lambda=0.35$、$m=(3,4)$、
> $\varepsilon=0.02$、每点 $300$–$500$ 次实现下，闭合与精确仿真在 $t\le12$ 上的平均
> 偏差随 $N$ 单调下降：$N=500,1000,2000,4000,8000$ 依次为
> $4.95,1.84,1.23,0.89,0.29$（单位 $10^{-3}$），到 $N=8000$ 已落到蒙特卡洛噪声底
> （$3.2\times10^{-4}$）以下、不可分辨。幂律拟合给出指数 $-0.92$（取不同子集为
> $-0.80$ 至 $-1.00$），与位形模型中短回路密度的 $O(1/N)$ 一致。注意每个测量值都含噪声
> 贡献，故是真实偏差的**上界**，真实衰减只会更陡。
>
> **但仍要诚实**：数值支持 $\ne$ 证明。上面三步是论证纲要而非严格证明，正文宜写成
> "在 (H1)(H2′) 下我们论证并数值验证了闭合的渐近精确性，严格证明（位形模型的局域弱
> 收敛 + 交换性）留待后续"。不要写成定理。

### 3.5 群内级联与爆发阈值

**这一小节是全文的支点。**

$\theta_a=1$、$\varepsilon\to0$ 时线性化。与有向情形的关键差别在于：一个成员被感染后，
群体在**其后续的群内级联期间**持续激活，因此"一个种子在一个群体内造成的感染数"
不是 $(m-1)T$，而要由一个小型 CTMC 决定。

**群内级联期望 $C(m,\lambda)$。** 考虑一个孤立的 $m$ 元群体，初始一个感染、$m-1$ 个
易感。记 $u(i,s)$ 为从"$i$ 个感染、$s$ 个易感"出发**此后**新增感染数的期望，则
$u(0,s)=u(i,0)=0$，且对 $i,s\ge1$（以 $\mu$ 归一时间）

$$
u(i,s)=\frac{\alpha_i}{\alpha_i+i}\bigl[1+u(i+1,s-1)\bigr]+\frac{i}{\alpha_i+i}\,u(i-1,s),
\qquad \alpha_i=\lambda s\,\mathbb{1}[i\ge\theta]\label{eq:cascade}
$$

$$
C(m,\lambda)=u(1,m-1)\label{eq:C}
$$

> **必须交代的文献关系（新颖性风险最高的一处）。** $C(m,\lambda)$ 正是经典
> **家庭模型**（households / two levels of mixing）中"户内最终规模"的期望，而
> $\rho(\mathsf{N})=1$ 的结构正是 Ball 等的 $R_*$ =（全局感染数）$\times$（户内
> 爆发规模）[1, 6]；带网络结构的重叠版本亦已有
> [2]，并被列为网络流行病学的公开挑战之一
> [7]。本文与之的区别必须写清楚且只能写这三条：
> (i) 群体不构成划分而是**任意重叠**，节点在层 $a$ 中属于 $k^{(a)}$ 个群体；
> (ii) **多型层结构**与联合层度分布 $P(\mathbf{k})$，阈值因而是矩阵谱半径而非标量；
> (iii) 群体阈值 $\theta_a\ge2$ 的激活规则。若不主动交代，审稿人会判定重复。

- $m=2$：$C=\lambda/(1+\lambda)=T$，恰为成对传递概率，**级联消失**。
- $m\ge3$：$C>(m-1)T$。数值上 $\lambda=1$ 时 $C/(m-1)T=1.11,1.22,1.31$
  （$m=3,4,5$）；$\lambda=0.5$ 时为 $1.10,1.20,1.30$。级联在 $m=5$ 已达三成。

**阈值矩阵。** 分支过程必须定义在**群体层面**：每个新感染节点"点燃"它的其余各群体，
每个被点燃的群体产出 $C$ 个新感染成员，这些成员再点燃各自的其余群体。
（若按节点层面朴素地数"某人直接传染了谁"，群内级联的后继会被重复计入。）
以"感染经由哪一层传来"为型，经由层 $a$ 感染的节点，其层 $b$ 群体数的
期望为 $\langle k^{(a)}k^{(b)}\rangle/\langle k^{(a)}\rangle-\delta_{ab}$（余度相减来自
\eqref{eq:psia}），每个这样的群体贡献 $C(m_b,\lambda_b)$ 个新感染，故

$$
\boxed{\;\mathsf{N}_{ab}=C(m_b,\lambda_b)\Bigl[\frac{\langle k^{(a)}k^{(b)}\rangle}{\langle k^{(a)}\rangle}-\delta_{ab}\Bigr]\;}\label{eq:N}
$$

> **命题 3.2（阈值）.** $\theta_a=1$ 时，爆发阈值由
> $$\rho(\mathsf{N})=1$$
> 确定，$\rho$ 为谱半径。$\lambda<\lambda_c$ 时 $\rho<1$，无宏观爆发。

**必做的退化检验（每一条都要写进正文）：**

| 退化 | 应回到 |
|---|---|
| $M=1$ | 单层超图 SIR 的已知阈值 |
| 全部 $m_a=2$ | $C=T$，\eqref{eq:N} 回到多层网络 SIR 的已知阈值 [4] |
| $M=1,m=2$ | $\lambda_c$ 由 $\langle k\rangle/(\langle k^2\rangle-\langle k\rangle)$ 给出 |

**$M=2$ 的闭式。** 记 $\mathsf{N}$ 的四个元素，

$$
\rho=\tfrac12\Bigl[\mathsf{N}_{11}+\mathsf{N}_{22}+\sqrt{(\mathsf{N}_{11}-\mathsf{N}_{22})^2+4\mathsf{N}_{12}\mathsf{N}_{21}}\Bigr]=1\label{eq:M2}
$$

代入 \eqref{eq:N} 并以 $\lambda$ 为未知量求解即得 $\lambda_c$。

> **命题 3.3（层间耦合永不抬高阈值）.** $\mathsf{N}$ 逐元非负，故
> $\rho(\mathsf{N})\ge\max_a\mathsf{N}_{aa}$，而 $\mathsf{N}_{aa}$ 恰是孤立层 $a$ 的
> 分支数。在 $\lambda_a=\lambda w_a$ 的共同标度下（$\rho$ 关于 $\lambda$ 单调增），
> 多重超图的阈值因而不高于任一单层的阈值；若 $\mathsf{N}$ 不可约且存在非零非对角元，
> 则严格更低。

**证明.** 非负矩阵的谱半径不小于任一对角元；不可约时 Perron–Frobenius 给出对元素的
严格单调性。$\square$

命题 3.3 的可检验推论，也是本文最有传播力的一句话：**两条各自亚临界的弱渠道，
合起来可以超临界**，且发生条件由 $\mathsf{N}_{12}\mathsf{N}_{21}$ 显式给出。
第 6 节要做出这个协同区的相图。已在 $P=\{(2,2),(3,3)\}$、$m_1=m_2=3$ 上验证：
$\lambda=0.13$ 时两层各自的 $\mathsf{N}_{aa}=0.386$（深度亚临界），
合起来 $\rho(\mathsf{N})=1.013$（超临界）。

### 3.6 与度均场的比较

度均场：设层度类 $\mathbf{k}$ 的感染概率 $i_\mathbf{k}$，群体激活概率取系综平均
$\Theta_a=1-(1-\phi_a)^{m_a-1}$，$\phi_a=\sum_{\mathbf{k}}k_aP(\mathbf{k})i_\mathbf{k}/\langle k^{(a)}\rangle$，

$$
\dot s_\mathbf{k}=-s_\mathbf{k}\sum_a k_a\beta_a\Theta_a,\qquad
\dot i_\mathbf{k}=s_\mathbf{k}\sum_a k_a\beta_a\Theta_a-\mu\, i_\mathbf{k}\label{eq:mf}
$$

线性化（$\Theta_a\simeq(m_a-1)\phi_a$）给出

$$
\mathsf{N}^{\rm MF}_{ab}=(m_b-1)\,\lambda_b\,\frac{\langle k^{(a)}k^{(b)}\rangle}{\langle k^{(a)}\rangle}\label{eq:Nmf}
$$

> **命题 3.4（均场的三处偏差）.** \eqref{eq:Nmf} 与 \eqref{eq:N} 的差别恰为三项，
> 且方向不同：
> 1. $\lambda_b$ 取代 $T_b=\lambda_b/(1+\lambda_b)$ —— 均场把整个感染期的传递概率
>    当成 $\beta/\mu$，忽略"群体一旦向某成员传递即对它用尽"，**高估**；
> 2. 缺少余度相减 $-\delta_{ab}$ —— 均场允许沿来路那个群体折返，**高估**；
> 3. 缺少级联因子 $C/[(m_b-1)T_b]>1$ —— 均场的线性化把各成员的贡献独立相加，
>    看不见群内级联延长激活，**低估**。

前两项在 $m=2$ 时退化为有向工作里"均场丢掉分母的 $-1$"的同一件事；第三项是高阶
特有的，且方向相反。**实测（$k\in\{3,5\}$、$M=2$、$m_1=m_2=m$）表明前两项占优，
净效应不变号**：$\lambda^{\rm MF}_c/\lambda_c$ 在 $m=2,3,4,6,8$ 处依次为
$0.76,0.85,0.88,0.90,0.91$——均场始终低估阈值，但差距随 $m$ 单调缩小，
因为级联项部分抵消。第 4.4 节把三项逐项量化，并检验"是否存在使净效应变号的
$(m,\lambda,P(\mathbf{k}))$ 区域"——**这是一个开放问题，不是已知结论**。

### 3.7 $\theta_a\ge2$

\eqref{eq:master} 中只有指示函数依赖 $\theta$，故方程组对任意 $\theta$ 形式不变。
但单个种子无法激活 $\theta_a\ge2$ 的群体，故 \eqref{eq:cascade} 给出 $C=0$，
该层在 $\mathsf{N}$ 中整行整列为零：

- 若**所有**层都有 $\theta_a\ge2$：$\mathsf{N}=0$，无线性失稳，转入 bootstrap /
  $k$-core 型渗流，须按渗流途径处理。
- 若只有**部分**层：$\mathsf{N}$ 退化为 $\theta=1$ 子层构成的子矩阵，阈值只由这些层
  决定；$\theta\ge2$ 的层在阈值处完全不起作用，却在超临界区放大终态规模。
  **这是多层特有的结论**：一条"高门槛渠道"不改变爆发与否，只改变爆发多大。
  **已验证**（$P=\{(2,2),(3,3)\}$、$m=(3,4)$、$\theta=(1,2)$）：加入 $\theta=2$ 的层
  后阈值由 $0.410619$ 变为 $0.410616$（相对移动 $9\times10^{-6}$，即二分精度内不变），
  而终态规模被放大 —— $\lambda=0.45,0.70,1.00$ 处依次为 $1.31,1.21,1.09$ 倍。
  放大在阈值上方最强、深超临界区衰减，这个非单调的形状值得作一小节并画成图。

### 3.8 三处已知的空缺（须在定稿前补上）

- **层内群体重叠**：全文只讨论层间重叠 $o_{ab}$，但同一层内两个群体共享多个节点同样
  破坏 (H1)。至少要给一个测量与一段讨论。
- **终态理论**：本大纲只做到阈值。家庭模型有成熟的最终规模自洽方程，且需要的是
  群内最终规模的**整个分布**而非仅其均值 $C$。这是审稿人一定会问的，应在第 3 节补一小节
  或明确声明留待后续。
- **爆发概率**：阈值由 $\rho(\mathsf{N})$ 决定，与种子无关；但**爆发概率**由多型分支
  过程的消亡不动点给出，依赖种子的型分布，与终态规模是两个不同的量。有向工作里正是
  这一区分给出了最有信息量的序参量。至少要写明这一点，最好给出不动点方程。


---

## 4. 数值校验

仿真为 \eqref{eq:rate} 所定义连续时间 Markov 过程的精确 Gillespie 抽样，$\mu=1$。

### 4.1 内部自洽（无需仿真）
- \eqref{eq:master} 求和是否给出 \eqref{eq:phidot}：残差随积分步长的收敛阶
- 恒等式 $x^{(a)}_{m_a-1,0,0}=\bigl[(1-\varepsilon)\psi_a(\boldsymbol{\Phi})\bigr]^{m_a-1}$
  （全体成员易感 $\Rightarrow$ 群体从未激活 $\Rightarrow$ 必未传递）
- **线性化的交叉验证**：直接线性化 \eqref{eq:master} 得到的增长矩阵，应与
  \eqref{eq:N} 逐元相同。这是本文最重要的一致性检验，因为两者路径完全不同
  （前者是 ODE 的 Jacobian，后者是分支过程计数）。
- $C(m,\lambda)$ 的两种算法：递推 \eqref{eq:cascade} vs 群体内直接 Gillespie
- 保联合度重排下 $\rho(\mathsf{N})$ 不变

### 4.2 含时复现与精度层级（**图 2**）
$M=2$、$m_1\ne m_2$ 的完整 $S(t),I(t)$：仿真 vs 群体闭合 vs 度均场。
偏差随 $N$ 的标度：闭合趋零（已验，见 3.4 的数值证据），均场应趋非零平台。
**报告纪律**：每个点必须同时给出蒙特卡洛噪声底，否则"误差重新上升"这类假象无法与
真实的闭合失效区分——初次测量正是在 $N=4000$ 处出现过这样的假上升。

**已覆盖的构型**（`tools/audit_theory.py`、`tools/scaling_asym.py`）：反相关度分布、
$\theta=(1,2)$、三层构型在 $N=4000$ 下均落在噪声内；非对称速率与群体规模
（$m=(3,5)$、$w=(1,0.4)$）在 $N=4000$ 处高出噪声 $2.7$ 倍，但标度到 $N=8000$ 即落入
噪声（$6.82\to0.39$，单位 $10^{-3}$），确认是有限尺寸效应而非闭合失效——群体越大，
同样 $N$ 下短回路越密。正文的规模扫描因而**必须按最大的 $m_a$ 选取 $N$ 的下限**。
阈值邻域须单独做有限尺寸标度，不能直接引用有限 $N$ 的值。

### 4.3 阈值的独立测量（**图 3**）
亚临界区单种子簇的平均总规模为 $\mathbf{1}^{\mathsf{T}}(\mathsf{I}-\mathsf{N})^{-1}\mathbf{v}_0$，
$\mathbf{v}_0$ 为种子的型分布。**注意不能写成 $1/(1-\rho)$**：那只在单型（$M=1$）或
$\mathbf{v}_0$ 恰为 Perron 向量时成立，$M\ge2$ 时比例常数依赖于 $\mathbf{v}_0$。
可用的是它在 $\rho\to1^-$ 时按 $(1-\rho)^{-1}$ 发散，故 $\varepsilon/R(\infty)$ 随 $\lambda$
线性趋零、零点即 $\lambda_c$。由此外推测 $\lambda_c$，与 $\rho(\mathsf{N})=1$ 比较
（给 $\sigma$ 数）。**全程不使用 $\mathsf{N}$**，是真正独立的测量。

### 4.4 三处均场偏差的逐项量化（**图 4**）
按命题 3.4 分别关闭三项（用 $\lambda$ 换 $T$、加/去 $-\delta_{ab}$、用 $(m-1)T$ 换 $C$），
在 $(m,\lambda)$ 平面上画出净偏差的符号与量级。

---

## 5. 结构生成与扫描工具

- 多重超图位形模型：按联合度序列 $\{\mathbf{k}(v)\}$ 配对各层成员槽，修复重复超边
- $\rho_{ab}$ 的调节：两条边际度序列固定下重排层度配对
- $o_{ab}$ 的调节：保联合度的定向重连
- 零模型基线与可动范围实测（见 2.3 的开放项）

---

## 6. 结构效应

1. **层间相关扫描** $\lambda_c(\rho_{12})$：把
   $\langle k^{(1)}k^{(2)}\rangle=\langle k^{(1)}\rangle\langle k^{(2)}\rangle+\rho_{12}\sigma_1\sigma_2$
   代入 \eqref{eq:N}，理论与仿真对照。协议与图可从有向工作直接复用。
2. **协同区相图**（主打结果）：两层各自亚临界，在 $(\lambda_1,\lambda_2)$ 平面上画出
   $\rho(\mathsf{N})=1$ 的曲线，标出协同区，用仿真验证边界。
3. **渠道不对称**：固定 $\sum_a w_a$，把权重在两层间转移，找最劣配置。
4. **群体规模不对称**：$m_1\ne m_2$ 时"少数大群体"与"多数小群体"哪个更危险——
   由 $C(m,\lambda)$ 的超线性增长，预期大群体更危险，需定量。
5. **层间群体重叠**：\eqref{eq:N} 只读取联合度分布，**看不见 $o_{ab}$**。
   与有向工作同款的可证伪推论：测得阈值随 $o_{ab}$ 移动，其幅度即为树状闭合失效的
   定量刻度。

---

## 7. 不可约推广：跨层剂量核

风险率取 $\beta f\bigl(\sum_a n_{e_a}\bigr)$，命题 2.1 失效。给出：
- 闭合如何修改：节点侧不再是各层独立乘积，需跨层联合状态
- 与可约情形的定量差异；若可忽略，就如实说"多层的效应主要来自结构而非动力学"

---

## 8. 真实数据

同一人群的多渠道群体数据。测量真实的 $\rho_{ab}$、$o_{ab}$ 与 $P(\mathbf{k})$，
代入 \eqref{eq:N} 预测阈值，与在真实结构上直接仿真比较。

---

## 9. 讨论与结论

- 命题 2.1 的意义：多层模型的价值在结构而非动力学，这对整个多层传播文献都成立。
- $\rho(\mathsf{N})=1$ 作为统一语言：单层、多层、成对、高阶都是其特例。
- $C(m,\lambda)$ 是高阶传播特有的量，成对极限下退化为传递概率。
- 局限：$\theta\ge2$ 的多层 bootstrap 未解；跨层剂量核只给了框架。

## 附录

A. 异质基数 · B. 非指数康复期 · C. \eqref{eq:master} 的完整推导 ·
D. $C(m,\lambda)$ 的性质与渐近 · E. 生成器与采样器 · F. 数值检验清单

---

## 附：可直接复用的资产（不属于论文正文）

| 资产 | 复用程度 |
|---|---|
| Gillespie 精确采样器 | 加层标签，逻辑不变 |
| 位形模型生成器 | 改为按联合度序列配对 |
| $r_{io}$ 扫描协议与图 | 换成 $\rho_{ab}$，协议原样 |
| 审计脚本框架 | 逐条可移植 |
| 绘图脚本与两道排版守卫 | 直接可用 |
| 构建工具链 | 直接可用 |

**不能复用**：有向工作的 (H2)（此处换为更弱的 (H2′)）；"没有余度相减"的论证
（此处相减回来了）；$R_0=\tau\kappa T$ 的形式（此处 $(m-1)T$ 换成 $C(m,\lambda)$）。

---

## 参考文献

1. F. Ball, D. Mollison, and G. Scalia-Tomba, *Epidemics with two levels of mixing*, The Annals of Applied Probability, 7 (1997), https://doi.org/10.1214/aoap/1034625252.
2. F. Ball, D. Sirl, and P. Trapman, *Threshold behaviour and final outcome of an epidemic on a random network with household structure*, Advances in Applied Probability, 41 (2009), pp. 765–796, https://doi.org/10.1239/aap/1253281063.
3. R. C. Barnard, I. Z. Kiss, L. Berthouze, and J. C. Miller, *Edge-Based Compartmental Modelling of an SIR Epidemic on a Dual-Layer Static–Dynamic Multiplex Network with Tunable Clustering*, Bulletin of Mathematical Biology, 80 (2018), pp. 2698–2733, https://doi.org/10.1007/s11538-018-0484-5.
4. G. Bianconi, *Epidemic spreading and bond percolation on multilayer networks*, Journal of Statistical Mechanics: Theory and Experiment, 2017 (2017), p. 034001, https://doi.org/10.1088/1742-5468/aa5fd8.
5. J. Fan, Q. Yin, C. Xia, and M. Perc, *Epidemics on multilayer simplicial complexes*, Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, 478 (2022), https://doi.org/10.1098/rspa.2022.0059.
6. L. Pellis, F. Ball, and P. Trapman, *Reproduction numbers for epidemic models with households and other social structures. I. Definition and calculation of R0*, Mathematical Biosciences, 235 (2012), pp. 85–97, https://doi.org/10.1016/j.mbs.2011.10.009.
7. L. Pellis, F. Ball, S. Bansal, K. Eames, T. House, V. Isham, and P. Trapman, *Eight challenges for network epidemic models*, Epidemics, 10 (2015), pp. 58–62, https://doi.org/10.1016/j.epidem.2014.07.003.
8. G. St-Onge, V. Thibeault, A. Allard, L. J. Dubé, and L. Hébert-Dufresne, *Master equation analysis of mesoscopic localization in contagion dynamics on higher-order networks*, Physical Review E, 103 (2021), https://doi.org/10.1103/physreve.103.032301.
9. H. Sun and G. Bianconi, *Higher-order percolation processes on multiplex hypergraphs*, Physical Review E, 104 (2021), https://doi.org/10.1103/physreve.104.034306.
10. J. Wan, G. Ichinose, M. Small, H. Sayama, Y. Moreno, and C. Cheng, *Multilayer networks with higher-order interaction reveal the impact of collective behavior on epidemic dynamics*, Chaos, Solitons &amp; Fractals, 164 (2022), p. 112735, https://doi.org/10.1016/j.chaos.2022.112735.
11. D. Zhao, L. Li, H. Peng, Q. Luo, and Y. Yang, *Multiple routes transmitted epidemics on multiplex networks*, Physics Letters A, 378 (2014), pp. 770–776, https://doi.org/10.1016/j.physleta.2014.01.014.
