# 摘要

现实中的高阶交互往往是有向的：施加影响的一方与接受影响的一方并不对称。我们在有向超图上建立 SIR 模型并求解其简单传播分支。每条超边分为不相交的尾集与头集，尾集中被感染的成员数达到阈值 $\theta$ 时超边激活，以速率 $\beta$ 向头集中每个易感节点传递；感染者以速率 $\mu$ 康复，控制参数 $\lambda=\beta/\mu$。$\theta=1$ 时的困难在于一条超边只有一个激活时钟却由 $\tau$ 个尾成员共同驱动：激活时长是各成员感染期之并，闭合变量不能逐成员相乘。以尾集中易感、感染、康复的成员计数为状态，可把"尚未传递"写成一条以传递事件为消灭机制的 Markov 链，得到 $\binom{\tau+2}{2}+1$ 维、与系统规模无关的闭方程组，其解给出完整的含时演化；终态在 $\tau=1$ 时另有闭式，$\tau\ge2$ 时须由方程组积分。线性化给出爆发阈值 $\lambda_c=1/(\tau\kappa-1)$，其中 $\kappa=\langle k_{\rm in}k_{\rm out}\rangle/\langle k_{\rm out}\rangle$；$\tau=\eta=1$ 时它退化为有向随机图上 SIR 的已知结果。阈值经 $\kappa$ 依赖于入出度相关，却与超边重叠无关——后者是树状闭合的一条可证伪推论。与精确 Gillespie 仿真比较，闭合对 $S(t)$ 的偏差随系统规模趋零，而同一模型的度均场把阈值压低到 $0.75\lambda_c$，在真实的亚临界区就给出宏观爆发，其偏差收敛到非零平台；阈值本身由亚临界终态外推独立测得，与解析值相差 $1.4\sigma$。模型另给出两个与直觉相反的结构结论：均匀单种子下的平均爆发规模几乎不能探测方向对称性破缺，须改用爆发概率与条件爆发规模；方向翻转下的奇宇称是破缺的必要而非充分条件。$\theta\ge2$ 时激活不再是单条超边的独立事件，模型转入 bootstrap 与 $k$-core 型渗流，另行处理。

**关键词.** 有向超图，高阶交互，SIR 传播，边基房室模型，空腔方法，生成函数，渗流，爆发阈值，方向对称性破缺

**MSC 分类号.** 05C65, 05C20, 05C80, 92D30, 60J28, 82B43, 37N25

# 1. 引言

许多传染并非沿边逐个发生，而是在一个群体内集体发生 [3, 4]，这类高阶交互在真实数据中普遍存在 [6]。但群体交互往往是有向的：施加影响的一方与接受影响的一方并不对称。代谢反应把一组底物转化为一组产物 [21]，比特币交易与引用数据同样是一组主体指向另一组 [27]。有向超图正是为刻画这种不对称而设的形式化——每条超边分为施加影响的尾集与接受影响的头集（二者不相交），这一结构在算法与运筹领域已使用三十余年 [14]。

已有的理论工具恰好绕开了这个交叉点。成对网络上的 SIR 已由两条彼此独立的路径解决：渗流给出阈值与终态 [16, 34, 35]——临界占据概率为 $\langle k\rangle/(\langle k^2\rangle-\langle k\rangle)$，二阶矩发散时趋零，无标度网络上的 SIS 模型给出过同一结论 [36]；方程给出峰值与时间演化——Volz 的低维闭合 [42] 经 Miller 等整理为边基房室模型（EBCM），以"沿一条边尚未传来感染"的概率为核心变量 [32]；二者可相互印证 [37]。方向性作为成对网络的属性也已处理成熟：有向随机图各巨分支的规模由入出度联合分布决定 [10]，阈值与分支规模可解析求出 [8, 31]，且再生数经 $\langle k_{\rm in}k_{\rm out}\rangle$ 依赖于入出度的相关 [2]。

高阶交互一侧同样成果丰富：单纯复形与超图上的传播可出现不连续相变与双稳 [12, 17, 38]，结构的影响也已刻画——链接度异质性抑制爆炸式相变 [24]，超核分解给出局域化的中心性 [30]，高阶连通分支的巨分支是单种子全局爆发的必要条件 [20]，最佳种子群体已被刻画 [39]，而超边重叠决定转变是否爆炸式且方向与直觉相反 [23, 28, 29]，超图渗流理论亦已建立 [7, 13]。两条线索各自成熟，却始终没有相交：处理方向性的工作都停留在成对网络，处理高阶交互的工作则都假设超边无向。

有向超图上，结构刻画的工具近年已迅速成型——真实数据的微观组织与高阶互惠性 [19, 27]、重叠测度 [25]、位形模型与保度随机化 [9, 33]、保入出度序列的均匀采样 [1, 22]。动力学一侧则少得多，据我们所知以 Li 等的工作为主：他们在有向超图上建立社会传播模型，发现双稳区间随有向强度减弱而收缩 [26]。该工作采用均场闭合，方向性由单一标量参数调节。均场忽略邻域相关，误差在阈值邻域最大；成对网络上正是 EBCM 补上了这一层 [32, 42]。因此，有向超图上的传播阈值如何依赖于结构，据我们所知尚无超出均场的含时理论，也无可供解析扫描的渗流理论。

本文建立并求解有向超图上的 SIR 模型。第 2 节定义模型：尾集中被感染的成员数达到群体阈值 $\theta$ 时超边激活，以速率 $\beta$ 向头集中每个易感节点传递，感染者以速率 $\mu$ 康复，控制参数为 $\lambda=\beta/\mu$；非指数康复期会显著改变成对网络的流行阈值 [41]，其推广见附录。第 3 节给出 $\theta=1$ 的并发型 EBCM 方程。技术困难在于一条超边只有一个激活时钟却由 $\tau$ 个尾成员共同驱动——激活时长是各成员感染期之并，因而闭合变量不能逐成员相乘。以尾集中易感、感染、康复的成员计数为状态，可把"尚未传递"写成一条以传递事件为消灭机制的 Markov 链——等价于给该链添一个"已传递"的吸收态——得到 $\binom{\tau+2}{2}+1$ 维、与系统规模无关的闭方程组。

爆发条件由线性化给出，$\lambda_c=1/(\tau\kappa-1)$，$\kappa=\langle k_{\rm in}k_{\rm out}\rangle/\langle k_{\rm out}\rangle$；$\tau=\eta=1$ 时退化为有向随机图的已知阈值 [2, 31]。阈值经 $\kappa$ 依赖于入出度相关 $r_{io}$，却与超边重叠无关——这不是疏漏，而是树状闭合的一条可证伪推论。$r_{io}$ 的作用方向依传播机制而异：有向网络的阈值型级联中，正的入出度相关反而提高系统稳健性 [43]，机制上的原因见第 3 节。与精确 Gillespie 仿真 [15] 比较，闭合的偏差随系统规模趋零，均场的不然：它把阈值压低四分之一，在真实的亚临界区就给出宏观爆发。

第 4 节转向渗流一侧，把成对网络上的 SIR–渗流映射及其半有向修正 [18, 35] 与有向渗流的生成函数机制 [8, 10] 推广到有向超图。两条路径互不依赖，因而可互为验证；我们进一步证明二者的自洽方程在 $\mu\to0$ 极限下同构。$\theta\ge2$ 时激活不再是单条超边的独立事件，模型转入 bootstrap 与 $k$-core 型渗流，其相变结构不同于普通渗流 [5, 11, 40]，亦在该节处理。

为使结构依赖可被受控地检验，第 2 节给出有向超边重叠的可计算定义：对每一对有序超边构造 $2\times2$ 列联矩阵，四个通道分别对应共享节点同为尾成员、同为头成员、或一为头一为尾。四通道之和精确等于无向重叠，因而结果可与无向文献直接对照 [23, 29]；互惠方向另采用文献中的强互惠判据 [27]。入出度相关 $r_{io}$ 与重叠 $\alpha$ 作用在生成流程的不同环节，因而可独立扫描；零模型基线取自已有的均匀采样器 [1, 22]。

方向对称性破缺的分析（第 II、VI 节）给出两个与直觉相反的结果。其一，均匀单种子下的平均爆发规模几乎不能用来探测破缺：可达有序对的总数在方向翻转下近乎守恒，翻转真正改变的是这些可达对如何分配——大爆发的概率由巨入分支决定、终态规模由巨出分支决定，而翻转恰好交换二者 [2, 18]。我们因而以爆发概率与条件爆发规模为序参量。其二，翻转奇是破缺的必要条件却远非充分：双度序列的不对称驱动的破缺是数十个标准差的效应，而同样翻转奇的重叠极性 $\Delta\alpha$ 效应始终不可分辨，我们只能给出上界。破缺因而由分支结构的不对称驱动，而非由重叠的极性驱动。

其余各节依次为：第 5 节的保度生成与定向重连工具，第 7 节在真实有向超图数据上的检验，第 8 节的讨论与结论。

---

## 参考文献

1. M. Abuissa, M. Riondato, and E. Upfal, *DiNgHy: null models for non-degenerate directed hypergraphs*, Data Mining and Knowledge Discovery, 40 (2026), https://doi.org/10.1007/s10618-026-01209-8.
2. A. Allard, C. Moore, S. V. Scarpino, B. M. Althouse, and L. Hébert-Dufresne, *The Role of Directionality, Heterogeneity, and Correlations in Epidemic Risk and Spread*, SIAM Review, 65 (2023), pp. 471–492, https://doi.org/10.1137/20M1383811.
3. F. Battiston, G. Cencetti, I. Iacopini, V. Latora, M. Lucas, A. Patania, J. G. Young, and G. Petri, *Networks beyond pairwise interactions: Structure and dynamics*, Physics Reports, 874 (2020), pp. 1–92, https://doi.org/10.1016/j.physrep.2020.05.004.
4. F. Battiston, E. Amico, A. Barrat, G. Bianconi, G. Ferraz de Arruda, B. Franceschiello, I. Iacopini, S. Kéfi, and et al., *The physics of higher-order interactions in complex systems*, Nature Physics, 17 (2021), pp. 1093–1098, https://doi.org/10.1038/s41567-021-01371-4.
5. G. J. Baxter, S. N. Dorogovtsev, A. V. Goltsev, and J. F. F. Mendes, *Bootstrap percolation on complex networks*, Physical Review E, 82 (2010), p. 011103, https://doi.org/10.1103/PhysRevE.82.011103.
6. A. R. Benson, R. Abebe, M. T. Schaub, A. Jadbabaie, and J. Kleinberg, *Simplicial closure and higher-order link prediction*, Proceedings of the National Academy of Sciences, 115 (2018), pp. E11221–E11230, https://doi.org/10.1073/pnas.1800683115.
7. G. Bianconi and S. N. Dorogovtsev, *Theory of percolation on hypergraphs*, Physical Review E, 109 (2024), p. 014306, https://doi.org/10.1103/PhysRevE.109.014306.
8. M. Boguñá and M. Á. Serrano, *Generalized percolation in random directed networks*, Physical Review E, 72 (2005), p. 016106, https://doi.org/10.1103/PhysRevE.72.016106.
9. P. S. Chodrow, *Configuration models of random hypergraphs*, Journal of Complex Networks, 8 (2020), p. cnaa018, https://doi.org/10.1093/comnet/cnaa018.
10. S. N. Dorogovtsev, J. F. F. Mendes, and A. N. Samukhin, *Giant strongly connected component of directed networks*, Physical Review E, 64 (2001), p. 025101, https://doi.org/10.1103/PhysRevE.64.025101.
11. S. N. Dorogovtsev, A. V. Goltsev, and J. F. F. Mendes, *k-core organization of complex networks*, Physical Review Letters, 96 (2006), p. 040601, https://doi.org/10.1103/PhysRevLett.96.040601.
12. G. Ferraz de Arruda, G. Petri, and Y. Moreno, *Social contagion models on hypergraphs*, Physical Review Research, 2 (2020), p. 023032, https://doi.org/10.1103/PhysRevResearch.2.023032.
13. G. Ferraz de Arruda, A. Aleta, and Y. Moreno, *Contagion dynamics on higher-order networks*, Nature Reviews Physics, 6 (2024), pp. 468–482, https://doi.org/10.1038/s42254-024-00733-0.
14. G. Gallo, G. Longo, S. Pallottino, and S. Nguyen, *Directed hypergraphs and applications*, Discrete Applied Mathematics, 42 (1993), pp. 177–201, https://doi.org/10.1016/0166-218X(93)90045-P.
15. D. T. Gillespie, *Exact stochastic simulation of coupled chemical reactions*, The Journal of Physical Chemistry, 81 (1977), pp. 2340–2361, https://doi.org/10.1021/j100540a008.
16. P. Grassberger, *On the critical behavior of the general epidemic process and dynamical percolation*, Mathematical Biosciences, 63 (1983), pp. 157–172, https://doi.org/10.1016/0025-5564(82)90036-0.
17. I. Iacopini, G. Petri, A. Barrat, and V. Latora, *Simplicial models of social contagion*, Nature Communications, 10 (2019), p. 2485, https://doi.org/10.1038/s41467-019-10431-6.
18. E. Kenah and J. M. Robins, *Second look at the spread of epidemics on networks*, Physical Review E, 76 (2007), p. 036113, https://doi.org/10.1103/PhysRevE.76.036113.
19. S. Kim, M. Choe, J. Yoo, and K. Shin, *Reciprocity in directed hypergraphs: measures, findings, and generators*, Data Mining and Knowledge Discovery (2023), https://doi.org/10.1007/s10618-023-00955-3.
20. J. H. Kim and K. I. Goh, *Higher-Order Components Dictate Higher-Order Contagion Dynamics in Hypergraphs*, Physical Review Letters, 132 (2024), p. 087401, https://doi.org/10.1103/PhysRevLett.132.087401.
21. S. Klamt, U. U. Haus, and F. Theis, *Hypergraphs and cellular networks*, PLoS Computational Biology, 5 (2009), p. e1000385, https://doi.org/10.1371/journal.pcbi.1000385.
22. Y. J. Kraakman and C. Stegehuis, *Uniformly sampling random directed hypergraphs with fixed degrees*, Discrete Mathematics, 349 (2026), p. 114961, https://doi.org/10.1016/j.disc.2025.114961.
23. S. Lamata-Otín, F. Malizia, V. Latora, M. Frasca, and J. Gómez-Gardeñes, *Hyperedge overlap drives synchronizability of systems with higher-order interactions*, Physical Review E, 111 (2025), p. 034302, https://doi.org/10.1103/PhysRevE.111.034302.
24. N. W. Landry and J. G. Restrepo, *The effect of heterogeneity on hypergraph contagion models*, Chaos, 30 (2020), p. 103117, https://doi.org/10.1063/5.0020034.
25. G. Lee, M. Choe, and K. Shin, *How Do Hyperedges Overlap in Real-World Hypergraphs? - Patterns, Measures, and Generators*, Proceedings of the Web Conference 2021 (2021), https://doi.org/10.1145/3442381.3450010.
26. J. Li, X. Wu, J. Lü, and L. Lei, *Enhancing predictive accuracy in social contagion dynamics via directed hypergraph structures*, Communications Physics, 7 (2024), https://doi.org/10.1038/s42005-024-01614-9.
27. Q. F. Lotito, A. Vendramini, A. Montresor, and F. Battiston, *The microscale organization of directed hypergraphs*, Communications Physics, 9 (2026), https://doi.org/10.1038/s42005-025-02472-9.
28. F. Malizia, A. Guzmán, I. Iacopini, and I. Z. Kiss, *Disentangling the Role of Heterogeneity and Hyperedge Overlap in Explosive Contagion on Higher-Order Networks*, Physical Review Letters (2025), https://doi.org/10.1103/z3d5-94zb.
29. F. Malizia, S. Lamata-Otín, M. Frasca, V. Latora, and J. Gómez-Gardeñes, *Hyperedge overlap drives explosive transitions in systems with higher-order interactions*, Nature Communications, 16 (2025), https://doi.org/10.1038/s41467-024-55506-1.
30. M. Mancastroppa, I. Iacopini, G. Petri, and A. Barrat, *Hyper-cores promote localization and efficient seeding in higher-order processes*, Nature Communications, 14 (2023), p. 6223, https://doi.org/10.1038/s41467-023-41887-2.
31. L. A. Meyers, M. E. J. Newman, and B. Pourbohloul, *Predicting epidemics on directed contact networks*, Journal of Theoretical Biology, 240 (2006), pp. 400–418, https://doi.org/10.1016/j.jtbi.2005.10.004.
32. J. C. Miller, A. C. Slim, and E. M. Volz, *Edge-based compartmental modelling for infectious disease spread*, Journal of the Royal Society Interface, 9 (2012), pp. 890–906, https://doi.org/10.1098/rsif.2011.0403.
33. K. Nakajima, K. Shudo, and N. Masuda, *Randomizing Hypergraphs Preserving Degree Correlation and Local Clustering*, IEEE Transactions on Network Science and Engineering, 9 (2022), pp. 1139–1153, https://doi.org/10.1109/TNSE.2021.3133380.
34. M. E. J. Newman, S. H. Strogatz, and D. J. Watts, *Random graphs with arbitrary degree distributions and their applications*, Physical Review E, 64 (2001), p. 026118, https://doi.org/10.1103/PhysRevE.64.026118.
35. M. E. J. Newman, *Spread of epidemic disease on networks*, Physical Review E, 66 (2002), p. 016128, https://doi.org/10.1103/PhysRevE.66.016128.
36. R. Pastor-Satorras and A. Vespignani, *Epidemic spreading in scale-free networks*, Physical Review Letters, 86 (2001), pp. 3200–3203, https://doi.org/10.1103/PhysRevLett.86.3200.
37. R. Pastor-Satorras, C. Castellano, P. Van Mieghem, and A. Vespignani, *Epidemic processes in complex networks*, Reviews of Modern Physics, 87 (2015), pp. 925–979, https://doi.org/10.1103/RevModPhys.87.925.
38. G. St-Onge, H. Sun, A. Allard, L. Hébert-Dufresne, and G. Bianconi, *Universal nonlinear infection kernel from heterogeneous exposure on higher-order networks*, Physical Review Letters, 127 (2021), p. 158301, https://doi.org/10.1103/PhysRevLett.127.158301.
39. G. St-Onge, I. Iacopini, V. Latora, A. Barrat, G. Petri, A. Allard, and L. Hébert-Dufresne, *Influential groups for seeding and sustaining nonlinear contagion in heterogeneous hypergraphs*, Communications Physics, 5 (2022), p. 25, https://doi.org/10.1038/s42005-021-00788-w.
40. H. Sun and G. Bianconi, *Higher-order percolation processes on multiplex hypergraphs*, Physical Review E, 104 (2021), p. 034306, https://doi.org/10.1103/PhysRevE.104.034306.
41. P. Van Mieghem and R. van de Bovenkamp, *Non-Markovian infection spread dramatically alters the susceptible-infected-susceptible epidemic threshold in networks*, Physical Review Letters, 110 (2013), p. 108701, https://doi.org/10.1103/PhysRevLett.110.108701.
42. E. Volz, *SIR dynamics in random networks with heterogeneous connectivity*, Journal of Mathematical Biology, 56 (2008), pp. 293–310, https://doi.org/10.1007/s00285-007-0116-4.
43. X. J. Xu, J. Y. Li, X. Fu, and L. J. Zhang, *Impact of directionality and correlation on contagion*, Scientific Reports, 8 (2018), https://doi.org/10.1038/s41598-018-22508-1.
