# Fig. 1 — Chaos/EBCM multiplex-hypergraph paper

Redrawn schematic for *"Edge-based compartmental theory for susceptible–infected–recovered
spreading on multiplex hypergraphs"* (Wang, You & Pan), in the visual idiom of
**Chaos, Solitons & Fractals** / Elsevier: Arial-metric sans-serif throughout, filled
hyperedge blobs, SIR-coloured nodes, bold `(a)`–`(d)` panel labels, one shared legend.

## Output

| file | width | use |
|---|---|---|
| `fig1_csf_190mm.pdf` | 190.0 × 158.0 mm | Elsevier / CSF double column |
| `fig1_csf_178mm.pdf` | 177.8 × 147.9 mm | REVTeX 4 two-column `\textwidth` (7 in), AIP *Chaos* |

Both are **fully vector** (no embedded raster), with selectable text and all base
lettering at **7.0–9.0 pt at final size** — the two widths differ in geometry only,
the type is identical, so neither needs rescaling on `\includegraphics`.

## Include

```latex
\begin{figure*}[t]
  \centering
  \includegraphics[width=\textwidth]{fig1_csf_178mm.pdf}  % 190mm variant for elsarticle
  \caption{\label{fig:framework} <caption below> }
\end{figure*}
```

## Caption

> **FIG. 1.** Structure, edge-based closure, and the two factors of the epidemic
> threshold. **(a)** A two-layer multiplex hypergraph. Layer $a$ is $m_a$-uniform
> (solid, $m_a=3$), layer $b$ is $m_b$-uniform (dashed, $m_b=4$), and the test node
> $u$ has layer-specific hyperdegrees $k^{(a)}(u)=2$, $k^{(b)}(u)=1$. The joint
> hyperdegree distribution $P(\mathbf{k})$ enters the theory through the generating
> functions $\Psi$, for a uniformly sampled node, and $\psi_a$, for a node reached
> through layer $a$ with its arrival hyperedge excluded.
> **(b)** Edge-based closure at finite prevalence. Transmission out of $u$ is
> suppressed, which makes its incident hyperedges conditionally independent on a
> locally tree-like incidence graph. For one representative layer-$a$ hyperedge $e$
> the theory retains the joint SIR composition $x^{(a)}_{sir}$ of the $m_a-1$ other
> members; those members are additionally infected from their remaining incident
> hyperedges at rate $h_a$, obtained from $\psi_a(\Phi)$. The cycle
> $x^{(a)}_{sir}\to A_a\to\Phi_a\to\psi_a(\Phi)\to h_a\to x^{(a)}_{sir}$ closes the
> dynamics and determines $S(t)$, $I(t)$, $R(t)$ and $R(\infty)$.
> **(c)** Rare-infection limit, within-hyperedge factor. A single seed in an
> otherwise susceptible hyperedge infects $(m-1)T$ members directly, but secondary
> infections lengthen the interval over which $e$ still contains an infectious node,
> so the expected total satisfies $C(m,\lambda,\theta)\geq(m-1)T$, with equality only
> at $m=2$. States are labelled by $(i,s)$ as in Eq.~(20).
> **(d)** Rare-infection limit, network factor. A newly infected node is a member of
> the hyperedge through which infection arrived — excluded from its offspring, giving
> the $-\delta_{ab}$ of Eq.~(23) — and of $B_{ab}$ further layer-$b$ hyperedges, each
> seeding $C_b$ new infections on average. The product $K_{ab}=B_{ab}C_b$ is the
> next-generation matrix, and $\rho(K)=1$ fixes the threshold $\lambda_c$.

## What changed relative to the previous Fig. 1

**Content**

1. The previous caption described panel (b) as the finite-prevalence edge-based
   construction yielding $S(t),I(t),R(t),R(\infty)$, but the panel actually showed the
   rare-infection quantity $C$ versus $(m-1)T$. Nothing in the old figure depicted
   $\Phi_a$ (Eq. 7), $x^{(a)}_{sir}$ (Eq. 9) or $h_a$ (Eq. 13) — the paper's headline
   contribution. Panel (b) is new and draws that closure explicitly, so the existing
   sentence in Sec. II A ("Figure 1 summarizes this construction and the subsequent
   rare-infection limit") is now accurate as written.
2. Two sub-diagrams in the old panel (b) — `direct` and `(i,s)=(1,3)` — were each
   missing a node: the hyperedge arc had a gap and a transmission arrow dangled with
   no arrowhead and no target, while the neighbouring `(2,2)` and `(3,1)` were complete
   4-node rings. All sub-diagrams here carry the full $m=4$ nodes, and `direct only`
   shows all three arrows.
3. The old strict inequality $C(m,\lambda,\theta)>(m-1)T$ is false at $m=2$, the very
   pairwise limit Eq. (22) discusses. Relaxed to $\geq$ with "equality iff $m=2$";
   verified numerically from the Eq. (20) recursion (see below).
4. In the old panel (c) the newly infected node floated between hyperedges, which reads
   as a pairwise link between groups — the projection the paper argues against. It is
   now drawn on the intersection of the three hyperedges it belongs to, matching how
   $u$ is drawn in panel (a).
5. Dark fill meant "test node $u$" in the old panel (a) but "infected" in (b) and (c).
   $u$ now has its own encoding (white, heavy ring) and the figure has a legend.
   A recovered node appears in panel (b), so all three SIR compartments are shown.
6. Panel (a) now shows $P(\mathbf{k})$, $\Psi$ and $\psi_a$, which the old caption
   claimed but the drawing omitted.
7. Layer indices unified to $a$/$b$ across all panels; each closure step carries its
   equation number.

**Production**

8. Old Fig. 1 was a **JPEG raster, 2595 × 864 px placed at 6.774 × 2.256 in = 383 dpi** —
   below AIP's 600 dpi for line art, lossy on black-on-white text, and the only raster
   figure in a paper whose Figs. 2–4 are vector. Now vector PDF.
9. Old in-figure lettering measured ~5 pt against a 9 pt caption (ink heights: caption
   "FIG. 1." 6.58 pt, "(i, s) = (1, 3)" 3.76 pt, "time with infectious nodes" 3.38 pt),
   with a ~3× internal size range. Now a single 7.0–9.0 pt scale.
10. Panel labels changed from `a. multiplex hypergraph` to `(a)`, matching Figs. 2–4.
11. Palette is grayscale- and deuteranopia-safe: node luminances 0.766 / 0.538 / 0.304
    for S / I / R form an even ladder, so the three compartments separate in
    black-and-white print as well as in colour.

## Text change still needed in the manuscript

The vocabulary of panel (c) — "direct", "secondary", $(m-1)T$, and the inequality —
appears nowhere in the body text. Suggested insertion after Eq. (22), Sec. II C:

> It is useful to compare $C(m,\lambda,\theta)$ with the number of infections the seed
> generates on its own. If secondary transmission inside the hyperedge is suppressed,
> the seed remains infectious for a mean time $1/\mu$ and each of the $m-1$ other
> members is infected with probability $T$, giving $(m-1)T$. Secondary infections both
> add infectious members and extend the interval over which the hyperedge satisfies its
> transmission condition, so $C(m,\lambda,1)\geq(m-1)T$, with equality only in the
> pairwise case $m=2$. The gap widens with hyperedge size: at $\lambda=0.5$ the ratio
> $C/[(m-1)T]$ is $1.20$ for $m=4$ and $1.60$ for $m=8$. This is the mechanism by which
> larger hyperedges lower the threshold at fixed hyperdegree distribution, and it is the
> quantity isolated in Fig.~2(b).

Values from the Eq. (20) recursion:

| $m$ | $\lambda$ | $C$ | $(m-1)T$ | ratio |
|---|---|---|---|---|
| 2 | 0.5 | 0.3333 | 0.3333 | 1.000 |
| 4 | 0.5 | 1.2019 | 1.0000 | 1.202 |
| 8 | 0.5 | 3.7367 | 2.3333 | 1.601 |

A second citation of Fig. 1 near Eq. (24) (e.g. "cf. Fig. 1(c),(d)") would also help;
the old figure was cited exactly once in the whole paper.

## Rebuild

```bash
pip install matplotlib                       # needs a Liberation Sans (Arial-metric) font
FIG_WIDTH_MM=190   python3 fig1.py           # Elsevier / CSF
FIG_WIDTH_MM=177.8 python3 fig1.py           # REVTeX / AIP Chaos
```

`style.py` holds the palette, type scale and the hyperedge-blob geometry (rounded
convex offset = Minkowski sum of the members' convex hull with a disk). Layout
coordinates in `fig1.py` are millimetres on a fixed 190 × 158 canvas; `FIG_WIDTH_MM`
rescales the geometry while leaving type at its point size.
