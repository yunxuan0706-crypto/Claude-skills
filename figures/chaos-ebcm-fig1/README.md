# Fig. 1 — Chaos/EBCM multiplex-hypergraph paper

Redrawn schematic for *"Edge-based compartmental theory for susceptible–infected–recovered
spreading on multiplex hypergraphs"* (Wang, You & Pan), in the visual idiom of
**Chaos, Solitons & Fractals** / Elsevier: Arial-metric sans-serif throughout, filled
hyperedge blobs, SIR-coloured nodes, bold `(a)`–`(d)` panel labels, one shared legend.

## Two variants

`fig1_lean_*` is the one to submit. `fig1_csf_*` carries the same drawing plus the
prose and equation numbers inline, which makes it readable with no caption at hand —
useful for a talk slide or a supplementary figure, too heavy for a journal page.

|  | text elements | characters | equation refs | size |
|---|---|---|---|---|
| `fig1_lean_190mm.pdf` | 125 | 336 | 0 | 190.0 × 138.0 mm |
| `fig1_csf_190mm.pdf`  | 251 | 1133 | 17 | 190.0 × 158.0 mm |

The lean variant keeps only labels that **name something drawn in the panel**
($u$, $e_1$, $\beta_a$, $h_a$, $x^{(a)}_{sir}$, $(i,s)$, $B_{aa}$, $K_{ab}=B_{ab}C_b$,
the closure cycle). Every sentence and every `[Eq. (n)]` moved into the caption, which
is why the lean caption below is longer than the full variant's.

Each variant is built at two widths:

| suffix | width | use |
|---|---|---|
| `_190mm` | 190.0 mm | Elsevier / CSF double column |
| `_178mm` | 177.8 mm | REVTeX 4 two-column `\textwidth` (7 in), AIP *Chaos* |

All four are **fully vector** (no embedded raster), with selectable text and all base
lettering at **7.0–9.0 pt at final size**. The two widths differ in geometry only —
the type is identical — so neither needs rescaling on `\includegraphics`.

## Palette

Elsevier's artwork guidance that bears on colour: vector formats preferred, 190 mm
double column, embedded fonts at a uniform size, and colour images that stay readable
for readers with impaired colour vision. Elsevier also asks for a usable black-and-white
version alongside the colour one, because automatic colour-to-grey conversion often
fails. **Whether colour in the printed version carries a charge is journal-specific —
check the CSF Guide for Authors.**

`style.py` ships three palettes, selected with the `PALETTE` environment variable.
`duo` is the default and the one to submit.

| `PALETTE` | hues | notes |
|---|---|---|
| `full` | 5 (blue, salmon, slate, cool + warm layer tints, brick, teal) | reads as an infographic; the warm layer-$b$ tint is redundant with the dashed outline, and the teal appears on one arrow |
| **`duo`** | **2 (blue, red) + neutrals** | **layers carried by line style alone, one neutral blob tint, coupling arrow neutral** |
| `mono` | 1 (red) + neutrals | most conservative, but grey S and grey R send the reader back to the legend |

`palette_compare.png` shows panel (d) and the legend under all three.

Colour never carries information on its own here. Node **state** is fill, node
**identity** ($u$) is a heavy ring plus a label, **layer** is line style, and every
category is named in the legend — so the figure survives greyscale printing intact.
The four node fills form an even luminance ladder, which is what makes that work:

| | $u$ | S | I | R | smallest gap |
|---|---|---|---|---|---|
| `full` | 1.000 | 0.766 | 0.538 | 0.304 | 0.228 |
| `duo`  | 1.000 | 0.743 | 0.513 | 0.276 | 0.230 |
| `mono` | 1.000 | 0.776 | 0.570 | 0.332 | 0.206 |

`check_cvd.png` renders the legend in greyscale, deuteranopia, protanopia and
tritanopia; all four node types stay distinct in each.

## Include

```latex
\begin{figure*}[t]
  \centering
  \includegraphics[width=\textwidth]{fig1_lean_178mm.pdf}  % 190mm variant for elsarticle
  \caption{\label{fig:framework} <caption below> }
\end{figure*}
```

## Caption (lean variant)

> **FIG. 1.** Structure, edge-based closure, and the two factors of the epidemic
> threshold. Node fill gives the SIR state and outline style gives the layer (legend).
> **(a)** A two-layer multiplex hypergraph: layer $a$ is $m_a$-uniform (solid,
> $m_a=3$, hyperedges $e_1,e_2$), layer $b$ is $m_b$-uniform (dashed, $m_b=4$,
> hyperedge $e_3$), and the test node $u$ has $k^{(a)}(u)=2$, $k^{(b)}(u)=1$. The joint
> hyperdegree distribution $P(\mathbf{k})$ enters through the generating functions
> $\Psi$, for a uniformly sampled node [Eq.~(3)], and $\psi_a$, for a node reached
> through layer $a$ with its arrival hyperedge excluded [Eq.~(4)].
> **(b)** Edge-based closure at finite prevalence. Transmission out of $u$ is
> suppressed ($\times$), which makes its incident hyperedges conditionally independent
> on a locally tree-like incidence graph. For one representative layer-$a$ hyperedge the
> theory retains the joint SIR composition $x^{(a)}_{sir}$ of the $m_a-1$ other members
> [Eq.~(9)]; that hyperedge transmits to $u$ at rate $\beta_a$, while the same members
> are infected from their remaining incident hyperedges at rate $h_a$ [Eq.~(13)],
> obtained from $\psi_a(\Phi)$. The cycle
> $x^{(a)}_{sir}\to A_a\to\Phi_a\to\psi_a(\Phi)\to h_a$ closes the dynamics
> [Eqs.~(11), (12), (4), (14)] and determines $S(t)$, $I(t)$, $R(t)$ and $R(\infty)$
> [Eqs.~(8), (17), (18)].
> **(c)** Rare-infection limit, within-hyperedge factor. A seed acting alone infects
> $(m-1)T$ members of an otherwise susceptible hyperedge. Secondary infections both add
> infectious members and prolong the interval over which the hyperedge still contains
> one, so the expected total obeys $C(m,\lambda,\theta)\geq(m-1)T$, with equality only
> at $m=2$. States are labelled $(i,s)$ as in Eq.~(20).
> **(d)** Rare-infection limit, network factor. A newly infected node belongs to the
> hyperedge through which infection arrived — excluded from its offspring, contributing
> the $-\delta_{ab}$ of Eq.~(23) — and to $B_{aa}$ further layer-$a$ and $B_{ab}$
> further layer-$b$ hyperedges, each seeding $C_a$ or $C_b$ new infections on average.
> The product $K_{ab}=B_{ab}C_b$ is the next-generation matrix [Eq.~(24)], and
> $\rho(K)=1$ fixes the threshold $\lambda_c$ [Eq.~(25)].

## Editable PowerPoint

`fig1_lean_editable.pptx` and `fig1_annotated_editable.pptx` are the same two figures
as **native PowerPoint shapes** — ellipses, freeform hyperedge blobs, arrows with real
arrowheads, rounded rectangles and text boxes with genuine sub/superscript runs. No
image is embedded in either deck (`0 pictures`), so every node can be dragged, every
label retyped, and every colour changed in PowerPoint.

| deck | slide | shapes |
|---|---|---|
| `fig1_lean_editable.pptx` | 190.0 × 138.0 mm | 127 |
| `fig1_annotated_editable.pptx` | 190.0 × 158.0 mm | 156 |

The slide is sized to the figure, so **File → Export → PDF gives the figure at its
correct dimensions** with no cropping or rescaling. Type is at its true 7–9 pt, which
looks small on screen — zoom in PowerPoint rather than resizing anything.

Two limits worth knowing:

- PowerPoint cannot stack a superscript over a subscript, so `x^{(a)}_{sir}` renders as
  a superscript followed by a subscript rather than the two aligned vertically. Every
  other symbol matches the PDF.
- Moving a node does not reshape the hyperedge blob around it — the blob is one
  freeform, not a container. Drag the blob's own outline to follow.

The decks are generated from the same layout code as the PDFs, not redrawn: the
drawing helpers in `style.py` record what they draw, `export_spec` writes that out as
JSON with y measured from the top, and `build_pptx.js` replays it through pptxgenjs.
Text boxes are placed from extents measured off the real matplotlib render, which is
what keeps the two formats aligned. `mathrun.py` converts the mathtext (`$\Phi_a$`,
`$x^{(a)}_{sir}$`) into PowerPoint runs with Greek characters and baseline offsets.

```bash
npm install pptxgenjs
EXPORT_SPEC=fig1_lean_spec.json FIG_WIDTH_MM=190 python3 fig1_lean.py
node build_pptx.js fig1_lean_spec.json fig1_lean_editable.pptx
```

Edit the Python and both formats follow. Edit the `.pptx` and only the deck changes.

## Rebuilding it by hand in PowerPoint

The `*_editable.pptx` decks above are editable but not *reproducible*: their hyperedges
are freeform shapes traced from a convex-hull offset, 70-odd points each, which nobody
is going to redraw by hand. `fig1_ppt_reproducible.pptx` is a third variant whose
geometry was redesigned so that **every element is one shape from PowerPoint's own
gallery**:

| element | PowerPoint shape |
|---|---|
| hyperedge | **Oval** (a circle, or an ellipse plus a rotation) |
| node | **Oval** (circle) |
| arrow, rule, exclusion cross | **Line**, arrowhead set on the line |
| $K_{ab}$ box | **Rounded Rectangle** |
| every label | **Text Box** |

Hyperedges are circles or ellipses with their member nodes placed on the boundary —
the construction the original figure used — so shared nodes fall naturally on
intersections. It contains **no freeform shapes and no Bezier arrows**: the closure
feedback loop is a three-segment elbow and the $h_a$ coupling arrow is straight.
python-pptx reports only `AUTO_SHAPE` and `TEXT_BOX`, 128 shapes, 0 pictures.

| file | what it is |
|---|---|
| `fig1_ppt_190mm.png` / `.pdf` | the target to reproduce (190 × 145 mm) |
| `fig1_ppt_reproducible.pptx` | the same figure as standard shapes — open it, copy pieces out |
| `BUILD_SHEET.md` | every shape's exact geometry in centimetres |

`BUILD_SHEET.md` is generated from the figure, not written by hand
(`python3 make_build_sheet.py fig1_ppt_spec.json`), so its numbers cannot drift from
the drawing. It gives the slide size, a style registry (node fills and diameters,
hyperedge line/fill and transparency, arrow colours, text sizes), and then per panel a
table of positions: X/Y of each shape's bounding box, width, height and rotation in
the exact units PowerPoint's Format Shape pane expects. Faded elements list their
pre-blended hex, since a PowerPoint line colour cannot carry alpha.

```bash
FIG_WIDTH_MM=190 python3 fig1_ppt.py                       # figure + PNG/PDF
EXPORT_SPEC=fig1_ppt_spec.json FIG_WIDTH_MM=190 python3 fig1_ppt.py
node build_pptx.js fig1_ppt_spec.json fig1_ppt_reproducible.pptx
python3 make_build_sheet.py fig1_ppt_spec.json             # BUILD_SHEET.md
```

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
pip install matplotlib                            # needs Liberation Sans (Arial-metric)
for w in 190 177.8; do                            # PALETTE=full|duo|mono
  FIG_WIDTH_MM=$w python3 fig1_lean.py            # submit this one
  FIG_WIDTH_MM=$w python3 fig1.py                 # annotated variant
done
```

Both scripts share `style.py`, which holds the palette, type scale and the hyperedge-blob geometry (rounded
convex offset = Minkowski sum of the members' convex hull with a disk). Layout
coordinates in `fig1.py` are millimetres on a fixed 190 × 158 canvas; `FIG_WIDTH_MM`
rescales the geometry while leaving type at its point size.
