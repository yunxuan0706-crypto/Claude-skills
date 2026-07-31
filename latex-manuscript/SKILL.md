---
name: latex-manuscript
description: Writes, formats, and troubleshoots LaTeX for academic manuscripts, theses, and grant documents — document structure, journal/conference templates, BibTeX/BibLaTeX citation management, equations, tables, and figures. Use this whenever the user is working with a .tex file, mentions LaTeX, needs to fit a manuscript into a specific journal's LaTeX template (e.g. Nature, IEEE, ACM, Springer, arXiv), asks about BibTeX/.bib files or citation commands, has a LaTeX compile error, or needs to convert a Word/plain-text draft into LaTeX. Also trigger for general academic document typesetting questions (equations, multi-column layout, cross-references) even if the user doesn't say "LaTeX" explicitly but the context (arXiv, journal submission, thesis) implies it.
---

# LaTeX Manuscript Preparation

## Why this matters

LaTeX is the default for most STEM journals and conferences because it
handles equations, cross-references, and citations reliably at scale — but
it fails loudly and unhelpfully, and the failure is often several lines
away from the actual mistake. The goal here is to get the document
compiling cleanly and matching the target venue's required format without
the user having to debug cryptic errors themselves.

## Workflow

1. **Identify the target venue and its template first.** Journals and
   conferences almost always provide a class file or template (e.g.
   `\documentclass{nature}`, an IEEEtran class, an ACM `acmart` class, a
   Springer `svjour3` class, arXiv's plain defaults). Ask the user which
   venue if it's not stated — using the generic `article` class for a
   submission that requires a specific class file will produce a document
   that gets desk-rejected or needs to be entirely reformatted later.

2. **Use the venue's required structure, not a generic one.** Most
   templates dictate section order, abstract placement, and where
   author/affiliation blocks go. Follow the template's own example file
   rather than inventing a structure, since journals often parse specific
   commands (e.g. author contribution blocks, competing-interest
   statements) that only work with their exact macros.

3. **Keep source organized for long documents.** For anything beyond a
   short paper, split into files with `\input{}`/`\include{}` (e.g.
   `intro.tex`, `methods.tex`, `results.tex`) rather than one giant `.tex`
   file — this makes diffs and collaboration much easier, and most
   journals don't require a single file until final submission.

## Citations (BibTeX / BibLaTeX)

- Maintain one `.bib` file with clean, consistent entries — every entry
  needs the fields the citation style actually uses (don't leave `year` or
  `journal` blank and expect the style to fill it in).
- Use consistent citation keys (`author2023keyword` or similar) so they're
  guessable and grep-able later.
- Match the citation package to what the venue expects: `natbib` for
  many journals (`\citep`, `\citet`), `biblatex` with a specified style
  for others, or a journal-specific package. Don't mix `natbib` and
  `biblatex` in the same document — pick the one the template uses.
- When converting a citation from a reference manager (Zotero, EndNote,
  Mendeley) export, check the exported `.bib` entry for garbage fields or
  malformed characters before compiling — exports are a common source of
  build failures.

## Common failure points and how to debug them

- **"Undefined control sequence"** — usually a missing `\usepackage{}` for
  the command being used, or a typo in the command name. Check the
  preamble first.
- **Citation shows as `[?]` or `??`** — the `.bib` file wasn't found, the
  citation key doesn't match exactly (case-sensitive), or `bibtex`/
  `biber` wasn't re-run after adding a new citation — LaTeX citation
  builds typically need multiple passes (`pdflatex` → `bibtex`/`biber` →
  `pdflatex` → `pdflatex`), not just one compile.
- **Figure/table not appearing where expected** — LaTeX's float placement
  is a request, not a command; `[htbp]` or `\FloatBarrier` (from the
  `placeins` package) or explicit float management helps but won't force
  exact placement — don't fight this by over-constraining unless the venue
  requires exact placement.
- **Overfull/underfull hbox warnings** — usually cosmetic (justification
  slightly off), but check if a long unbreakable string (URL, inline code,
  long equation) is causing it and needs `\url{}`, `\texttt{}` with
  breaking allowed, or manual line-breaking.
- **Package version mismatches after copying a template** — if a template
  was written for an older/newer LaTeX distribution, missing packages
  will fail cleanly with a "file not found" error naming the exact
  package to install.

## Equations and notation

- Use semantic commands over ad hoc formatting: define custom commands
  (`\newcommand{\vect}[1]{\mathbf{#1}}`) for notation used repeatedly, so
  a later change (e.g. switching vector notation) is a one-line edit
  instead of a find-and-replace across the document.
- Number only the equations that are referenced elsewhere in the text;
  suppress numbering (`\nonumber` or starred environments) for one-off
  display equations to reduce visual clutter.
- Keep consistent notation for the same quantity throughout — a variable
  that means one thing in Methods and something else in Results is a
  common and confusing error in math-heavy papers.

## Before calling it done

- Compile clean from scratch (delete `.aux`/`.bbl`/build artifacts and
  rebuild) to make sure the document isn't only working because of stale
  cached auxiliary files.
- Check the reference list renders in the exact style the venue requires
  (numbered vs. author-year, specific punctuation) — this is often
  checked manually by editorial staff and is an easy desk-reject trigger
  if wrong.
- Verify all cross-references (`\ref`, `\eqref`, `\cite`) resolve — no
  stray `??` anywhere in the compiled PDF.
