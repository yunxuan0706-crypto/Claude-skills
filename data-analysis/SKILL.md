---
name: data-analysis
description: Helps choose appropriate statistical methods, run and interpret analyses correctly, and report results the way a research audience expects (effect sizes, confidence intervals, appropriate tests, honest treatment of uncertainty). Use this whenever the user has research/experimental data and needs to analyze it, asks which statistical test to use, wants help interpreting a p-value or model output, is checking whether an analysis is valid, needs to report results for a paper's Results section, or is deciding on a study's analysis plan. Also trigger when the user's proposed analysis shows signs of a common pitfall (multiple comparisons without correction, small sample overinterpretation, correlation framed as causation) even if they didn't ask about it directly.
---

# Data Analysis for Research

## Why this matters

A statistically invalid analysis doesn't just risk a wrong number — it
risks the paper's core claim, and it's one of the most common reasons
reviewers reject or demand major revisions. Getting the method right,
and reporting it honestly, protects both the science and the author's
credibility.

## Choosing a method

Match the method to the question and the data, not to what's familiar or
what a colleague used on a different problem:

- **What's the outcome variable's type?** Continuous, count, binary,
  ordinal, and time-to-event outcomes call for different models
  (linear regression vs. Poisson/negative binomial vs. logistic vs.
  ordinal regression vs. survival analysis). Using linear regression on a
  binary outcome, for instance, is a common and avoidable error.
- **What's the actual comparison being made?** Two groups vs. more than
  two, paired vs. independent samples, one predictor vs. many — this
  determines t-test vs. ANOVA vs. regression vs. mixed-effects models, and
  getting it wrong (e.g. running several pairwise t-tests instead of one
  ANOVA) inflates false-positive risk.
- **Are observations actually independent?** Repeated measures on the same
  subject, clustered data (patients within hospitals, students within
  classrooms), or time-series data violate the independence assumption
  most basic tests rely on — these need mixed-effects/hierarchical models
  or clustered standard errors instead.
- **Check assumptions before trusting the test**, and know the fallback:
  normality (for small samples), homogeneity of variance, linearity — when
  violated, either transform the data, use a robust/nonparametric
  alternative, or use a model that doesn't assume it, rather than ignoring
  the violation.

## Avoiding the common pitfalls

- **Multiple comparisons.** Running many tests (many genes, many outcome
  measures, many subgroups) inflates the false-positive rate unless
  corrected for (Bonferroni, Benjamini-Hochberg/FDR, or a pre-registered
  primary-outcome approach). If the user is running several tests and
  hasn't mentioned correction, raise it.
- **p-hacking / analysis flexibility.** Trying multiple analysis
  approaches and reporting the one that "worked" invalidates the p-value's
  meaning. If an analysis plan was pre-specified, stick to it and report
  deviations explicitly; if it wasn't, be transparent that the analysis
  was exploratory rather than confirmatory.
- **Small samples and overinterpretation.** A striking effect size from
  n=6 is much less trustworthy than the same effect size from n=200 — pair
  point estimates with confidence intervals so the uncertainty is visible,
  and be cautious about strong causal language from underpowered studies.
- **Correlation vs. causation.** Observational associations support
  causal claims only with a design or method that addresses confounding
  (randomization, instrumental variables, well-justified adjustment,
  natural experiments). Flag when the proposed language claims more than
  the design supports.
- **Post-hoc storytelling.** Distinguish confirmatory analyses (testing a
  pre-specified hypothesis) from exploratory ones (finding patterns in the
  data) in both the analysis and the write-up — exploratory findings need
  to be reported as hypothesis-generating, not confirmed.

## Reporting results

Research audiences expect more than a bare p-value:

- Report **effect sizes** (mean difference, odds ratio, Cohen's d, R²,
  etc.) alongside significance — a significant but tiny effect and a
  significant, large effect tell very different stories, and the p-value
  alone doesn't distinguish them.
- Report **confidence/credible intervals**, not just point estimates —
  this communicates precision, which matters as much as the estimate
  itself.
- State the **exact test used and why**, including software/package and
  version if precision matters for reproducibility.
- Report **sample sizes** at each stage, especially if data was excluded
  — undisclosed exclusions are a common source of reviewer skepticism.
- Use **exact p-values** (p = 0.032) rather than only threshold statements
  (p < 0.05) except for very small values, where "p < 0.001" is standard.
- When results are null/non-significant, report them as such rather than
  omitting them or reframing as "trending toward significance" — a
  non-effect is still a finding.

## Working with code for analysis

When writing analysis code (R, Python/pandas, etc.):

- Keep the analysis script reproducible end-to-end from raw data to final
  numbers/figures — a reviewer or collaborator should be able to rerun it
  and get the same result.
- Set and record a random seed for anything stochastic (bootstrapping,
  cross-validation splits, simulations).
- Separate data cleaning from analysis from visualization into distinct,
  clearly ordered steps rather than one tangled script, so each stage can
  be checked independently.
- Sanity-check intermediate outputs (row counts before/after filtering,
  summary statistics) rather than trusting the final number blindly — the
  most common source of a wrong result is a silent data-wrangling bug
  upstream of the statistical test itself.
