---
name: research
description: Conducts thorough, multi-source research and produces a well-organized, credible synthesis rather than a shallow single-source answer. Use this whenever the user asks a factual question that benefits from checking multiple sources, wants a research report, market/competitive/literature overview, background briefing, fact-check, or asks to "look into," "research," "find out about," or "dig into" a topic — especially when accuracy, currency, or source credibility matters (numbers, dates, claims about people/companies/events, technical comparisons). Also trigger when the user pushes back on an answer with "are you sure?" or "can you verify that?"
---

# Research

## Why this matters

A single search result or a single confident-sounding source is not
research — it's a guess with a citation. The value of doing this properly
is catching the cases where the first plausible-looking answer is wrong,
outdated, or one-sided, before that error reaches the user.

## Process

1. **Clarify the actual question before searching.** A vague topic
   ("research X") usually hides a more specific need — a decision the user
   is trying to make, a claim they want checked, a comparison they need.
   If it's ambiguous and the stakes seem meaningful, ask; otherwise infer
   the most useful scope and say what you assumed.

2. **Gather from multiple, independent sources.** One source is an
   anecdote. Aim for enough independent sources that the picture isn't
   resting on a single author's framing — independent means different
   authors/organizations, not five articles that all cite the same
   original press release.

3. **Evaluate each source before trusting it:**
   - **Who is saying this, and do they have a stake in the answer?**
     (vendor claiming their own product is best, advocacy group, primary
     research vs. aggregator vs. opinion piece)
   - **How recent is it, and does recency matter here?** Fast-moving topics
     (pricing, product specs, personnel, live events) need current sources;
     stable facts (historical events, established science) don't.
   - **Is this a primary source (the study, the filing, the transcript) or
     someone's summary of one?** Prefer primary when the specific numbers
     or quotes matter; a summary is fine for general context.

4. **Cross-verify anything load-bearing.** If a specific number, date,
   quote, or causal claim will drive the user's decision or appear as a
   headline fact, confirm it against at least one more source before
   presenting it as settled. If sources disagree, that disagreement is
   itself a finding — report it rather than silently picking one.

5. **Separate fact from interpretation from speculation** as you take
   notes, and keep that separation visible in the final synthesis. A
   claim like "the market is projected to grow" is an estimate someone
   made with assumptions attached, not a fact of the same kind as "the
   company was founded in 2015."

6. **Synthesize, don't just list.** The point of gathering from several
   sources is to produce a single coherent picture — what's well
   established, what's disputed, what's uncertain — not to dump source A's
   view, then source B's view, unintegrated.

## Report structure

For anything beyond a quick factual answer, organize the output like this
(adapt headers to fit a short chat reply vs. a longer written report):

- **Direct answer / key takeaway first.** Don't make the reader dig for the
  bottom line.
- **Findings**, organized by sub-topic or by question, with the load-bearing
  claims attributed to their sources (inline citation or link, not a bare
  list of URLs at the end with no mapping to claims).
- **Caveats and disagreements** — explicitly call out where sources
  conflict, where data is stale, or where confidence is genuinely low.
  Smoothing over disagreement to produce a cleaner-sounding answer is a
  failure mode, not a feature.
- **Sources**, listed so the user can follow up or verify independently.

## Failure modes to actively avoid

- Treating the first search result as the answer.
- Citing a source for a claim it doesn't actually make (check that the
  source really says what you're about to attribute to it).
- Presenting an estimate, projection, or someone's opinion with the same
  confidence as a verified fact.
- Silently dropping information that complicates a clean narrative — if
  something contradicts the emerging story, that's worth surfacing, not
  omitting.
- Stopping after confirming the answer once when the claim is significant
  enough to warrant a second, independent check.
