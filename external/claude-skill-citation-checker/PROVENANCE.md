# Provenance

This directory is a **vendored, locally-modified copy** of a third-party
Claude skill, not an unmodified mirror.

- Upstream: https://github.com/PHY041/claude-skill-citation-checker
- Vendored from commit: 9911ca60750079f8ffc8382e7cf0fe565bcc2950 (2026-03-22)
- Upstream declared no LICENSE file at the time of vendoring.

## Why it is vendored instead of a git submodule

We needed to modify `scripts/citation_checker.py` to add rate-limit
resilience (429 exponential backoff with a capped Retry-After), an
OpenAlex/CrossRef polite-pool `mailto`, and an optional Semantic Scholar
API key read from the environment. Those changes must be committed to this
repository, and we cannot push them to the upstream repo, so the skill was
converted from a submodule to a vendored copy.

## Local modifications

See the "LOCAL MODIFICATIONS" banner near the top of
`scripts/citation_checker.py`. Behavioural changes only; the matching /
confidence logic is unchanged from upstream.

Env vars are read through `_clean_env`, which discards any value still
containing angle-bracket placeholder text (e.g. `<your-key>`) so a
copy-pasted doc placeholder is never sent as a real credential.

A `rate_limited` verification status was added (upstream has only
verified / suspicious / not_found / error). When a citation ends up in no
source *only because a source returned HTTP 429*, it is reported as
"could not verify" rather than "not_found" — so a transient rate limit on
a shared IP never mislabels a real reference as a hallucination. The
matching / confidence thresholds are still unchanged; this only reclassifies
the empty-result case based on whether sources were actually reachable.
