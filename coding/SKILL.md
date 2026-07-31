---
name: coding
description: Guides writing, reviewing, and revising code toward simplicity, correctness, and maintainability rather than cleverness or premature abstraction. Use this whenever the user asks to write, refactor, review, or debug code in any language, wants feedback on code quality, asks "is this good code" or "how would you improve this," or is deciding between implementation approaches. Also trigger before finishing a coding task to run a pre-submit self-review pass.
---

# Coding

## Why this matters

Code is read far more often than it's written, and it's modified by people
(including a future version of the author) who don't have the current
context in their head. The practices below exist to keep that future
reader's job easy: minimize what they have to hold in their head to trust
that the code does what it says.

## Core principles

- **Solve the actual problem, not the anticipated future one.** Don't add
  configurability, abstraction layers, or generality for requirements that
  don't exist yet. Three similar lines of code are better than a premature
  helper function that has to guess at what the next caller will need.
  YAGNI is a default, not an exception.
- **Prefer the boring solution.** A clever one-liner that requires the
  reader to pause and decode it is worse than three plain lines, unless the
  cleverness is idiomatic and every reader in this codebase will recognize
  it instantly.
- **Name things for what they are, not how they're used today.** Good
  names carry meaning without needing a comment. If a name needs a comment
  to explain what it holds, the name is wrong.
- **Comments explain why, never what.** Well-named code already says what
  it does. A comment earns its place only when it captures something the
  code can't: a non-obvious constraint, a workaround for a specific bug, an
  invariant that would surprise the next editor. If deleting the comment
  wouldn't make the code more confusing, delete it.
- **Handle errors at boundaries, not everywhere.** Validate and defend at
  the edges of a system — user input, network calls, external APIs. Trust
  your own internal functions and language/framework guarantees rather
  than wrapping everything in defensive checks for states that can't
  actually occur.
- **Don't build compatibility shims for code you can just change.** Inside
  a codebase you control, rename the thing and update the callers rather
  than keeping an old alias around "just in case." Shims accumulate and
  nobody remembers why they exist.
- **Match the surrounding style.** Consistency within a file or module
  matters more than any individual stylistic preference — a codebase that
  reads uniformly is easier to trust than one where every file reflects a
  different author's taste.

## Correctness and testing habits

- **Understand the actual requirement before writing code**, including the
  edge cases that matter (empty input, concurrent access, the boundary
  values) — not every conceivable edge case, just the ones that can really
  happen given how the code will be called.
- **Write tests that would catch a real regression**, not tests that just
  exercise the code for coverage's sake. A test that can't fail if the
  logic breaks isn't testing anything.
- **Run the code (or the tests) before declaring the task done.** Type
  checking and a passing test suite verify the code is well-formed and
  matches the spec you wrote — they don't verify the feature actually works
  for a real user. Exercise the actual path when you can.
- **When fixing a bug, understand the root cause before patching the
  symptom.** A fix that makes the specific failing case pass without
  understanding why it failed is likely to leave the underlying issue
  intact for the next input.

## Pre-submit review checklist

Before treating a change as finished, scan for:

1. **Scope creep** — did the change stay focused on what was asked, or did
   it drag in unrelated refactoring, renames, or "while I'm here" cleanup?
2. **Dead code** — anything left behind that's now unused (old branches,
   unused imports/variables, commented-out code)?
3. **Duplicated logic** — is the same non-trivial logic now written in two
   places that could drift apart? (Note: a *little* duplication is fine —
   don't over-correct into premature abstraction; this is about real,
   nontrivial repeated logic.)
4. **Unhandled edge cases that can actually occur** — not every
   hypothetical, but the ones reachable from real callers/inputs.
5. **Naming** — would someone unfamiliar with this change understand each
   name without needing the surrounding context explained to them?
6. **Comments** — does every comment explain a *why* that isn't obvious
   from the code, and nothing else?
7. **Security basics at any boundary touched** — injection (SQL, shell,
   command), unsanitized input reaching a sink, secrets in code/logs,
   missing auth checks on a new endpoint.

If a change fails one of these, fix it before calling the task done rather
than leaving a known issue for later — "later" rarely comes back to small
things like this.
