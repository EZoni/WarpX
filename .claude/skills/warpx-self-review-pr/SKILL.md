---
name: warpx-self-review-pr
description: Review the changes on the current branch relative to `development`, as a WarpX maintainer would review a pull request.
disable-model-invocation: true
---

# Self-Review a Pull Request

Review the changes on the current branch relative to the `development` branch, as if you were a WarpX maintainer reviewing a pull request.
This is a first pass that helps the author catch issues before requesting a review from other WarpX developers.
It does not replace the author's own critical review.

For best results, this skill should be run in a *fresh* session (not the one that wrote the code): an assistant asked to review its own work in the same session tends to defend earlier choices rather than scrutinize them.
The author should commit their work and make sure their branch is up to date with `development`, so that the diff the assistant reviews matches what reviewers will see.

## Review procedure

Review the changes on my current branch relative to the `development`
branch, as if you were a WarpX maintainer reviewing my pull request.
Do not make any changes yet; report your findings first.

Start by reading AGENTS.md for the project conventions, then run
`git fetch origin development` followed by
`git diff origin/development...HEAD` to see my changes. Fetching first
ensures the diff is taken against the latest upstream `development`,
not a stale local copy.

Check the following and report concrete issues with file and line references:

1. Correctness: logic errors, off-by-one/index mistakes, uninitialized
   values, incorrect physics or units, wrong sign conventions.
2. Algorithmic scaling: flag book-keeping or data-structure logic with
   worse asymptotic complexity than necessary, e.g. an O(N^2) loop over
   lists where an O(N) or O(N log N) approach exists.
3. Dimensionality: does the code handle all relevant builds
   (1D, 2D/XZ, 3D, RZ) correctly, including the compile-time macros?
4. GPU/CPU portability: any particle-to-grid deposition, scatter-add,
   histogram, or shared-counter loop must use `amrex::For` (not
   `amrex::ParallelFor`). Flag atomics that do not actually make a
   `ParallelFor` safe. See Docs/source/developers/portability.rst.
5. Backward compatibility: if a user-facing input parameter was removed
   or renamed, is there a guard in the relevant BackwardCompatibility()?
6. Testing: is there a test covering the new feature or bug fix? Is it
   fast enough for a 2-core CI runner and written portably?
7. Style: does the diff follow the C++/Python style in AGENTS.md, and
   does it avoid reformatting unrelated code?
8. Auto-generated files: flag any manual edits to `.pyi` stubs,
   `dependencies.json`, or `Regression/Checksum/benchmarks_json/*.json`.
9. Documentation: are new user-facing parameters or features documented?
10. Scope: is anything unrelated to the stated purpose of the PR included?

For each finding, state the severity (blocking / should-fix / nit) and
suggest a concrete fix. End with a short summary of whether this PR looks
ready for human review.

## Presenting the findings

Treat the output as a checklist of things the author should verify themselves, not as a verdict: it may raise false positives or miss real problems.
Confirm each finding against the code before acting on it.
