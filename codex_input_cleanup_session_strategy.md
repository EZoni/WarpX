# Codex Session Strategy for the WarpX Input-File Cleanup

## Purpose

This playbook explains how to use Codex to execute
`input_file_style_cleanup_plan.md` while keeping all work local. The user will review and push the
changes to GitHub manually.

The strategy uses fresh, bounded sessions and a repository progress document for handoff. This keeps
each change set reviewable and reduces context drift during a repository-wide migration. It follows
the official OpenAI recommendation to give coding agents explicit context, constraints, approval
boundaries, success criteria, and stopping conditions:

<https://developers.openai.com/api/docs/guides/latest-model>

## Recommended number of sessions

Plan for approximately 12-16 fresh Codex sessions:

| Session | Purpose | Expected count |
| --- | --- | ---: |
| 1 | Inventory and style convention | 1 |
| 2 | Conservative style checker | 1 |
| 3 | Representative pilot migration | 1 |
| 4 onward | One reviewable migration batch per session | 8-12 |
| Final | Completion and semantic-safety audit | 1 |

The exact total should be recalculated after Session 1. A reasonable migration batch is approximately
5-10 directories or 20-40 files. Base files and their dependent override files must stay in the same
batch, even when that changes the batch size.

Do not use one session for the entire project. Start a fresh session at each phase boundary and for
each migration batch. Within a session, use one main prompt and add follow-up prompts only to resolve
specific failures, findings, or decisions.

## Persistent handoff artifacts

The following repository files carry context between sessions:

- `AGENTS.md`: repository-specific working and validation instructions;
- `input_file_style_cleanup_plan.md`: scope, style goals, phases, and acceptance criteria;
- the developer-facing input style document created in Session 1; and
- `input_file_style_cleanup_progress.md`: inventory, decisions, completed batches, validation, skipped
  coverage, exceptions, and the next recommended batch.

Every session must read these applicable files before editing and update the progress document before
finishing. The progress document is the authoritative handoff between fresh sessions; do not rely on a
previous chat being available.

## Git and GitHub boundary

Append the following block to every session prompt:

```text
Work only in the local checkout.

Do not use GitHub or any other remote service. Do not run `git push`,
`git pull`, `git fetch`, `gh`, or any command/API that creates or changes
pull requests, issues, comments, releases, branches on a remote, or remote
repository state. Do not install or request a GitHub plugin. Do not change
git remotes.

Local inspection, file edits, builds, tests, and local git commits are
authorized. Never modify, discard, stage, amend, or commit unrelated
worktree changes. I will review and push changes manually.
```

Local commits are recommended because they make each completed batch independently reviewable. Codex
must not amend user commits or include unrelated worktree changes.

If the user also wants to create all local commits manually, replace the final paragraph of the block
with:

```text
Local inspection, file edits, builds, and tests are authorized. Do not run
`git add`, `git commit`, `git rebase`, or `git reset`. Leave all intended
changes unstaged for my review. Never modify or discard unrelated worktree
changes.
```

## Session 1 prompt: inventory and convention

Use this prompt once in a fresh session:

```text
Read AGENTS.md and input_file_style_cleanup_plan.md completely.

Execute only Phases 1 and 2 of the plan: create the inventory and dependency
map, then add the proposed developer-facing input-file style documentation.
Do not implement the checker or migrate input files yet.

Create input_file_style_cleanup_progress.md as the durable handoff document
for future Codex sessions. Record:
- inventory totals and classifications;
- FILE inheritance groups;
- exceptional syntax and semantic-ordering risks;
- mapping to relevant CTest tests;
- files changed;
- validation performed;
- open questions and recommended migration batches.

Follow the plan's safety invariants. Inspect existing documentation and
representative files before deciding details. Run appropriate documentation
and formatting checks. At the end, review the full diff, report results and
remaining risks, and make one local commit if all relevant checks pass.

[APPEND THE GIT AND GITHUB BOUNDARY]
```

### User review gate after Session 1

Before continuing, manually review:

- the canonical section order;
- native versus PICMI/Python rules;
- the comment policy;
- inheritance and exceptional-syntax findings;
- proposed migration batches; and
- whether the documentation is specific enough without overconstraining specialized tests.

Correct the convention now rather than after it has been applied broadly.

## Session 2 prompt: conservative checker

Use this prompt once in a fresh session:

```text
Read AGENTS.md, input_file_style_cleanup_plan.md, and
input_file_style_cleanup_progress.md completely. Inspect the current local
git history and worktree before editing.

Execute only Phase 3 of the plan. Implement and test a conservative
input-style checker. It must enforce only semantic-neutral, objectively
detectable rules. It must not reorder parameters, rewrite expressions,
require comment density, or treat Python syntax as native input syntax.

Add focused tests or fixtures for normal assignments, quoted expressions,
inline comments, continuations, parser expressions, FILE inheritance, and
Python inputs. Do not migrate the Examples tree merely to make the new
checker pass; if necessary, introduce the checker in reporting/non-blocking
mode until migration is complete.

Update input_file_style_cleanup_progress.md with design decisions, commands,
results, and unresolved exceptions. Review the diff and run relevant checks.
Make one local commit if the work is complete and validated.

[APPEND THE GIT AND GITHUB BOUNDARY]
```

### User review gate after Session 2

Confirm that the checker:

- does not reorder or rewrite files by default;
- distinguishes Python from native inputs;
- handles quoted strings, comments, continuations, and parser expressions;
- produces actionable file-and-line diagnostics; and
- is non-blocking until the existing tree has been migrated, if it does not yet pass globally.

## Session 3 prompt: representative pilot

Use this prompt once in a fresh session:

```text
Read AGENTS.md, input_file_style_cleanup_plan.md,
input_file_style_cleanup_progress.md, and the new style documentation
completely. Inspect the current local history and worktree.

Execute only Phase 4. Select the representative pilot described by the plan,
or use the pilot already selected in the progress document. Keep base files
and every affected override together.

Before editing, capture the effective assignment order, repeated assignments,
FILE inheritance, CMake invocations, and available associated tests. Apply
the documented style without changing parameter names, values, expressions,
quoting with semantic significance, ordering dependencies, or Python
execution order.

Run the checker, pre-commit checks on changed files, and all associated
runnable CTests. Investigate every non-checksum failure. Do not update
checksums or analysis tolerances.

Update the progress document with every changed file, test result, skipped
configuration, exception, and any convention/checker refinement. Review the
diff for semantic changes. Make one local commit if validated. Do not begin
the broader migration.

[APPEND THE GIT AND GITHUB BOUNDARY]
```

### User review gate after Session 3

Review the pilot line by line. In particular, compare values and expressions with the parent version,
verify override precedence, and decide whether the resulting comment density and headings are useful.
Do not begin batch migration until the pilot establishes a trustworthy pattern.

## Migration-batch prompt

Reuse this prompt in a fresh session for every Phase 5 batch:

```text
Read AGENTS.md, input_file_style_cleanup_plan.md, the input-file style
documentation, and input_file_style_cleanup_progress.md completely. Inspect
the current local history and worktree.

Execute exactly one next uncompleted Phase 5 migration batch. Choose the next
batch recorded in the progress document. If none is specified, select a
coherent group of approximately 5-10 directories or 20-40 files. A FILE base
file and all affected dependent overrides must remain in the same batch, even
if that requires adjusting the size.

Do not revisit completed batches unless validation reveals a concrete
problem. Preserve every semantic-safety invariant from the plan. Add comments
only when they explain test intent, physics, units, a non-default choice, CI
scaling, or a production alternative.

For this batch:
1. inspect its CMake declarations, base/override relationships, README files,
   and analysis scripts;
2. record effective assignments and ordering-sensitive constructs;
3. perform the style cleanup;
4. run the style checker and relevant pre-commit checks;
5. run every associated CTest supported by the local build;
6. inspect both ordinary and whitespace-aware diffs;
7. verify that no checksum, tolerance, generated, or unrelated file changed.

Update input_file_style_cleanup_progress.md with the exact scope, validation,
skipped tests, exceptions, and recommended next batch. Make one local commit
if all relevant checks pass. Stop after this single batch.

[APPEND THE GIT AND GITHUB BOUNDARY]
```

### User review gate after each migration batch

Before starting the next session:

1. Read the Codex completion report and progress entry.
2. Inspect the local commit or unstaged diff.
3. Confirm that only the declared batch and progress file changed.
4. Spot-check parameter values, expressions, inheritance, and comments.
5. Confirm that test failures and skipped configurations are explained.
6. Fix or revert an unsatisfactory batch before allowing later batches to build on it.

The instruction to stop after one batch is intentional. It establishes a stable review point and
prevents a session from silently expanding across the entire tree.

## Final-session prompt: completion audit

Use this prompt once after the progress document says every migration batch is complete:

```text
Read AGENTS.md, input_file_style_cleanup_plan.md, the input-file style
documentation, and input_file_style_cleanup_progress.md completely. Inspect
the full series of local commits and the current worktree.

Execute only Phase 6 and the plan's final acceptance audit. Do not perform
broad cleanup or redesign during this session.

Recompute the complete inputs inventory and verify that every file is either
compliant or has a documented narrow exception. Run the style checker, the
existing input naming/registration checker, applicable pre-commit checks,
and the broadest practical CTest matrix supported by the local build.

Audit the complete change series for:
- changed parameter names, values, expressions, or list ordering;
- changed override or duplicate-assignment precedence;
- altered Python execution order;
- checksum or analysis-tolerance modifications;
- missing inputs or CMake registration changes;
- unrelated modifications;
- incomplete or inaccurate documentation.

Fix only clearly in-scope defects found by the audit, validate those fixes,
and update input_file_style_cleanup_progress.md with final results and
unavailable test coverage. Provide a concise completion report. Make a final
local audit commit only if files needed correction. Do not push anything.

[APPEND THE GIT AND GITHUB BOUNDARY]
```

## Handling incomplete or failed sessions

If a session stops before completing its assigned phase or batch, do not start the next planned
session as though the work succeeded. Start a recovery session with this prompt:

```text
Read AGENTS.md, input_file_style_cleanup_plan.md, the input-file style
documentation if it exists, and input_file_style_cleanup_progress.md. Inspect
the current worktree and recent local commits.

The previous Codex session did not complete its assigned work. Diagnose the
current state without discarding changes. Identify completed, incomplete,
unvalidated, and unrelated modifications. Continue only the interrupted
phase or batch, run its required validation, and update the progress document.
Do not begin the next phase or migration batch.

[APPEND THE GIT AND GITHUB BOUNDARY]
```

If the incomplete work is unsafe or cannot be understood, stop and report the exact files and risks
instead of resetting or overwriting it.

## Final manual workflow

After the completion audit:

1. Inspect `git status` and the complete local commit series.
2. Review the diff against the branch point from `development`.
3. Run any additional build matrix available outside the Codex environment.
4. Squash or reorganize local commits only if desired; preserve reviewable logical boundaries.
5. Push the branch and create the pull request manually.
6. Include skipped platforms/features and the semantic-preservation strategy in the pull request
   description.

Codex is not authorized to perform steps 4-6 unless the user later gives separate, explicit
authorization. Under the default workflow in this playbook, all remote GitHub activity remains manual.
