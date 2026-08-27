# Plan: Standardize WarpX Example Input Files

## Objective

Standardize the organization, whitespace, headings, and comments of test input files under
`Examples/` without changing simulation behavior. Document the convention, migrate the files in
reviewable batches, and add narrowly scoped automated checks that prevent the style from drifting.

This is a style and maintainability project. Do not change physics, test coverage, parameter values,
test tolerances, checksums, or runtime behavior as part of this work.

## Repository context and constraints

- WarpX input files use AMReX `ParmParse` syntax. Parameter ordering can be significant when a file
  contains repeated assignments, parser constants, or overrides inherited through `FILE = ...`.
- Shared `inputs_base*` files and their dependent files must be reviewed as a unit.
- Python inputs include both PICMI simulations and preparation/wrapper scripts. They need a convention
  appropriate for Python rather than a literal copy of the native-input layout.
- Tests span dimensionalities and optional backends/features. Only run tests supported by the active
  build, and record tests that could not be run locally.
- Ignore checksum-only failures when running tests, as required by the repository guidance, but
  investigate numerical, parsing, runtime, or analysis failures.
- Do not modify assertion tolerances merely to make a test pass.
- Do not modify generated `.pyi` files, `dependencies.json`, or regression checksum JSON files.
- Preserve unrelated changes in the worktree.
- Consult `Docs/source/glossary.rst` when writing physics-facing comments.

An inventory taken on 2026-08-27 found approximately:

- 405 files whose basename starts with `inputs` in 100 directories;
- 332 native WarpX/AMReX input files;
- 73 Python/PICMI or preparation scripts;
- 30 shared `inputs_base*` files; and
- 123 native files containing `FILE = ...` inheritance.

Recompute this inventory before implementation; these numbers are planning context, not fixed
acceptance criteria.

## Desired style

### General principles

1. Optimize files for readers trying to understand what the test exercises.
2. Keep related parameters together, but never reorder them until semantic dependencies and
   overrides have been checked.
3. Use comments to explain intent, units, physics, non-obvious choices, reduced test settings, and
   production alternatives. Do not comment every assignment or merely restate a parameter name.
4. Preserve useful existing comments. Rewrite only when needed for accuracy, clarity, or consistency.
5. Keep diffs limited to the input-style work. Do not opportunistically modernize test logic.

### Native input files

Use the following section order when the sections are present and when reordering is demonstrably
safe:

1. File inheritance and parser constants
2. Run control
3. Geometry and mesh
4. Mesh refinement
5. Boundary conditions
6. Numerical algorithms and field solvers
7. Applied fields and lasers
8. Particles and fluids
9. Collisions, reactions, and other interactions
10. Diagnostics
11. Reduced diagnostics

The list is a semantic ordering guide, not a requirement to add empty sections. Specialized tests may
add clearly named sections where appropriate.

Formatting rules:

- Use a simple heading such as `# Geometry and mesh`.
- Do not use boxed or decorative ASCII headings made from repeated `#` characters.
- Use exactly one space on each side of the assignment operator for ordinary assignments.
- Do not align assignment operators by inserting variable amounts of whitespace.
- Separate top-level sections with one blank line.
- Use blank lines within a section only to separate meaningful subgroups, such as individual species.
- Prefer a full-line comment above an assignment when an explanation does not fit comfortably inline.
- Keep inline comments short and separate them from the value with two spaces.
- Avoid broad normalization of quoting, capitalization, numeric notation, expression spelling, or
  parameter values. These changes make semantic review harder and may alter parsing.
- Preserve line continuations and expression strings unless a targeted change has been verified.
- Do not impose a rigid line-length rule on long parser expressions. Improve readability only when the
  input grammar supports the chosen continuation syntax and the associated test validates it.

### Comment policy

Good comments answer at least one of these questions:

- What physical quantity or unit is represented?
- Why is this non-default setting needed by the test?
- What behavior or code path is the test intended to cover?
- Why is a deliberately coarse/short configuration used in CI?
- What production value or configuration would normally differ?
- Why must this parameter appear before or after another assignment?

Avoid comments such as `# number of cells` immediately above `amr.n_cell`, unless they add information
that is not already apparent. Do not add comments solely to reach a target comment count.

### Python and PICMI input files

- Continue to follow the repository's Ruff configuration for ordinary Python formatting.
- Organize imports, physical constants, distributions, grids/solvers, species, interactions, lasers,
  diagnostics, simulation construction, and execution in a consistent logical order when applicable.
- Preserve command-line parsing, preparation-script setup, callback registration, and execution order.
- Use brief comments for test intent and non-obvious WarpX/PICMI extensions.
- Use docstrings for reusable helpers where they add information; do not add boilerplate docstrings to
  trivial test-local functions.
- Do not force preparation scripts into the same section structure as PICMI simulation scripts.

## Non-goals

- Renaming input files, tests, species, diagnostics, or parser constants.
- Replacing native input files with PICMI, or vice versa.
- Changing defaults, numerical methods, resolutions, time steps, diagnostic contents, or output formats.
- Consolidating or splitting tests unless separately approved.
- Updating expected checksums.
- Automatically sorting parameters alphabetically.
- Reformatting all Python or documentation outside the input-file scope.

## Safety invariants

The following must remain unchanged unless a separate bug is discovered, documented, and explicitly
approved:

- The ordered sequence of effective parameter assignments after recursively processing `FILE = ...`.
- Repeated-assignment precedence.
- Parser constant definitions appearing before any use for which ordering matters.
- All parameter names, values, expression strings, list ordering, and quoting with potential semantic
  significance.
- CMake test names, dimensions, process counts, input arguments, runtime parameters, and analysis
  commands.
- Python import side effects, callbacks, command-line arguments, and simulation execution order.

When unsure whether a reordering is safe, retain the existing order and improve only headings,
whitespace, and comments.

## Implementation phases

### Phase 1: Inventory and dependency map

1. Enumerate all files under `Examples/` whose basename starts with `inputs`.
2. Classify each as native input, PICMI simulation, Python wrapper, preparation script, or shared base
   file.
3. Parse or conservatively scan native files for `FILE = ...`, repeated assignments, parser constants,
   continuations, and unusual syntax.
4. Build a map from every base file to its direct and transitive dependents.
5. Map each `inputs_test_*` file to its `add_warpx_test()` declaration and note dimensions, optional
   features, process counts, runtime overrides, and analysis scripts.
6. Record exceptional files that cannot safely follow the default section order.

Keep generated inventory artifacts temporary unless they are useful as maintainable inputs to the
eventual checker.

### Phase 2: Document the convention

1. Add a concise developer-facing style document, preferably under `Docs/source/developers/`, and add
   it to the appropriate documentation toctree.
2. Include one native example and one PICMI example showing the desired result.
3. Document the comment policy and the rule that semantic preservation takes priority over canonical
   ordering.
4. Document how base files and override files should be formatted together.
5. Get agreement on the convention before converting the full tree.

### Phase 3: Implement a conservative checker

Add or extend a repository check for rules that can be enforced without interpreting WarpX physics.
Prefer a dedicated checker if extending `.github/workflows/source/check_inputs.py` would mix naming and
formatting concerns excessively.

Initially enforce only objective rules such as:

- no decorative ASCII section banners;
- no trailing whitespace;
- normalized spacing around ordinary assignment operators where the grammar is unambiguous;
- no excessive consecutive blank lines;
- recognized simple section-heading syntax; and
- Python files passing the existing Ruff checks.

Do not make the checker:

- reorder assignments;
- require every canonical section;
- require comments for every parameter;
- enforce a comment-density threshold;
- rewrite quoted expressions; or
- apply a simplistic assignment regex to Python, comparison operators, or expression contents.

Provide actionable diagnostics containing the file, line, violated rule, and expected form. If an
automatic fixer is added, keep it opt-in and restrict it to provably semantic-neutral transformations.

### Phase 4: Pilot conversion

Select a small but representative pilot containing:

- one standalone native input;
- one base file with multiple override files;
- one PICMI input; and
- one specialized file with long expressions or unusual syntax.

For each pilot file:

1. Capture the original effective assignment order and CMake invocation.
2. Apply the proposed organization, headings, spacing, and comment policy.
3. Review the diff specifically for changed values, expression text, precedence, or execution order.
4. Run the input checker, relevant formatting checks, and all associated CTests available in the local
   build.
5. Refine the written convention and checker based on genuine edge cases found during the pilot.

Do not proceed to broad conversion until the pilot demonstrates that the convention is clear and the
validation process detects accidental semantic changes.

### Phase 5: Batched migration

Migrate by example directory or tightly related directory group. Keep each pull request or commit small
enough for line-by-line review; roughly 5-15 directories is a reasonable initial target, adjusted for
file size and dependency complexity.

For every batch:

1. Update base files and all affected override files in the same batch.
2. Preserve the safety invariants above.
3. Add only useful, technically accurate comments.
4. Run the style checker and repository pre-commit checks on changed files.
5. Run every associated test supported by the configured build.
6. Review the diff with whitespace ignored and again with whitespace visible.
7. Confirm that no checksums, analysis tolerances, or unrelated files changed.
8. Record tests skipped because their dimension, backend, or optional dependency is unavailable.

A suggested migration order is:

1. Small standalone native tests.
2. Small PICMI-only tests.
3. Families with simple base-file inheritance.
4. Large physics applications used in documentation.
5. Complex restart, multi-level, PSATD, embedded-boundary, QED, and collision families.
6. Preparation scripts, wrappers, and remaining exceptional cases.

### Phase 6: CI enforcement and completion audit

1. Run the final checker across all of `Examples/`.
2. Run `.github/workflows/source/check_inputs.py` to confirm naming and registration remain valid.
3. Run the broadest practical CTest matrix across supported dimensionalities and features.
4. Verify every inventoried file is either compliant or has a documented, narrowly scoped exception.
5. Enable the style checker in CI only after the migration is complete.
6. Document how contributors can run the checker locally and how to request an exception.

## Validation strategy

### Static validation

- Compare the set of input filenames before and after each batch.
- Compare native parameter names, values, and their effective order, accounting for `FILE` inheritance.
- Flag added, removed, or reordered duplicate assignments for manual review.
- Compare CMake test declarations and command-line overrides.
- Run Ruff/pre-commit on changed Python files.
- Run the existing input-name and registration checker.

The comparison tool must understand comments, quoted strings, continuations, and inheritance well enough
to avoid claiming semantic equivalence based only on sorted `key=value` text. If robust comparison is
not possible for a file, use a manual review plus runtime validation.

### Runtime validation

- Configure/build inside the repository root according to `AGENTS.md`.
- Use `ctest --test-dir build -R <test-regex> --output-on-failure` for associated tests.
- Do not create or search for build directories outside the repository/worktree root.
- If `cmake` or `ctest` is unavailable, activate the `warpx-cpu-mpich-dev` environment as described in
  `AGENTS.md`.
- Treat parser errors, changed analysis results, crashes, hangs, and non-checksum test failures as
  blockers for the batch.
- Treat checksum-only failures as platform-sensitive evidence to inspect, not as permission to update
  checksum files.

## Review checklist for each changed file

- [ ] The file's purpose is apparent from its headings and a small number of useful comments.
- [ ] Sections follow the preferred order where safe.
- [ ] No empty or irrelevant canonical section was added.
- [ ] No decorative banner remains.
- [ ] Assignment and blank-line formatting follow the convention.
- [ ] Comments add intent, units, rationale, or test/production context.
- [ ] Parameter names, values, quoting, expressions, and list order are unchanged.
- [ ] Override and duplicate-assignment precedence is unchanged.
- [ ] Base-file dependencies were reviewed together.
- [ ] Python execution order and side effects are unchanged.
- [ ] The relevant static checks pass.
- [ ] Associated runnable tests pass, apart from understood checksum-only differences.
- [ ] No generated checksum, tolerance, or unrelated file changed.

## Acceptance criteria

The project is complete when:

1. A documented native-input and Python/PICMI style convention exists and is linked from the developer
   documentation.
2. Every `inputs*` file under `Examples/` has been reviewed and is either compliant or covered by a
   documented exception.
3. Decorative headings and inconsistent assignment spacing have been removed where safe.
4. Sparse files contain enough explanation to communicate test intent, without boilerplate or a
   comment-density quota.
5. Shared base files and overrides have consistent organization and unchanged precedence.
6. A conservative checker produces actionable diagnostics and passes across the tree.
7. Existing input naming/registration checks pass.
8. Relevant tests pass across the available build matrix, with skipped configurations documented.
9. No regression checksum, analysis tolerance, or unrelated source file was changed.
10. The final diff can be justified entirely as documentation/style-preserving maintenance.

## Instructions for the implementing coding assistant

- Start with Phase 1; do not immediately rewrite all files.
- Before editing, inspect `AGENTS.md`, the relevant test `CMakeLists.txt`, any `README.rst`, analysis
  script, base file, and dependent override files.
- Make reasonable progress without asking about trivial formatting choices, but stop and report any
  ambiguity that could change simulation behavior.
- Use small patches and report the files, checks, tests, and skipped coverage after every batch.
- Never infer semantic equivalence merely because CTest reaches a checksum comparison.
- If an existing comment appears inaccurate, verify it against the parameter documentation or source
  before rewriting it.
- If the work exposes an unrelated test bug, record it separately and keep it out of the style patch
  unless explicitly authorized.
