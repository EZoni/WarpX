.. _developers-llm:

LLM-Assisted WarpX Development
==============================

Large Language Models (LLMs) can assist with WarpX development tasks such as navigating the codebase, understanding existing implementations, writing new features, debugging, adding tests, and preparing pull requests.
This guide documents how the WarpX repository is configured for LLM-based coding assistants and how to get the most out of them.

.. note::

   LLM-generated code should always be reviewed carefully.
   LLMs can hallucinate APIs, miss physics constraints, or produce code that compiles but is incorrect.
   The configurations described below help by providing accurate, up-to-date context, but developer judgment remains essential.

   The WarpX community thus urges you to perform a careful, manual review of all LLM-generated code and documentation before asking for a review of your pull-request.
   This is important, otherwise you risk to waste the valuable time of our most proficient developers that will need to review your LLM-generated code.
   (Be considerate that WarpX developers can prompt an LLM just as efficiently as you can. Your critical thinking skill to make sense of the LLM-generated code and make it sensible for review and maintainable for the long term is what is needed!)

.. note::

   This section is not understood as an endorsement of any of the listed (or unlisted) coding assistants or MCP services.
   Contributions to this section documenting further services, clients, skills, etc. are encouraged.

.. tip::

   When working with LLM coding assistants, keep in mind that *"most best practices are based on one constraint: [the] context window fills up fast, and performance degrades as it fills"* (`Claude Code Best Practices <https://code.claude.com/docs/en/best-practices>`__).
   Keep instructions concise, use plan modes, break complex tasks into focused sessions, and provide targeted context rather than overwhelming the assistant with information.


AGENTS.md / CLAUDE.md
---------------------

The repository includes an ``AGENTS.md`` file at its root (as well as a ``CLAUDE.md``, which directly points to ``AGENTS.md``.)
These files are automatically loaded by LLM coding assistants (Claude Code reads ``CLAUDE.md``; other tools such as OpenAI Codex CLI read ``AGENTS.md``) to provide project-specific instructions.

The file contains in a compressed form instructions for an LLM agent:
With this file present, an LLM assistant working inside the WarpX repository will automatically know how to build, test, and style code without being told each time.

To update these instructions, edit ``AGENTS.md``.
Keep this file under 300 lines to preserve LLM context.


Skills
------

WarpX defines reusable *skills* in the ``.claude/skills/`` directory.
Skills are scripted workflows that an LLM assistant can execute on demand, automating multi-step tasks that follow a fixed procedure.

All WarpX skills use the prefix ``warpx-`` for easy discovery (start typing ``/warpx`` to see them).

Currently available skills:

``/warpx-answer-user-question``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Drafts a response to a WarpX user question from a GitHub issue, discussion, or email.

Usage (in Claude Code):

.. code-block:: text

   /warpx-answer-user-question https://github.com/BLAST-WarpX/warpx/discussions/1234

The skill will:

#. Fetch the question from the provided URL (or use a pasted question directly).
#. Categorize the question (installation, input parameters, diagnostics, physics, etc.).
#. Search the WarpX source code, documentation, and past issues for relevant information.
#. Draft a response in the style of an experienced WarpX developer.
#. Present the draft for review before posting.

``/warpx-new-paper-highlight``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adds a new paper to the WarpX science highlights documentation (``Docs/source/highlights.rst``).

Usage (in Claude Code):

.. code-block:: text

   /warpx-new-paper-highlight https://doi.org/10.1103/PhysRevLett.133.045002

The skill will:

#. Fetch the paper metadata (authors, title, journal, DOI) from the provided URL.
#. Choose the appropriate highlights section (e.g., Plasma-Based Acceleration, HPC and Numerics).
#. Format the entry in the RST style used in the file.
#. Create a branch, commit the change, and optionally open a pull request.

``/warpx-self-review-pr``
^^^^^^^^^^^^^^^^^^^^^^^^^

Reviews the changes on your current branch relative to ``development``, as a WarpX maintainer would review a pull request (see the `Self-Reviewing a Pull Request`_ section below).

Usage (in Claude Code):

.. code-block:: text

   /warpx-self-review-pr

The skill will:

#. Read ``AGENTS.md`` for the project conventions.
#. Fetch the latest upstream ``development`` and diff your branch against it.
#. Review the diff against a checklist (correctness, algorithmic scaling, dimensionality, GPU/CPU portability, backward compatibility, testing, style, auto-generated files, documentation, scope).
#. Report concrete findings with file and line references, each with a severity and a suggested fix.

To add new skills, create a directory under ``.claude/skills/<skill-name>/`` containing a ``SKILL.md`` file that describes the step-by-step procedure.


Self-Reviewing a Pull Request
-----------------------------

As emphasized above, you should carefully review all LLM-generated code *yourself* before requesting a review from other WarpX developers.
A coding assistant can help with this first pass: it can re-read your own diff with fresh context and flag issues you may have missed.
This does not replace your own critical review, but it makes that review more effective.

Before using the prompt below, commit your work and make sure your branch is up to date with ``development``, so that the diff the assistant reviews matches what reviewers will see.
Run it in a *fresh* session (not the one that wrote the code): an assistant asked to review its own work in the same session tends to defend earlier choices rather than scrutinize them.

The review is a static reading of the diff: the assistant is asked not to compile the code or run the test suite, since the CI checks do that once the pull request is open.
This keeps the pass fast and avoids spending a long local build on something CI reports anyway.
The prompt says so explicitly because ``AGENTS.md`` documents how to build and test WarpX, and an assistant reading it for the project conventions would otherwise be inclined to do both.

In Claude Code, the ``/warpx-self-review-pr`` skill runs this review for you.
For other assistants, copy and adapt the following prompt:

.. The checklist in the code-block below is duplicated verbatim in
   .claude/skills/warpx-self-review-pr/SKILL.md (which drives Claude Code).
   Keep both copies in sync.

.. code-block:: text

   Review the changes on my current branch relative to the `development`
   branch, as if you were a WarpX maintainer reviewing my pull request.
   Do not make any changes yet; report your findings first.

   This is a read-only review of the source. Do not configure or compile
   the code (no `cmake`, no `pip install`) and do not run the test suite
   (no `ctest`, no analysis scripts): compilation, tests, and checksums are
   covered by the CI checks that run once the pull request is open. Judge
   the diff by reading it, not by executing it. The `Build Commands` and
   `Testing` sections of AGENTS.md describe development work in general and
   do not apply to this review: ignore them here.

   Start by reading AGENTS.md for the project conventions (style,
   portability, dimensionality, backward compatibility), then run
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
   6. Testing: is there a test covering the new feature or bug fix? Judging
      from the input file and analysis script alone (without running it),
      does it look fast enough for a 2-core CI runner and written portably?
   7. Style: does the diff follow the C++/Python style in AGENTS.md, and
      does it avoid reformatting unrelated code?
   8. Auto-generated files: flag any manual edits to `.pyi` stubs,
      `dependencies.json`, or `Regression/Checksum/benchmarks_json/*.json`.
   9. Documentation: are new user-facing parameters or features documented?
   10. Scope: is anything unrelated to the stated purpose of the PR included?

   For each finding, state the severity (blocking / should-fix / nit) and
   suggest a concrete fix. End with a short summary of whether this PR looks
   ready for human review.

Treat the assistant's output as a checklist of things to verify yourself, not as a verdict: it may raise false positives or miss real problems.
Confirm each finding against the code before acting on it.


Documentation Context via MCP Servers
--------------------------------------

LLM assistants work best when they can query up-to-date project documentation.
The :ref:`AI-Assisted Input File Design <ai_input_design>` workflow page describes how to set up `Model Context Protocol (MCP) <https://modelcontextprotocol.io>`__ servers for this purpose.
That setup is equally useful for development tasks: the same documentation context that helps write input files also helps an assistant understand WarpX internals, AMReX and pyAMReX APIs, and conventions when writing C++ or Python code.

See :ref:`ai_input_design` for general MCP setup instructions with Context7.

WarpX and Dependency Documentation on Context7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since WarpX builds on top of `AMReX <https://amrex-codes.github.io/amrex/>`__, `pyAMReX <https://pyamrex.readthedocs.io>`__, and `openPMD-api <https://openpmd-api.readthedocs.io>`__, providing documentation for these dependencies alongside WarpX documentation gives the assistant much richer context for development tasks.

The following documentation is available through Context7:

- **WarpX**: `context7.com/blast-warpx/warpx <https://context7.com/blast-warpx/warpx>`__
- **AMReX**: `context7.com/amrex-codes/amrex <https://context7.com/amrex-codes/amrex>`__
- **pyAMReX**: `context7.com/amrex-codes/pyamrex <https://context7.com/amrex-codes/pyamrex>`__
- **openPMD-api**: `context7.com/openpmd/openpmd-api <https://context7.com/openpmd/openpmd-api>`__
- **openPMD-viewer**: `context7.com/openpmd/openpmd-viewer <https://context7.com/openpmd/openpmd-viewer>`__
- **PICSAR-QED**: `context7.com/ecp-warpx/picsar <https://context7.com/ecp-warpx/picsar>`__
- **PICMI**: `context7.com/picmi-standard/picmi <https://context7.com/picmi-standard/picmi>`__
- **pybind11**: `context7.com/pybind/pybind11 <https://context7.com/pybind/pybind11>`__

When Context7 connected, the assistant can look up any of those when needed:
AMReX data structures (e.g., ``MultiFab``, ``ParticleContainer``, ``Geometry``), pyAMReX and pybind11 binding patterns, and openPMD I/O APIs directly, which is especially helpful when working, for instance, on:

- Field solver implementations that use AMReX mesh data structures
- Particle routines built on ``amrex::ParticleContainer``
- Python bindings that wrap C++ classes via pybind11 and pyAMReX
- I/O and diagnostic code that interacts with AMReX plotfiles or openPMD

For instructions on configuring Context7 as an MCP server in your coding assistant (Claude Code, Cursor, VS Code, Windsurf, Codex CLI, and others), see the `Context7 client documentation <https://context7.com/docs/resources/all-clients>`__ and the :ref:`ai_input_design` page.
