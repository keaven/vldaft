# Agent Instructions

The authoritative project guide is `inst/site/llms.txt`. Read it before
making package changes; it describes the model, layout, build system,
backend contracts, and validation expectations.

Draft repo-local skills live under `.agents/skills/`:

- `$vldaft-development` for R/C/Rust package implementation, backend parity,
  tests, and package checks.
- `$vldaft-model-development` for turning statistical survival-model prompts
  into implementation plans, backend changes, simulation checks, and docs.

Planned:

- `$vldaft-docs-site` for roxygen, Rd, README, vignettes, pkgdown, and
  rendered-site verification. Add this skill before relying on the name.

Keep this file short. Put durable project guidance in `inst/site/llms.txt`,
and put task-specific reusable workflows in the relevant skill.
