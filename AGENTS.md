# Project instructions

## Purpose
This is a scientific research codebase. But we ultimately
want something very fast and efficient.

## Referring
- within a branch, consider spec.md and code-review.md or similar files as intermediates
- These will not make it to the main branch
- Do not refer to spec.md in README, or to code-review.md

## Development rules
- Prefer simple, modular implementations.
- Do not introduce unnecessary dependencies.
- Preserve existing command-line interfaces if there is any.
- If using libraries, prefer header-only libraries.
- Implement things in lightweight manner. Should be both fast and memory efficient.

## Scientific correctness
- Do not invent formulas or statistical assumptions.
- Flag ambiguities in algorithm definitions.
- Use fixed random seeds in tests.
- Compare numerical outputs against known baselines.
- Do not modify ground-truth data to make tests pass.

## Verification
- Write tests for every new feature.
- Run relevant tests after each change.
- Report what passed, failed, or could not be tested.

## Workflow
- Do not modify unrelated files.
- Show a concise plan before significant changes.
- Summarize modified files and remaining limitations.
- Never overwrite experimental results without approval.
