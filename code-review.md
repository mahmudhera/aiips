# Review: features/ci-cd against origin/main

Comparison base: `1461cdeb5b8e8084f4a4c2e1779dfd7e3147c6f8`.

## 1. [P1] Adapt existing tests before enabling the Windows job

Location: `.github/workflows/ci.yml:67–68`.

The new Windows job runs tests that necessarily fail there. In `tests/cli.rs`,
`sketch_force_replaces_and_cleans_staging` expects successful forced replacement,
which `commands.rs` explicitly refuses on non-Unix platforms. The unguarded
escaped-path tests also create filenames containing control characters prohibited
by [Windows filename rules](https://learn.microsoft.com/en-us/windows/win32/fileio/naming-a-file#naming-conventions).

Add platform-specific test expectations and guards while retaining Windows
coverage; otherwise every workflow run fails.

## 2. [P1] Preserve the output spelling in the Windows smoke test

Location: `scripts/smoke.sh:74–76`.

Passing an unsuffixed `out/sample` makes the CLI call `with_file_name`, which
inserts the platform separator according to the
[Rust implementation](https://doc.rust-lang.org/src/std/path.rs.html#1589-1595).
On Windows, `report::label` then escapes that backslash, so stdout differs from
the hard-coded forward-slash expectation. CRLF normalization cannot fix this.

Pass the already-suffixed `out/sample.aiips`, which preserves the supplied
spelling, or make the expected label platform-aware.

## 3. [P2] Remove the README link to the intermediate specification

Location: `README.md:176`.

This introduces a `SPEC.md` reference despite [AGENTS.md:8–10](AGENTS.md#referring),
which identifies specification files as branch intermediates and explicitly
prohibits referencing them in README.

Remove the link and keep the CI description self-contained.

## Validation

Local macOS validation passed: **124 tests**, release build, smoke tests,
formatting, and Clippy. The full test run required execution outside the sandbox
because an existing test creates a Unix socket.

Windows findings are based on static inspection; Windows execution was
unavailable. No implementation files were changed during the review.
