# CI specification: build-and-test workflow

## Goal
Every push to `main` must show that `aiips` builds, that all of its subcommands work,
and that all tests pass on Linux, macOS and Windows. A failure on any platform marks
the push as failed.

## Scope
- One workflow: `.github/workflows/ci.yml`.
- Triggers: `push` to `main` only. `workflow_dispatch` is allowed for manual reruns.
  Pull-request triggers, releases, artifact publishing and deployment are out of scope.

## Platforms
Matrix (`fail-fast: false`, so every platform reports its own result):
- `ubuntu-latest` (x86_64 Linux)
- `macos-latest` (arm64 macOS)
- `windows-latest` (x86_64 Windows, MSVC)

## Toolchain
- Install Rust 1.99.1 with rustup, plus `rustfmt` and `clippy`. This matches the pin in
  `environment.yml`. Do not use conda in CI.
- Every Cargo command uses `--locked`, so `Cargo.lock` is authoritative.
- Cache the Cargo registry and `target/`, keyed on OS, toolchain and `Cargo.lock`.
- No new crate dependencies. Third-party actions are pinned to a commit SHA or a full
  version tag.

## Required checks (each one is a failing step)
1. Lint (Linux only, since it does not depend on the platform):
   `cargo fmt --check` and `cargo clippy --locked --all-targets -- -D warnings`.
2. Build: `cargo build --locked --release` on every platform.
3. Tests: `cargo test --locked` on every platform, covering all unit, integration and
   doc tests. No test may be skipped or filtered in CI. Tests marked `#[cfg(unix)]` or
   `#[cfg(windows)]` run only on their target OS, which is expected.
4. Subcommand smoke tests on every platform, in `scripts/smoke.sh BINARY`. The script
   copies the committed fixtures into a temporary work directory as
   `tests/fixtures/`, runs from there, and passes only relative `/`-separated paths,
   so printed labels are the same on every platform:
   - `aiips --help`, `aiips sketch --help` and `aiips compare --help` exit 0.
   - `aiips sketch tests/fixtures/sample.fa -k 3 -o out/sample` exits 0. Its stdout
     is exactly:
     ```
     input	records	bases	non_acgt_bases	k	valid_kmers
     tests/fixtures/sample.fa	2	13	2	3	7
     output	out/sample.aiips
     ```
     It also creates an empty file `out/sample.aiips`. The counts are the known
     baseline from `tests/cli.rs` (`sketch_format!`).
   - `aiips compare tests/fixtures/sample.fa tests/fixtures/sample.fq.gz -k 3` exits 0.
     Both statistics rows show `2	13	2	3	7`, and all six metric values are `-1`.
   - `aiips compare out/sample.aiips tests/fixtures/sample.fa -k 3` exits 0. This is
     the sketch round trip, and the statistics row for the sketch is `NA` throughout.
   - Re-running the same `sketch` without `-f` exits 1, leaves `out/sample.aiips`
     unchanged, and prints an error that suggests `-f`.
   - An unknown command, such as `aiips bogus`, exits 1.
   - Compare stdout exactly against the expected text, after normalizing CRLF line
     endings on Windows.
   - Use the `bash` shell on every OS (available on Windows runners), so the same
     script runs everywhere.

## Scientific-correctness constraints
- Use only committed fixtures. Never modify `tests/fixtures/` or `data/`, and never
  write smoke-test output into the repository tree.
- Expected values come from the existing tests. Do not change a baseline to make CI
  pass.
- Tests stay deterministic (fixed seeds, no network access).

## Known risks and ambiguities
- Windows has never been tested. Forced replacement (`-f`) is refused on non-Unix
  platforms by design, as the README documents. If Windows exposes other failures,
  report them and fix them in the code. Do not drop Windows from the matrix without
  approval.
- On Windows, path labels in stdout may use `\` separators, depending on how the
  path is passed. The smoke script avoids this by passing only relative `/` paths.
- Git's line-ending conversion on Windows could change fixture bytes and break the
  known counts. CI disables `core.autocrlf` before checkout.
- `README.md` documents Unix paths (`./target/release/aiips`). On Windows the binary
  is `target/release/aiips.exe`.
- `macos-latest` is arm64 only. Intel macOS is not covered.

## Acceptance
- A push to `main` runs 3 platform jobs and 1 lint job, and all of them pass.
- Breaking the build, a test, formatting or Clippy, or any smoke-test assertion fails
  the workflow.
- When the workflow is added, the README gets a short CI section with a status badge
  and a summary of what CI checks.
