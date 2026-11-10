# Complete the placeholder sketch/compare CLI

## Scope

Replace the top-level sequence-summary command with `sketch` and `compare`.
Reuse streaming FASTA/FASTQ parsing and existing k-mer semantics. This stage creates
empty placeholder sketches and unavailable metric values; actual sketching,
serialization, and estimation are future work.

The user explicitly approved removal of the legacy command and stdin support.
This is an authorized exception to AGENTS.md's general CLI-preservation rule;
the subcommand interface is the compatibility baseline going forward.

## Command contract

```text
aiips sketch INPUT -o OUTPUT [-k LENGTH] [-f]
aiips compare INPUT_A INPUT_B [-k LENGTH]
```

- The old `aiips INPUT...` interface is removed. Missing/unknown commands and
  invalid arguments fail with usage guidance and exit status 1.
- Top-level and command-specific `-h`/`--help` print help and exit 0. For a valid
  subcommand, help before `--` takes precedence over other invalid arguments.
- `--` ends option parsing for positional filenames. Only filesystem paths are
  accepted; stdin (`-`) is rejected. Output values beginning with `-` must use an
  explicit path prefix, such as `./`, to distinguish them from missing values.
- `-k`/`--kmer-size` defaults to 21 and accepts positive values fitting `usize`.
  Repeating either spelling, including mixed spellings, is an error. Help before
  `--` still takes precedence over that error.
- No speculative algorithm options or metric-selection flags are introduced.

## Sketch

- Accept exactly one plain/gzip FASTA/FASTQ file, detected by content.
- Require `-o`/`--output`. Append `.aiips` unless the final extension already
  matches case-insensitively. Print the actual destination on success.
- Treat the dotfile `.aiips` as already suffixed. Reject spellings ending in a
  separator, `.` or `..`, before appending any suffix. Reject a directory at the
  final (suffixed) destination; a directory named like the unsuffixed spelling
  does not block it. Directory errors must not suggest using force, and no error
  suggests force when force was already given.
- Recognize every separator accepted by the target platform (both `/` and `\`
  on Windows) when rejecting directory-like output spellings.
- Preflight output existence, type, aliases, and destination-parent validity
  before parsing sequence content. Recheck output safeguards before installing
  the result; preflight must not create output or bypass exclusive creation.
- Fully parse and summarize the input before creating a zero-byte output.
- Without force, create the empty output with atomic exclusive file creation;
  require no hard-link support and never fall back to a replacing rename.
- Refuse existing output by default, explaining `-f`/`--force` in the error.
- Support force only on `sketch`: stage replacement in the destination directory
  and rename after successful validation. Preserve existing output on parsing or
  replacement failure; clean up staging files.
- Reject output resolving to the input by Unix device/inode identity, including
  case aliases, hardlinks, and symlinks. Do not create
  destination directories. Forced replacement replaces symlinks themselves.
- Preserve permission bits when replacing a regular file. The new inode does not
  preserve ownership, group, ACLs, xattrs, or timestamps; replacing a symlink does
  not copy its target's metadata. Refuse forced replacement of existing paths on
  non-Unix platforms until equivalent file-identity safeguards are implemented.
- Keep an original-output backup until stdout writes and flush succeed. On a
  reported output failure, remove newly created output or restore the original
  file/symlink. If restoration fails, retain the original backup and report its
  path. Partial stdout cannot be retracted. Abrupt termination/crash recovery and
  concurrent external path changes are outside this rollback guarantee.
- Once stdout writes and flush succeed, commit the new output. If deleting the
  old backup then fails, keep both files, warn on stderr with the backup path,
  and exit successfully; do not restore old data over reported-success output.
- Print the statistics table followed by `output<TAB>actual-path` to stdout.

## Compare

- Accept exactly two inputs, detecting each independently. A final `.aiips`
  extension or whole filename, case-insensitively, identifies a sketch. All other
  inputs are FASTX. Reuse one metadata lookup for type and placeholder-length
  checks, rather than statting a sketch twice.
- Support FASTX/FASTX, sketch/sketch, and mixed pairs.
- Sketch inputs must be readable zero-byte regular files. Reject nonempty files
  until a real format is implemented.
- Validate both inputs before printing output. Print one statistics row per input;
  use `NA` for all five statistic fields, including `k`, for sketches.
- After a blank line, print a comparison table with these columns:

```text
input_a input_b jaccard containment_a_in_b containment_b_in_a ani af_a_to_b af_b_to_a
```

- Fields are tab-separated. All six metrics are `-1`, meaning unavailable.
  Containment and AF are reported in both directions; Jaccard and ANI have one
  column each. No real estimator or statistical assumptions are introduced.
- `compare` accepts neither output-path nor force options.

## Shared output and errors

Statistics columns remain `input`, `records`, `bases`, `non_acgt_bases`, `k`, and
`valid_kmers`. Count overlapping occurrences, including repeats. A/C/G/T in either
case are valid; other bases break runs, and records never share windows. Qualities
are validated for length but do not filter k-mers.

Success output goes to stdout; errors go to stderr and identify the relevant path
where applicable. Return status 0 on success and 1 on error. Validate inputs before
printing success output. Empty/malformed sequence files, missing files, and
non-regular inputs fail. Escape backslash, tab, carriage return, and newline in
TSV filenames. Keep the dependency set unchanged.
Every stdout write/flush error, including compare statistics, names `stdout:`.

## Acceptance checks

- Retain existing library unit tests and the exhaustive naive k-mer baseline.
- Add at least 20 distinct integration scenarios per command, covering valid
  formats, gzip, known counts, argument errors, missing/malformed files, help,
  input combinations, placeholder validation, output generation, suffix behavior,
  overwrite refusal, force, preservation on failure, and alias protection.
- Test legacy interface rejection. Use Cargo's test executable and isolated
  temporary directories; never alter the downloaded genome or existing results.
- Run locked Cargo tests, formatting checks, Clippy, a release build, and
  representative command invocations. Update README examples and limitations.
- Cover every CODE_REVIEW.txt finding: aliases, suffix edge cases, directory
  diagnostics, permissions, stdout write/flush rollback, help precedence, sketch
  validation, output-option validation, and resilient temporary workspaces.
  Workspace creation must retry collisions; cleanup must not panic during unwind.
- Both unit and integration filesystem tests use the shared isolated workspace
  helper, so retained backups and directories are cleaned up even on assertion
  failure. Inject backup-cleanup and stdout failures for deterministic coverage.

`.aiips` identifies an empty placeholder only and establishes no future binary
format. Real sketch metadata and estimator definitions remain out of scope.
