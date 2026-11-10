# Complete the placeholder sketch/compare CLI

## Scope

Replace the top-level sequence-summary command with `sketch` and `compare`.
Reuse streaming FASTA/FASTQ parsing and existing k-mer semantics. This stage creates
empty placeholder sketches and unavailable metric values; actual sketching,
serialization, and estimation are future work.

## Command contract

```text
aiips sketch INPUT -o OUTPUT [-k LENGTH] [-f]
aiips compare INPUT_A INPUT_B [-k LENGTH]
```

- The old `aiips INPUT...` interface is removed. Missing/unknown commands and
  invalid arguments fail with usage guidance and exit status 1.
- Top-level and command-specific `-h`/`--help` print help and exit 0.
- `--` ends option parsing for positional filenames. Only filesystem paths are
  accepted; stdin (`-`) is rejected. Output values beginning with `-` must use an
  explicit path prefix, such as `./`, to distinguish them from missing values.
- `-k`/`--kmer-size` defaults to 21 and accepts positive values fitting `usize`.
- No speculative algorithm options or metric-selection flags are introduced.

## Sketch

- Accept exactly one plain/gzip FASTA/FASTQ file, detected by content.
- Require `-o`/`--output`. Append `.aiips` unless the final extension already
  matches case-insensitively. Print the actual destination on success.
- Fully parse and summarize the input before creating a zero-byte output.
- Refuse existing output by default, explaining `-f`/`--force` in the error.
- Support force only on `sketch`: stage replacement in the destination directory
  and rename after successful validation. Preserve existing output on parsing or
  replacement failure; clean up staging files.
- Reject output resolving to the input, including symlink aliases. Do not create
  destination directories. Forced replacement replaces symlinks themselves.
- Print the statistics table followed by `output<TAB>actual-path` to stdout.

## Compare

- Accept exactly two inputs, detecting each independently. A final `.aiips`
  extension, case-insensitively, identifies a sketch. All other inputs are FASTX.
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

`.aiips` identifies an empty placeholder only and establishes no future binary
format. Real sketch metadata and estimator definitions remain out of scope.
