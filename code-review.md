# Code review

1. **[P2] Avoid overflow in the sketch read limit** — `src/commands.rs:45`

   `max_file_bytes` saturates to `u64::MAX` for sufficiently large declared FracMinHash sizes, but `read_sketch` then evaluates `limit + 1` using ordinary addition. A tiny file declaring `size=18446744073709551615` reproducibly panics with exit 101 in the debug CLI instead of reporting its invalid value count. In release builds the addition wraps to zero and skips reading any remaining body. Use checked or saturating arithmetic consistently and add a CLI regression for a saturated size bound; the existing direct-parser test does not exercise this reader.

2. **[P2] Describe bias under the implemented empty-sketch policy** — `README.md:249–253`

   The statement that the reported raw containment is biased downward by `1 - (1 - 1/S)^|A|` applies the wrong conditioning. [Hera et al. (2023), Theorem 1](https://genome.cshlp.org/content/33/7/1061), analyze an estimator multiplied by an indicator for a nonempty query sketch, effectively assigning zero to empty-query outcomes. This implementation instead rejects comparisons whenever either sketch is empty. The downward-bias formula therefore does not describe the values it reports.

   Exact enumeration under the ideal independent-retention model illustrates the difference: for A={x,y}, B={x,z}, and retention probability 1/2, true containment is 0.5. Assigning zero to empty-query outcomes gives mean 0.375, consistent with the paper's factor. Conditioning on both sketches being nonempty, as the code does, gives mean 0.6. This conditional calculation is a review derivation, not a claim from the paper. Revise the README to distinguish these policies rather than asserting downward bias for reported comparisons; align the estimator documentation if needed. This follows the requirement to state statistical assumptions and ambiguities ([AGENTS.md:20–21](AGENTS.md#scientific-correctness)).

Validation: all 179 tests passed with `cargo test --locked --quiet` outside the sandbox, and all smoke tests passed. The malformed-size panic and the eight-outcome bias example were reproduced separately. The earlier MinHash scientific-documentation finding is resolved. Only this review document was added; source code and experimental results were unchanged. Other platforms were not tested.
