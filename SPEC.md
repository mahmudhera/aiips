## Spec name

In this spec/branch, we'll complete the cli. Hopefully this behavior will persist till the end of development.

## Sub-commands we want

- sketch: create and store a sketch. for now, simply create an empty placeholder file.
- compare: will compare two sketches, and estimate the similarity/distance. For now, we will again print -1 values. We plan to support the following metrics: Jaccard, containment, ANI, AF (alignment fraction).

Intended behavior: `aiips sketch` or `aiips compare`, etc. Both of the programs can, for now, read the input files and print the kmer statistics (the current behavior using `./target/release/aiips`)

The current behavior of the command: remove that.

## Supported parameters/arguments

- input files: for sketch, one. for compare, two
- k: length of kmers we will use in the sketches/comparisons
- any other param name

## Tests

- Add at least 10/20 tests for each of these commands: sketch and compare. Cover exiting gracefully on missing files, output file getting generated, etc. 

- Ensure the current tests are arewritten so that they no longer expect `target/release/aiips` to work.

## Update documentation

Update accordingly.