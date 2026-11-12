Anchor-Probe Sketching (APS): Implementation Specification
Status: proposed research prototype.

0. Goal and scope
Implement a standalone command-line tool that builds a compact sketch of each assembled genome (multi-contig FASTA or multi read FASTQ supported), compares sketches without accessing FASTA again, and reports:
1. Global k-mer-set Jaccard similarity.
2. Directional k-mer-set containment, A in B and B in A.
3. Local-homology-aware approximate average nucleotide identity (ANI).
4. Estimated directional alignment fraction (AF), A and B separately.

MVP must not implement base-level alignment or DP chaining. Use deterministic FracMinHash sampling, hash joins, ordered arrays, and simple collinearity checks. Report uncertain ANI/AF as -1, with diagnostic coverage counters, rather than fabricating confident numbers. Distinguish set similarity from nucleotide-level similarity everywhere.

1. Working definitions
- Anchor: selected canonical k1-mer, sampled genome-wide with hash h1 at probability p1=1/c1. Anchors locate correspondence between genomes. Keep occurrence positions and orientation.
- Probe: selected canonical k2-mer, sampled genome-wide with independent hash h2 at probability p2=1/c2. Probes measure k-mer overlap and local sequence similarity. Keep occurrence positions.
- Island: interval from one anchor to the next on the same contig. Each probe occurrence is owned by the last anchor at or before its start position, with a START sentinel for probes before the first anchor. A terminal END sentinel closes the last island.
- Shared anchor: same anchor hash present in both sketches. For MVP, use only hashes occurring exactly once in each whole genome for positional mapping.

2. Initial parameters
These are prototype defaults, not optimized values:
- k1=21, c1=125, hash seed h1 fixed and versioned.
- k2=15, c2=125, independent hash seed h2 fixed and versioned.
- if any other parameter needed, ask me about them

3. Hashing and sampling rules
For each A/C/G/T-only k-mer occurrence:
1. Compute the reverse complement.
2. Choose canonical representation canon=min(forward, reverse_complement) lexicographically (or fixed 2-bit encoding equivalent).
3. Record strand = 0 if forward is canonical, 1 if reverse complement is canonical. Odd k avoids strand ambiguity from perfect reverse-complement palindromes.
4. Calculate h1(canon) for k1-mers or h2(canon) for k2-mers. These are fixed, stable 64-bit hashes. Hashes must not change across processes, builds or genomes.
5. Select anchors iff h1 < floor(2^64 / c1); select probes iff h2 < floor(2^64 / c2). Beware integer overflow at 2^64: use a safe threshold calculation. c=1 selects everything.
6. Skip k-mers containing N or other ambiguous bases. Do not create k-mers across contig boundaries.
Important: sample probes globally, BEFORE assigning them to islands. Never reset/reseed probe sampling at an anchor. Keep every selected probe occurrence, including those overlapping anchors. Exclude selected anchor-overlapping probes only later when calculating local ANI (see section 8). K-mer hash collisions are presumed negligible for a 64-bit prototype, but retain a collision warning and make a 128-bit alternative possible later.


4. Sketch storage format
Come up with a format that complements but stands out from other storage formats.

5. BUILD algorithm (FASTA -> reusable sketch)
1. Parse input FASTA and record contig names/lengths. Uppercase bases. Any k-mer with ambiguous bases is excluded.
2. For each contig, scan all k1-mer starts and sample anchors with h1. Append (hash, pos, strand); positions increase.
3. Scan all k2-mer starts and sample probes with h2. Append (hash, pos, strand, role) in increasing position order.
4. Assign every probe to the immediately preceding real anchor by start position. Use a linear merge of the sorted anchor/probe arrays, not a per-probe binary search. For each anchor, store the half-open range of its owned probes. Keep pre-first-anchor probes under START. Alternatively, u can do this in a single pass.
5. Set every anchor's next_anchor_idx. END marks the last anchor.
6. Add every selected probe hash to a genome-wide list; sort and deduplicate for Jaccard and containment. Store this separately from positional occurrences.
7. Write sketch.
8. Print construction diagnostics: contigs, bases, anchor occurrences, probe occurrences, distinct global probe hashes, masked/ambiguous kmers, and file bytes.

Expected retained occurrence counts for long low-repeat genomes: anchors about L/c1; probes about L/c2. Distinct global probe count about D_k2/c2, where D_k2 is the number of distinct valid canonical k2-mers. These counts differ when a genome contains repeats. A probe belongs to exactly one preceding-anchor bucket but remains usable across several intermediate anchor buckets when comparing a distant shared-anchor pair.
BUILD pseudocode
function build_sketch(fasta, cfg):
    S = new Sketch(cfg)
    for seq in fasta.contigs:
        A = []
        P = []
        for (pos, kmer) in valid_kmers(seq, cfg.k1):
            canon, strand = canonicalize(kmer)
            hv = hash64(canon, cfg.seed1)
            if selected(hv, cfg.c1):
                A.append(Anchor(hv, pos, strand))
        for (pos, kmer) in valid_kmers(seq, cfg.k2):
            canon, strand = canonicalize(kmer)
            hv = hash64(canon, cfg.seed2)
            if selected(hv, cfg.c2):
                role = hash64(canon, cfg.seed3) & 1
                P.append(Probe(hv, pos, strand, role))
        assign_probe_ranges_by_linear_merge(A, P, START, END)
        S.contigs.append(Contig(seq.name, len(seq), A, P))
        S.global_probe_hashes.extend(probe.hash for probe in P)
    S.global_probe_hashes = sorted_unique(S.global_probe_hashes)
    save(S)
    
6. COMPARE: preliminary metrics from global probe sets
Require EXACTLY matching k2, c2, hash algorithm and seed2 (and canonicalization rules) across both sketches. Require matching anchor parameters for local matching. Never compare incompatible sketches silently.

Let SA and SB be sorted unique global probe hashes, a=|SA|, b=|SB|, x=|SA intersect SB|, u=a+b-x. Count x in O(a+b) by a two-pointer merge.
- Jaccard of distinct genomic k2-mer sets: J_hat = x/u, when u>0; else null.
- Containment of A in B: C_AinB = x/a, when a>0; else null.
- Containment of B in A: C_BinA = x/b, when b>0; else null.
- Max-containment (screening): C_max=x/min(a,b), if min(a,b)>0.
- Rough putative ANI for candidate screening, not final ANI: ANI_putative = C_max^(1/k2) only under a strong homologous-contained-sequence/no-large-repeats model. Do not assert it is valid for incomplete or rearranged genomes.
Use the whole global set for J/C (both role buckets). These J/C estimators are ratio estimates and their uncertainty can be substantial for small sketches. An empty sample produces unknown/null, not zero similarity. Optional coarse screening can apply a lower h2 threshold to the already-selected global hashes (nested FMH); do not add another sketch initially.

The approximate relationship for equal-length fully homologous sequences is q~ANI^k and J~q/(2-q), hence ANI_rough~(2J/(1+J))^(1/k). This is not a generic ANI estimator for incomplete genomes or unequal homologous fractions.
7. COMPARE: matching positional anchors and constructing candidate blocks
1. Build hash -> occurrence count for anchors of each genome.
2. Select anchor hashes with exactly ONE occurrence in A and ONE occurrence in B. Hash join to produce exact matches (contigA, posA, strandA, contigB, posB, strandB, hash).
3. Infer relative orientation from strandA XOR strandB. For reverse-orientation candidates, convert B anchor positions into forward-increasing oriented coordinates: y_oriented = contigB.length - (posB+k1). All subsequent gap checks use oriented coordinates.
4. Group by (contigA, contigB, relative_orientation). Sort each group by posA. Contig boundaries must never be crossed.
5. Scan each group in increasing posA. Form candidate pairs from consecutive shared anchors whose y_oriented coordinates also increase and whose inter-anchor gaps are positive and reasonably comparable (configured max_gap_ratio). Skip/break on large gaps, order reversals and implausible spacing; log reasons.
6. Greedily coalesce adjacent valid candidate pairs until the block span reaches min_block_span, or a contig/orientation/gap boundary occurs. Keep blocks disjoint within a genome. No dynamic programming.
7. For each block save its coordinate intervals [left_anchor.pos, right_anchor.pos+k1) in A and the corresponding oriented interval in B. The interior for scoring excludes anchor-overlapping k-mers.
Notes:
- Shared anchors c,d need not be adjacent in either genome's original anchor list.
- A simple greedy ordered scan is deliberately conservative. It may miss valid inversions, duplicated regions and ambiguous rearrangements.
- Matching only unique anchor hashes is a correctness-first MVP policy, not a final solution for repeated genomic elements.
- For sketches with many contigs, coordinate (contig_id, pos) is mandatory; a global concatenated coordinate alone is not enough.
MATCH pseudocode
function candidate_blocks(SA, SB, cfg):
    unique_A = anchor_hash_to_occurrence_if_count_one(SA)
    unique_B = anchor_hash_to_occurrence_if_count_one(SB)
    shared = hash_join(unique_A, unique_B)
    groups = group_by_contig_pair_and_relative_orientation(shared)
    blocks = []
    for group in groups:
        hits = sort_by_query_pos(group)
        run = []
        for neighboring pair (left, right) in hits:
            gapA = right.posA - left.posA
            gapB = right.y_oriented - left.y_oriented
            if gapA <= 0 or gapB <= 0:
                close_run(run); run = []; continue
            if max(gapA,gapB) > cfg.max_anchor_gap:
                close_run(run); run = []; continue
            if max(gapA,gapB) / min(gapA,gapB) > cfg.max_gap_ratio:
                close_run(run); run = []; continue
            extend_run_with_pair(run, left, right)
            if run.length >= cfg.min_block_span:
                blocks.append(run_as_block(run)); run = []
        close_run(run)  # include remainder only if enough evidence later
    return resolve_overlaps_conservatively(blocks)
For a production implementation, ensure close_run() does not discard a valid trailing run without either emitting it or explicitly reporting it as unresolved. Avoid assigning the same genomic bases to multiple accepted blocks.
8. COMPARE: validating a block using MAPPING probes
For each candidate block:
1. Retrieve all probes between the left and right anchors on each contig. This can be a contiguous slice from the position-ordered probe arrays; it is equivalent to following next-anchor pointers across all intermediate islands. Binary search in a sorted probe array is allowed to find the two endpoints, but avoid building large temporary copies.
2. Select only role=MAPPING probes. For positional matching, mask probe hashes appearing more than once inside either block (or, for stricter mode, more than once in either whole genome). This avoids one-to-many positional mappings.
3. Build a hash -> positionB table for unique B mapping probes. Scan A mapping probes in genomic order and look up matching hashes. Record pairs (xA, yB_oriented). For probes k2 on reverse B, use yB_oriented = contigB.length - (posB+k2).
4. Count mapping hits H, fraction of hits that preserve monotone order, and coarse displacement consistency (y-x). A simple and fast baseline is: sort hits by xA (already scanning in order); count adjacent y values increasing; identify a dominant diagonal bin only when indel variation is small. Do not assume a single fixed diagonal across long intervals with indels.
5. Accept only when enough matching probes support a plausible local mapping, or when adjacent shared anchor pairs themselves provide strong unambiguous positional context. Otherwise mark unresolved. Tune min_mapping_hits and diagonal/order cutoffs from simulations.
6. Optionally, segment a long block into shorter intervals based on rolling probe-match support. This avoids calling a mixed homologous/novel block completely homologous.
7. Never call all-zero mapping hits 'ANI=0'; it may simply be insufficient sampled evidence.
The MVP may use (H >= min_mapping_hits and >=80% of mapped hits preserve order), together with anchor uniqueness and length plausibility, as a provisional heuristic; make every threshold adjustable. It must not be described as a calibrated statistical significance test. A single matching probe or a perfect match out of 1-2 probes is not enough to declare homology.

For very small intervals, merge adjacent anchor intervals BEFORE applying the mapping test; do not insist every 150-bp island have several matching probes. For regions around 80% ANI, exact 15-mer survival is only 0.8^15~0.035, so sparse sampling is noisy even in true homologs.
9. COMPARE: ANI from independent SCORING probes
Use only role=SCORING probes within accepted homologous blocks. Using a different deterministic subset from MAPPING reduces upward bias caused by selecting blocks using their own probe-match scores. Independence is approximate: sequence heterogeneity and conditioning on shared anchors still matter.
1. Exclude scoring probe occurrences whose k2-mer spans overlap any shared anchor occurrence used to delimit/establish the accepted block. These anchors were selected because they match exactly, so overlapping probe matches would be artificially enriched.
2. Exclude ambiguous high-copy/low-complexity scoring k-mers consistently in both genomes. Keep explicit counts of excluded probes. For a minimal initial implementation, ignore scoring hashes with multiple occurrences inside either accepted block. When counting matches, use a hash set or sort hashes explicitly; the stored probe arrays are ordered by position, not by hash.
3. Let nA be count of remaining eligible SCORING probe hashes/occurrences in A's accepted block; nB analogous in B. Let x be the number of uniquely matched scoring hashes between the two corresponding accepted regions. For the MVP, unique hashes within each block make x a straightforward set intersection.
4. Compute directional local qA=x/nA and qB=x/nB when denominators nonzero. For one symmetric block statistic compute q=2*x/(nA+nB). Then ANI_block=q^(1/k2). This is a model-based approximation to ANI, not an exact nucleotide alignment identity.
5. Do NOT average block ANI values unweighted. Sum all eligible scoring counts and matches over accepted disjoint blocks: XA=sum(x_block), NA=sum(nA_block), NB=sum(nB_block). Then q_global=2*XA/(NA+NB) and ANI_hat=q_global^(1/k2). Return directional q and ANI estimates as optional diagnostics.
6. If NA+NB is too small, or accepted aligned-region support is inadequate, return ANI=null with explicit status and retained count fields. Never return zero merely because no exact k2-mer matches are seen.
Model assumptions: homologous sequence, negligible accidental matches and repeats, approximately independent substitutions, and minor effects of indels. In reality, indels destroy overlapping k-mers, recombination creates heterogeneous ANI, and sparse region selection can bias the estimate. Treat the reported ANI as uncalibrated until benchmarked against base-level alignments.
Score pseudocode
function score_ani(accepted_blocks, SA, SB):
    X = NA = NB = 0
    for block in accepted_blocks:
        PA = eligible_scoring_probes(SA, block.regionA,
                  exclude_shared_anchor_overlap=true,
                  exclude_repetitive=true)
        PB = eligible_scoring_probes(SB, block.regionB,
                  exclude_shared_anchor_overlap=true,
                  exclude_repetitive=true)
        na = len(PA)
        nb = len(PB)
        x = hash_set_intersection_count(PA, PB)  # arrays are sorted by POSITION, not hash
        X += x; NA += na; NB += nb
    if NA+NB == 0 or not enough_support(X,NA,NB,accepted_blocks):
        return null, {X,NA,NB,status:'insufficient_evidence'}
    q = 2*X/(NA+NB)
    return q**(1.0/k2), {X,NA,NB,q,status:'estimated'}
If desired later: compute ANI for fixed genomic windows then aggregate (skani uses query chunks); windowing allows heterogeneity analysis and optional bootstrap. For MVP, pool counts directly.
10. COMPARE: estimated alignment fraction (AF)
AF is NOT directly recoverable from global unordered FracMinHash sets. It must be estimated from positional candidate blocks, and the MVP will output a proxy rather than a validated base-level alignment fraction.
1. Map each ACCEPTED homologous block to a half-open base interval in A and a corresponding interval in B. Each must be on a single contig. Avoid extending blocks across unsupported long gaps.
2. Sort the intervals per genome and merge overlaps; sum their union lengths: coveredA and coveredB.
3. Define AF_span_A = coveredA / total_contig_length_A; AF_span_B = coveredB / total_contig_length_B. Include N-containing bases in denominators consistently with contig lengths. This is a conservative-to-optimistic proxy depending on terminal coverage and whether large accepted blocks contain internal novel regions: it is NOT guaranteed to be a lower or upper bound.
4. A block passing an order check but containing a long internal nonhomologous segment can overestimate AF. Optional local window segmentation should eventually exclude unsupported sections. Conversely, segments outside the outermost shared anchors will be uncounted and can underestimate AF.
5. Also report a purely sketch-based cross-check, AF_containment_A ~= C_AinB / (ANI_hat^k2), and analogously for B, under a simplified model where fraction f of k-mers is truly homologous and matched fraction among homologous k-mers is ANI^k2. This estimates k-mer homologous fraction, not true aligned-base fraction; it can exceed 1 or be biased under repeats/incompleteness. Return as a diagnostic, preferably unclipped plus warning, never replace AF_span silently.
6. Report AF_span_A, AF_span_B separately; optional AF_mean=(AF_A+AF_B)/2. Do not use Jaccard or a single containment as AF.
7. If no trustworthy blocks are established, set AF_span to null with status 'no_supported_mapping', not to an asserted zero. Optionally also report raw accepted_span_bases=0 as an observable.
11. COMPLETE COMPARE pseudocode
function compare(sketchA, sketchB, cfg):
    verify_compatible(sketchA, sketchB)

    # Always available if global probe sketches are nonempty.
    x_global = sorted_intersection_size(
        sketchA.global_probe_hashes, sketchB.global_probe_hashes)
    a = len(sketchA.global_probe_hashes)
    b = len(sketchB.global_probe_hashes)
    J = safe_div(x_global, a+b-x_global)
    C_A = safe_div(x_global, a)
    C_B = safe_div(x_global, b)
    putative_ANI = safe_pow(safe_div(x_global, min(a,b)), 1.0/k2)

    # Positional matching from anchors.
    candidates = candidate_blocks(sketchA, sketchB, cfg)
    accepted = []
    diagnostics = Counter()
    for block in candidates:
        mapping_info = validate_block_using_mapping_probes(block, sketchA, sketchB, cfg)
        diagnostics.add(mapping_info)
        if mapping_info.status == 'supported':
            accepted.append(block_or_supported_subwindows(block,mapping_info))

    ANI, ani_diag = score_ani(accepted, sketchA, sketchB)
    AF_A, AF_B, af_diag = estimate_span_AF(accepted, sketchA, sketchB)
    AF_C_A = C_A / ANI**k2 if ANI is valid and ANI>0 else null
    AF_C_B = C_B / ANI**k2 if ANI is valid and ANI>0 else null

    return {
        jaccard: J,
        containment_A_in_B: C_A,
        containment_B_in_A: C_B,
        putative_ani: putative_ANI,
        ani: ANI,
        af_span_A: AF_A,
        af_span_B: AF_B,
        af_containment_proxy_A: AF_C_A,
        af_containment_proxy_B: AF_C_B,
        counts: {global_a:a, global_b:b, global_shared:x_global,
                 candidate_blocks:len(candidates),
                 accepted_blocks:len(accepted), ...ani_diag, ...af_diag},
        status: {jaccard:..., containment:..., ani:..., af:...},
        warnings: [...]
    }
12. Complexity targets
At build time, scan O(L) k-mer positions with rolling hashes/2-bit encoding, sort global distinct probe hashes once, and write O(L/c1 + L/c2) positional records in typical low-repeat genomes (with all retained positions counted). Hashing a full k-mer string every step is acceptable for a correctness prototype but should be replaced by rolling packed encoding for speed.

At pairwise compare time, global J/C uses O(D_A+D_B) where D are distinct stored probe counts; anchor hash join O(A_A+A_B) expected; each position-ordered probe should ideally be visited at most once across disjoint candidate blocks, yielding O(P_A+P_B) expected processing. Sorting matched anchor groups can cost O(M log M); building many candidate pairs, processing repeated regions or rescanning overlapping blocks can exceed linear time. Report actual per-stage durations and counts rather than asserting formal worst-case linearity.

Keep memory light: arrays, compact binary records, no full DNA strings, no per-base alignment matrix, no full seed-to-seed DP table. Hash tables should be local to one candidate block and released before the next.
13. CLI contract
conform to existing

14. Testing checklist (implement before optimization)
Unit tests:
- FASTA with single and multiple contigs; ambiguous bases N; windows never cross contig ends.
- Canonical k-mer representation invariant under reverse complement.
- Deterministic hash and sampling across independent runs; serialization round-trip exact.
- All selected probes assigned exactly once, including START/END buckets.
- Reconstruct probe slice between nonadjacent anchors across several internal anchor buckets; ensure nothing is omitted or counted twice.
- Same data from both files: J=1, C_A=C_B=1 (if sample nonempty); ANI near 1 and AF near 1 when adequate anchor support exists. Terminal gaps may give AF below 1; this is expected.
- Completely different random sequences: J and containment near 0 with sufficient sample size; ANI/AF null without homologous blocks.
- A proper subsequence (with enough unique shared anchors): C_sub_in_full near 1, reverse containment much smaller; J should reflect the length difference. AF_sub should be higher than AF_full. Identity of matched regions near 1.
- Reverse-complement of a contig: same global J/C and positional ANI/AF after coordinate transformation.
- Destroy one intermediate anchor while preserving two outer anchors: candidate must still be retrievable and include probes across the lost anchor.
- Insert a new anchor between two conserved anchors: candidate and probe retrieval still succeed.
- Repeated anchor hashes: ignored for positional matching in MVP; no crash, status/warning logged.
- Nonmatching contig / insertion / rearrangement / inversion: conservative block rejection or segmentation; do not return artificially precise ANI/AF.
- Small or empty samples: null with status; never NaN/inf/divide by zero.
- Distinct kmer set counts vs duplicate occurrences: global J/C use unique hashes, positional probes use occurrences.
Simulation tests:
- Generate a sequence with configurable substitutions (ANI 0.80, 0.85, 0.90, 0.95, 0.99, 0.999).
- Independently test insertion/deletion rates, 10%-50% sequence removal, fragmented contigs, duplication, and inversion.
- Compare estimated J/C with exact full-kmer-set statistics and ANI/AF with ground-truth base-level alignments.
- At each setting track candidate anchor coverage, block acceptance, ANI bias/MSE, AF bias/MSE, zero-probe probability, runtime and memory.
- Sweep k1, c1, k2, c2, block size and mapping-hit threshold; report failure/null rate as well as accuracy conditional on reporting.
15. Build order for coding agent
Milestone 1: Global FMH correctness. FASTA reader, canonical k-mer hashing, one global probe sketch, persisted file, Jaccard + directional containment; tests vs exact kmer sets on small sequences.

Milestone 2: Positional sketch format. Anchor sampling, sorted anchor/probe occurrences, per-anchor ranges and next pointers, START/END sentinels; test nonadjacent anchor retrieval.

Milestone 3: Lightweight anchor mapping. Unique anchor matching, orientation, sorted candidate intervals, gap filtering, nonoverlap resolution; tests for lost/new intermediate anchors and reverse-complement contigs.

Milestone 4: Local positional validation. Mapping probe hash join, order/gap checks, accepted/unresolved blocks, diagnostics; keep deliberately conservative.

Milestone 5: ANI and AF. Independent scoring-bucket probes, boundary-anchor-overlap masking, pooled ANI estimate, merged accepted-span AF and optional global-containment AF cross-check. Return null when unsupported.

Milestone 6: Benchmark and tune. Simulations plus skani comparisons, sensitivity/accuracy/runtime/memory; only then decide whether to add diagonal voting, sliding-window segmentation, faster hashing or C++/Rust rewrite.
16. Non-goals for MVP and common mistakes
- Do not implement Smith-Waterman, BLAST-like alignment or dynamic-programming chaining.
- Do not implement a learned ANI correction initially.
- Do not derive AF directly from Jaccard, or confuse k-mer containment with sequence alignment fraction.
- Do not use sampled-position occurrence counts as though they were distinct k-mer counts for global J/C.
- Do not use different hash seeds, k sizes or thresholds between sketch files.
- Do not accept all-probe-match regions with just 1-2 probes as proven homologous.
- Do not reject true 80%-ANI regions because a tiny window contains zero matching k2-mers.
- Do not treat nonmatches due to genome incompleteness as substitutions for ANI estimation.
- Do not use scoring probes to select the same region they are scored on (mapping/scoring split).
- Do not double count overlapping candidate intervals.
- Do not hide uncertainty: return diagnostics and null values when mapping evidence is insufficient.
17. References motivating estimators (not claiming novelty)
- Shaw and Yu (2023), 'Fast and robust metagenomic sequence comparison through sparse chaining with skani', Nature Methods. https://www.nature.com/articles/s41592-023-02018-3
- 'Estimating similarity and distance using FracMinHash', https://pmc.ncbi.nlm.nih.gov/articles/PMC12082993/
- 'Deriving confidence intervals for mutation rates across a wide range of evolutionary distances using FracMinHash', https://genome.cshlp.org/content/33/7/1061
This document is an implementable research baseline. In particular, local ANI is a model-based estimate, and AF_span is an empirical positional-coverage proxy. Do not describe either as fully validated against alignment truth until experiments establish accuracy.