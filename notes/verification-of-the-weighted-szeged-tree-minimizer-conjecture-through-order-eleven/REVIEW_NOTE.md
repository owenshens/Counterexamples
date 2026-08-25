# Verification of the Weighted Szeged Tree-Minimizer Conjecture Through Order Eleven

`verification-of-the-weighted-szeged-tree-minimizer-conjecture-through-order-eleven`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **83 checks, all passing**:

    VERDICT: ALL 83 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    18081259928ccf0846b5efe12127d098e7897b3c17d754324487acc728120567

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE: what this program re-ran and what it did not.
> SCOPE: exhaustive over EVERY isomorphism class of connected graph for n = 3,4,5,6,7,8 (generated here) and n = 9 (covering enumeration); both minima are therefore exact for n <= 9.
> SCOPE: exhaustive over EVERY free tree for n = 3..11, so the whole tree-minimum column of Table 1 is exact.
> SCOPE: for n = 10 the non-tree search was exhaustive only over connected graphs with at most 13 of the 45 possible edges (spanning tree plus k <= 4 extra edges).
> SCOPE: for n = 11 the non-tree search was exhaustive only over connected graphs with at most 13 of the 55 possible edges (spanning tree plus k <= 3 extra edges).
> SCOPE: NOT re-run: the part of the census with at least 14 edges at n = 10 and at least 14 edges at n = 11, which is the bulk of the 1018690327 classes and needs nauty geng; the class counts themselves are instead re-derived exactly by cycle index plus inverse Euler transform.
> SCOPE: the dense tail is partly closed by proof rather than by search.  Since n_u >= 1 and n_v >= 1 on every edge, wSz(G) >= sum_v deg(v)^2, which over degree sequences of sum 2m is minimised by the balanced sequence:
> SCOPE:   n = 10: every connected graph with m >= 39 edges has wSz >= 610 > 600, so the only edge counts neither searched nor bounded are 14 <= m <= 38 (of the full range 9..45).
> SCOPE:   n = 11: every connected graph with m >= 47 edges has wSz >= 806 > 778, so the only edge counts neither searched nor bounded are 14 <= m <= 46 (of the full range 10..55).
> SCOPE: the per-edge-count minima of the full n = 7 and n = 8 census are NOT monotone in m -- measured in this run, the strict descents are n=7: m 11->12 min 262->244, m 14->15 min 320->306, m 15->16 min 306->252, m 19->20 min 378->340, m 20->21 min 340->252; n=8: m 12->13 min 368->362, m 15->16 min 402->392, m 16->17 min 392->382, m 20->21 min 518->480, m 21->22 min 480->392, m 25->26 min 608->592, m 26->27 min 592->522, m 27->28 min 522->392 -- so the truncation above is a scope limit, not a proof; the residual windows just named are the exact gap between this program and the paper's claim.
> SCOPE: all eleven rows of the paper's witness table are re-decoded, re-encoded byte for byte and re-evaluated here, and the structural descriptions of the tree witnesses at n = 7..11 (spiders, and two adjacent degree-3 centres with legs) are checked by isomorphism.
> SCOPE: the paper's shard replications (mod 64, 256, 240, 512) are not reproduced; two distance kernels and, for trees, a third edge-split kernel are.  The per-edge-count shard SIZES at n = 10 (the OEIS A054924 comparison) are not recomputed: only the total class count and the set of realizable edge counts are.
