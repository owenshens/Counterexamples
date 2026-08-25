# An Exhaustive Census for the n=20 Case of the\\ Aouchiche–Caporossi–Hansen Tricyclic-Energy Conjecture

`an-exhaustive-census-for-the-n-20-case-of-the-aouchiche-caporossi-hansen-tricyclic`

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
passes. The recorded run reports **64 checks, all passing**:

    VERDICT: ALL 64 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    5ac349bdfcbbb536f0088e0b6d2b691a2116a4aae1bf54539397311a685a905b

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT REPRODUCED (stated explicitly):
> * The paper's energy scan over all 3288208176 isomorphism classes is NOT re-run here.
> Regenerating those classes and computing about 3.3e9 twenty-by-twenty
> symmetric eigensolves is thousands of CPU-hours; this program instead
> (i) re-derives the census TOTAL 3288208176 exactly, by Burnside/Polya counting
> plus the multiset logarithm, with no use of geng, and (ii) re-runs the
> energy maximisation exhaustively over two slices that do fit the budget:
> the class T(20) of connected 20-vertex 22-edge graphs with all degrees in
> {2,3} (8157 isomorphism classes, containing both P and Q), and the complete
> one-edge-swap neighbourhood of P (3363 connected graphs, most of them outside
> T(20)).  Classes of the census outside those two slices are not rescanned.
> * The paper's shard bookkeeping (997 geng shards, and the confirmation run
> with 512 shards) is not reproducible from a single-process stdlib program,
> and no per-shard record counts are published to compare against.
> * The paper's binary64 figures - the counts 0 and 1 of records within 1e-7 of
> the maximum, and the 5.15e-13 agreement between two eigensolvers on 20000
> sampled graphs - concern that floating-point pipeline itself; this program
> replaces floating point by exact rational enclosures for the graphs it does
> examine, so those two numbers are corroborated in spirit, not re-measured.
> * The paper prints no graph6 string, so none of the three graph6 strings used
> here is a quotation from it.  P's and Q's are pinned to the paper's own
> words by checks A3 and B2 respectively; the third is an alternative
> labelling of P whose claimed nauty/labelg canonicity is NOT verified (nauty
> is not used), so checks B8, B9 and B10 exercise this program's graph6 and
> isomorphism code on a relabelled copy of P and corroborate nothing printed
> in the paper.
> * The conjecture for n >= 22 is outside the paper's scope and is not tested.
