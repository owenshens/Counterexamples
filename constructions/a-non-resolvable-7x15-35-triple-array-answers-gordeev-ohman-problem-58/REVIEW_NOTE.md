# A Non-Resolvable (7 x 15, 35)-Triple Array

`a-non-resolvable-7x15-35-triple-array-answers-gordeev-ohman-problem-58`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |
| `REVIEW_NOTE.md` | this file |

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package, no external data file
and no network. The program prints one line per check and a closing verdict, and exits 0
only if every check passes. The recorded run reports **29 checks, all passing**:

    VERDICT: ALL 29 CHECKS PASS

Its only input is the 7 x 15 array exhibited in Section 2 of the paper, transcribed into the
program; every quantity the paper states about that array is re-derived, in exact integer
and `Fraction` arithmetic, with no floating point anywhere. What it re-derives: the five
clauses of Definition 1 (binary in rows and columns, e = rc/v = 3, lambda_rc = 3 on all 105
row-column pairs, lambda_rr = 5 on all 21 row pairs, lambda_cc = 1 on all 105 column pairs);
non-triviality and non-extremality; that the cell admits resolutions at all (k = c/e = 5,
rk = v); the 35 row-sets, their 31 distinct values, the fibre-size multiset {1^27, 2^4} and
the four named coincidences; the failure of Definition 2, both from that multiset and,
independently, by brute force over all C(35,5) = 324,632 five-element sets of symbols, none
of which has a common row-set; the linear algebra of the resolvable row-set pattern
(N N^T = 2I + J and the unique solution of N m = 3.1); and that the column design is a
2-(15,3,1) design with exactly 3 parallel classes, no resolution, and automorphism group of
order 6, hence not PG(3,2). It also checks that the array is not a quad array.

Two controls are included so that the resolvability decider is not a program that can only
say "no": on a synthetic row-set map made of five copies of each Fano line the same decider
answers YES, and on a synthetic map with legal class sizes whose seven row-sets are not a
2-(7,3,1) design it answers NO *for that reason*.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by
an exit status, both written by the run harness. The header records the SHA-256 of the
program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    4831224b97f2cc4cc559a3994e0ccded24b5fe666fa1701a2c76dedeaa85296b

## Scope

The program's own closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the CP-SAT search that PRODUCED this array is not re-run and is not needed --
> the array above is the input, and every claim of the paper is re-derived from it by the
> checks above.

> NOT RE-RUN: no minimality, uniqueness or counting claim. This program does not enumerate
> the (7 x 15, 35)-triple arrays, does not decide how many non-resolvable ones exist, and
> does not identify the column design inside any published catalogue of the 80 STS(15)
> isomorphism classes -- it only computes that design's own invariants.

> NOT RE-RUN: nothing about the PREVIOUSLY known (7 x 15, 35)-triple array. Its array is not
> printed here, so the statements that its column design is PG(3,2) and that it is resolvable
> are quoted from the literature, not verified.

> NOT RE-RUN: no bibliographic claim. Whether this is the first non-resolvable triple array
> on any non-extremal parameter set is a literature question and is outside this program.

> NOT RE-RUN: Problems 56 and 57 and Questions 59 and 60 of the same section, and the
> resolvable-array enumeration of the source paper.

The paper claims no priority and no census. In particular it does not assert that no
non-resolvable array on these parameters was previously recorded: the full text of
J. L. Yucas, *The structure of a 7 x 15 triple array*, Congr. Numer. **154** (2002) 43-47 --
a whole paper on this one parameter set -- was not accessible to us, Congressus Numerantium
not being online, so no such statement is made. Nor is anything claimed about the number of
STS(15) classes that carry a non-resolvable array.

The verification is in any case readable by hand. The deciding fact -- that no five of the
35 symbols share a row-set -- can be read straight off the row-set list printed in the
paper, so a referee does not need to run anything.
