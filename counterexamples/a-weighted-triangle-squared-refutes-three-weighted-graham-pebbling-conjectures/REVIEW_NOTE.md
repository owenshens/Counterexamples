# A Nine-Vertex Counterexample to the Weighted Graham Pebbling Conjectures as Stated in the E-print of Herscovici, Hester and Hurlbert

`a-weighted-triangle-squared-refutes-three-weighted-graham-pebbling-conjectures`

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
It runs in about ten seconds. The program prints one line per check and a closing verdict, and
exits 0 only if every check passes. The recorded run reports **60 checks, all passing**:

    VERDICT: ALL 60 CHECKS PASS

Its inputs are the objects exhibited in the paper — the three edge weights of the triangle, the
distribution `D`, the distance vector and the graph6 string of the product's skeleton — and it
re-derives every quantity the paper states from them, in exact integer and rational arithmetic
with no floating-point decision. In particular it checks the unsolvability of `D` twice by
structurally different routes: once by enumerating all 229 distributions reachable from `D` with
no use of distances, boxes or potentials, and once by mechanising the paper's own potential
argument (every move of a successful play costs at most 1/9, only 16 of the 36 directed moves
do, and the 112 states reachable from `D` by those 16 never place a pebble on the root). It also
reproduces five pebbling numbers published by other authors, Chung's Theorem 3 on a 144-instance
grid of weighted `K_2` products, and the weighted `K_4` cover-pebbling example of the source
paper itself; and it runs anti-controls on which it must, and does, report no violation.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    9be96339491d8a6145963995743d53c8c1eded0fa8b6cbe1fd353b135139eca5

## Scope

The decisive inequality of the paper, `pi(G x G,(1,1)) >= 10 > 9`, is proved by hand in two
cases and needs no program; the program re-derives that proof but the paper does not rest on it.
The **equality** `pi(G x G,(1,1)) = 10` is a machine result and is labelled as one in the paper,
since equality additionally requires every distribution of ten pebbles to be solvable, which no
step of the hand proof performs.

The statements refuted are the ones quoted from the e-print source; the printed text of
*Discrete Math.* **312** (2012) 2286–2293 was not accessible, so the paper's title, abstract and
corollary claim the refutation for the e-print formulations only.

The program's closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the search that FOUND this witness -- a census of 47,336 tight weighted instances
> over six cells -- is not repeated here, and no minimality claim of any kind is made or checked.
> The smallest counterexample, over weighted graphs or even over weighted triangles, is not
> determined by this program.

> NOT RE-RUN: the printed text of Discrete Math. 312 (2012) 2286-2293 is paywalled and was never
> read, so the printed CONJECTURE NUMBERS are unverified; the paper quotes the statements from
> the e-print source by LaTeX label and line range instead, and this program checks mathematics
> only, never bibliography.

> NOT RE-RUN: nothing here bears on Graham's (unweighted) conjecture, on the unweighted rooted
> conjecture, or on Chung's proved region beyond the 144-instance grid check above.

The program also checks two further witnesses, `K_3(2,3,2) x K_3(2,4,3)` and `K_3(2,4,3)^2`,
which the paper does not discuss: their rooted values are `13 > 12` and `18 > 16`, so each
refutes the rooted conjecture, but the hand proof of the paper's Section 2 does not cover them
and they do **not** refute the unrooted statement, since `pi(K_3(2,4,3)) = 5` raises the
unrooted right-hand sides to `3 * 5` and `5 * 5` while their unrooted pebbling numbers are 14
and 19. The unrooted claim, and hence the answer to Sieben's question, attaches to the witness
of the paper alone.
