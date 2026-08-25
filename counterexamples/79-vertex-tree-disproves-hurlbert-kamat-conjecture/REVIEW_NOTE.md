# A 79-Vertex Counterexample to Conjecture 4.1 of Hurlbert and Kamat

`79-vertex-tree-disproves-hurlbert-kamat-conjecture`

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
passes. The recorded run reports **25 checks, all passing**:

    VERDICT: ALL 25 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    d9c430af1f862dfbea910f5a55844b98888529ee1af28b18e3ee48e12cfe608c

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> scope: the exhibited tree, all its hypotheses, all four asserted counts and all three polynomial identities are recomputed from the edge list alone by rooted dynamic programming; each of the four asserted counts is then recomputed a second time by the paper's algebraic route -- formula (5) for [z^m]S^q, the closed form for [z^m]H_t, and coefficient convolution, binomial arithmetic only, touching neither the graph nor the polynomial primitives the DP uses -- and the two routes agree digit for digit.
> NOT RE-RUN: the statement being refuted. This program reads no literature, so it cannot confirm that the conjecture quoted in the paper -- every tree with no vertex of degree two is r-HK for every r with 1 <= r <= alpha -- is the source statement, and in particular it cannot rule out a restriction on r in the source. That gap is load-bearing here: mu(T)=26, so the exhibited failure rank r=18 exceeds mu(T)/2=13, and the ranks at which T fails to be r-HK are exactly [18, 19, 20, 21, 22, 23]. Under any reading that admits only r <= mu(T)/2, the threshold recurring in the Holroyd-Talbot literature, every count above remains correct and none of them contradicts the conjecture. Also not tested: minimality of the order 79 over all trees, which the paper does not claim, and the census sweep, which covers only branch sizes in {0,3,...,11}.
