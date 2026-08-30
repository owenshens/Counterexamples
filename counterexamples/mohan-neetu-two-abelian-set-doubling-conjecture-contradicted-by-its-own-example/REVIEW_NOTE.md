# A Four-Element Counterexample to Conjecture 6.1 of Mohan and Neetu on Two Disjoint Abelian Sets

`mohan-neetu-two-abelian-set-doubling-conjecture-contradicted-by-its-own-example`

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
passes. The recorded run reports **244 checks, all passing**:

    VERDICT: ALL 244 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    5c362f18ac6ac3fc6e63f119d5c9e01108def048c2ba7f1e9e223f333ad860e1

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: torsion-freeness of K is confirmed only on the 81-element coordinate box r,n in [-4,4] and only for exponents 1..200; the unbounded statement is the hand proof of Lemma 1.
> NOT RE-RUN: the family of equation (2) is machine-checked for 1 <= m <= 40 only; the statement for all m is Theorem 3, proved by hand.
> NOT RE-RUN: no search for a witness of smaller cardinality (in particular nothing is decided about |S| = 3), for a smaller doubling, or in any host group other than K = BS(1,-1).
> NOT RE-RUN: this program does not fetch the e-print, so the byte locator, the line numbers and the SHA-256 of section 1 of the paper are NOT re-checked here; they were read off the source file directly.
> NOT RE-RUN: nothing here bears on Conjecture 6.2 or Conjecture 6.3 of [MN], nor on whether any repaired form of Conjecture 6.1 is true.
> NOT RE-RUN: the provenance question of whether the line-724 object is inherited from the closed-access companion paper is not a computation and is not addressed.
