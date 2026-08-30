# A Multilinear Polynomial with Constants Whose Image on M_2(K) Is Not a Vector Space

`an-interior-rank-one-coefficient-refutes-the-generalized-kaplansky-lvov-conjecture`

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
passes. The recorded run reports **89 checks, all passing**:

    VERDICT: ALL 89 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    b0a1f7a2b994cee2d08fe2118a6e4e99af59c0131adff22562561200b790e8d0

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE -- WHAT THIS PROGRAM DOES NOT COVER
> NOT RE-RUN: the PUBLISHED wording of the conjecture.  Every locator is the e-print arXiv:2312.13865v1 (conc.tex, lines 21-23).  The Wiley text of Mathematika 71 (2025) e70031 returned HTTP 403 and was never read, so whether the published sentence restricts the coefficient positions is settled neither here nor in the paper.
> NOT RE-RUN: the generalized-commutator cell of the original row -- omega(x,y) = A x y - B y x over an algebraically closed field of characteristic 0.  Parts 5(d) and 5(e) show only why a rank obstruction cannot reach that shape; nothing here decides it, and the paper claims nothing about it.
> NOT RE-RUN: any infinite, characteristic-0 or algebraically closed field.  Every enumeration runs over GF(2) or GF(3).  The refutation does not depend on them: Part 1 is an exact computation valid over an arbitrary field, and Parts 2-4 verify the structural identities that carry it to every field and every n >= 2.
> NOT RE-RUN: n >= 5; at n = 4 only q = 2, and only r in {1, 2}.  The rank-r family claim is PROVED in the paper by rank factorization -- Part 4 is a finite confirmation, not a proof by exhaustion.
> NOT RE-RUN: the two AWS finite-field censuses filed with the original row (census_kl_job1.py, census_kl_job2.py, numpy-based, 78 s and 1982 s).  Their stdout was never captured to a file and the claim they supported was struck in review as off-cell.  Their decisive negative-polarity integers -- |image| = 10, 50, 33, 339 -- are recomputed from scratch in Part 3 here, with no numpy and no floats.
