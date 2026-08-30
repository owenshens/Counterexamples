# A 72-Vertex Counterexample to a Spectral-Gap Conjecture of Greaves and Zhu

`a-72-vertex-pancake-counterexample-to-the-greaves-zhu-spectral-gap-equality`

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
passes. The recorded run reports **109 checks, all passing**:

    VERDICT: ALL 109 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    eb7512e187240f9652d5e34524311b2891815f75ce79bbd8cb9704f59ca050f9

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the exact value psi_{6,2} = 1.  This program proves only psi_{6,2} <= 1,
> NOT RE-RUN:   which is all the refutation needs; the matching lower bound is Blanco's
> NOT RE-RUN:   published theorem (arXiv:2608.15398v1) and is quoted, not verified here.
> NOT RE-RUN: lambda_2(A(P_6(2))) and lambda_2(A(P_6(3))) are never computed -- only
> NOT RE-RUN:   bounded below, by the exhibited eigenvector and by lambda_max(T).
> NOT RE-RUN: the cells m = 3, 4, 5 at n = 2.  Section E shows only that the separator
> NOT RE-RUN:   test does not fire there.  The test bounds psi_{m,2} above and beta_{m,2}
> NOT RE-RUN:   below, so a negative is NOT evidence that psi_{m,2} = beta_{m,2}; those
> NOT RE-RUN:   cells are left undecided in both directions.
> NOT RE-RUN: MINIMALITY.  Nothing here shows (6,2) is the smallest counterexample. The
> NOT RE-RUN:   cell (2,3), on 48 vertices, is not decided by this program at all, and
> NOT RE-RUN:   m = 2 uses a different connection set (q_2 = 1) that is not built here.
> NOT RE-RUN: cells with n >= 4, and m > 20 at n = 2, are not examined; the family claim
> NOT RE-RUN:   for all m >= 6 rests on section C's identities, not on section D's table.
> NOT RE-RUN: the asymptotic expansion beta_{m,2} = 3 eps - eps^2/4 + O(eps^3).
