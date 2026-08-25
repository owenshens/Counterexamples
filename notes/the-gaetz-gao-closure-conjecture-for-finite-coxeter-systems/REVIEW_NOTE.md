# The Gaetz–Gao Closure Conjecture for Finite Coxeter Systems

`the-gaetz-gao-closure-conjecture-for-finite-coxeter-systems`

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
passes. The recorded run reports **53 checks, all passing**:

    VERDICT: ALL 53 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    2f2f768706b0fa9530522f1bd9707110ff1d15009f118d5c0ba9a50dd9bc0593

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE: the paper's H_4 census (Proposition 5) is re-derived in full (14400 elements, 230400 BP tests, 201991 pairs, 0 failures).
> NOTE Theorem 1 is additionally verified outright for 25 further finite Coxeter systems, among them I_2(5)xH_4 (144000 elements) and A_2xH_4 (86400), plus I_2(m) for every m <= 200 (which includes G_2 = I_2(6)).
> NOTE NOT re-derived here: Theorem 1 as a whole.  Its general type-A/D/E cases and its components B_n, F_4, E_7, E_8 rest on the Gaetz-Gao theorem for finite Weyl groups, which the paper cites and this program does not reprove; only A_1..A_7, D_4..D_6, E_6 and products of these are re-verified.
> NOTE NOT re-derived here: I_2(m) for m > 200, the types B/F needing 2cos(pi/4) (outside Z[phi]), E_7 and E_8.
