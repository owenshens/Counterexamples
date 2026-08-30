# Polynomial Certificates for the First Two Partial Quotients of ( (3)-S_N)/2, and the Third for N 2

`polynomial-certificates-pin-the-first-three-partial-quotients-of-the-zeta-three-tail`

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
passes. The recorded run reports **95 checks, all passing**:

    VERDICT: ALL 95 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    cfcf324171ef9310a73f056ff83b675648ac5a0d15cb0153514d376ebf9ea5d3

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> scope: what this program does NOT cover
> NOT RE-RUN: the trailing ellipsis of the source's display. Nothing here
>   addresses a_4, a_5 or their quasi-periods, so the source's Conjecture
>   prop:gcf as a whole is NOT verified and NOT claimed; only the three
>   quotients its equation eq:gcf displays are settled.
> NOT RE-RUN: any Lean formalisation. No Lean was written and none was
>   inspected; the source's development is not public, so its single
>   `sorry` could not be diffed against any artifact. What is discharged is
>   the mathematical statement the source describes, not a machine proof.
> NOT RE-RUN: the search that found the antidifferences psi and psi3.
>   Nothing depends on it -- a positive-coefficient certificate is
>   self-checking once the polynomial is exhibited.
> NOT COVERED: irrationality, transcendence or normality of zeta(3), and
>   the growth of the partial quotients. A decimal of zeta(3) enters this
>   program in exactly one place, the bracket control above, and no proved
>   statement rests on it.
> NOT COVERED: the source's Theorem thm:main. It concerns sqrt(R_{N_k}),
>   a different real number, and is used here only as a forced-positive
>   control on the continued-fraction routine; it is not reproved.
> NOT COVERED: minimality of the certificate. No claim is made that the
>   five-term enclosure is the shortest that works, only that the
>   four-term one measurably does not.
