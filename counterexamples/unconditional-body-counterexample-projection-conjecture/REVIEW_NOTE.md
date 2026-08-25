# A Counterexample to Conjecture 3 of Fradelizi–Manui–Meyer–Ndiaye

`unconditional-body-counterexample-projection-conjecture`

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
passes. The recorded run reports **49 checks, all passing**:

    VERDICT: ALL 49 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    884fe18b77cf755c673b600ff41d90fbeeb367b63357671e3aa70ae3d11c80da

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE (1): the positive cases n = 1 and n = m-1 of the classification are quoted from the cited work and are not recomputed.
> NOT RE-RUN HERE (2): the case n = m is supported only through the orthant reduction above, not by a 4-dimensional volume computation.
> NOT RE-RUN HERE (3): the product identities P_E K = P_{E_0}K_0 x [-1,1]^{n-2} and P_E(K cap R_+^m) = P_{E_0}K_+ x [0,1]^{n-2} are verified as sets only for the single pair (m,n) = (5,3), and there only on a rational grid; the all-(m,n) statement above is the arithmetic consequence 2^{n-2}*(1399/336) > 2^n, enumerated for m <= 20 only.
> NOT RE-RUN HERE (4): checks whose detail line says 'grid' are finite samples, not proofs over R^m; the load-bearing orthant reduction and the K_0 encoding also have exact, sample-free counterparts above.
> NOT RE-RUN HERE (5): the hypothesis class is taken from the paper. Proved above: K_0 is invariant under all 16 coordinate sign changes (hence unconditional) and under no nonidentity coordinate permutation, so it satisfies the sign-invariance hypothesis and not the stronger permutation-invariant one. That sign invariance is the whole hypothesis of the conjecture being refuted is transcribed from the cited work, which is not part of this material and is not read by this program.
> NOT RE-RUN HERE (6): no computation here touches bibliographic data -- the preprint identifier and version, the author names, and the pinpoint citations to the two propositions that supply the positive cases n = 1 and n = m-1. The negative half of the if-and-only-if classification (every m >= 4, 2 <= n <= m-2) is recomputed above; the positive half rests on those citations together with the elementary case n = m.
