# A Proof of the Grinberg–Lafrenière Lower Bound for the Random-to-Below Shuffle

`a-proof-of-the-grinberg-lafreniere-lower-bound-for-the-random-to-below-shuffle`

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
passes. The recorded run reports **72 checks, all passing**:

    VERDICT: ALL 72 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    9943242ec03abba0658c109c30dd72b185c3ee70eb6e3f42660e6b3f167d9c05

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the exact rational S(n) is formed for n <= 80 only; for 81 <= n <= 400 the quantities S(n)/n, D_n, L_n and M(n) are computed by the outward-rounded interval arithmetic of Step 1 rather than as exact Fractions (the enclosures are rigorous, so no inequality above is weakened, but the huge-numerator exact form is not built).
> NOT RE-RUN: no statement about n > 400 is verified NUMERICALLY here. The tail for n >= 36 is settled by the monotonicity of (IV) and (V), which Steps 6 and 7 verify only on n = 36..1999; beyond that the note relies on the term-by-term monotonicity proofs of Sections 3 and 4, which this program does not symbolically differentiate.
> NOT RE-RUN: inf_{n>=2} M(n) is NOT determined. M(7) = 0.2674107869 is the minimum over the cell 2..35 only; Step 8 verifies M(n) > 0.0034553 on 36..400 but no exhaustive minimisation is performed there, so the note does not claim -- and this program does not check -- that 0.2674107869 is the optimal constant.
> NOT RE-RUN: D_n -> 1 - gamma is NOT verified. Step 5 bounds D_n only from ABOVE, so the corollary that the authors published upper bound is asymptotically sharp to within (log 2 - gamma) n is outside this program.
> NOT RE-RUN: the probabilistic content. That S(n) IS the expected value of the authors strong stationary time is their Theorem, quoted in Section 1 and not reproved here; this program verifies the inequality about the sum.
