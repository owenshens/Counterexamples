# Counterexamples to the Hypercube Maximum Propagation-Time Conjecture

`hypercube-zero-forcing-sets-disprove-propagation-time-conjecture`

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
passes. The recorded run reports **22 checks, all passing**:

    VERDICT: ALL 22 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    35f5d85df2f2b66ee35671a94088a2897a5b098044038798db13c4860b15cd0e

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> info: NOT RE-RUN -- the lower bound Z(Q_d) >= 2^(d-1) for d = 5,6,7 is taken from the literature, not recomputed here; only the upper bound is reproved above, and only d <= 4 is settled exhaustively.
> info: NOT RE-RUN -- no complete census of the minimum zero forcing sets of Q_5, Q_6 or Q_7 is attempted (C(32,16), C(64,32) and C(128,64) candidates), so the exact values of PT(Q_6) and PT(Q_7) are not determined; only the lower bounds 18 and 43 are established, which is exactly what the paper claims.
> info: NOT RE-RUN -- the separate interval assertion of the original conjecture is outside the paper's claim and is not examined.
