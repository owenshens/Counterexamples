# An order-six counterexample to the Brualdi–Liu even-permanent conjecture

`order-six-counterexample-brualdi-liu-even-permanent`

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
passes. The recorded run reports **29 checks, all passing**:

    VERDICT: ALL 29 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    4ddc78b506b394dd3b2457b306a8add002c13d4e9baa44f5c45b5036e22c9a0c

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: no exhaustive census over all 4-element subsets of the 360 even permutation matrices of order 6 (nor of the 60 of order 5) was attempted; the order-6 sweep is over the weights on the four exhibited vertices only, and the order-5 refutation sweeps the weights on one exhibited 4-vertex set. The order-4 refutation is exhaustive over 4-vertex sets at denominator 10.
