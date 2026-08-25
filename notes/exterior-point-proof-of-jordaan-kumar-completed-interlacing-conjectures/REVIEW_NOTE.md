# An Exterior-Point Proof of Two Completed-Interlacing Conjectures

`exterior-point-proof-of-jordaan-kumar-completed-interlacing-conjectures`

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
passes. The recorded run reports **12 checks, all passing**:

    VERDICT: ALL 12 CHECKS PASS.

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    1c7388f2105a60160a1ddd66d7f0d81f654cf5f49c1fcab6391bc2cf1bd7e468

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE: Lemma 2.1 is tested on random RIGID RATIONAL instances (integer zeros for G, rational interlacing zeros for Q, c forced by the degree cancellation).  It is NOT proved here for general real data, and the actual Meixner-Pollaczek / pseudo-Jacobi instances have irrational zeros; those are covered separately, and only pointwise, by checks 11 and 12.
