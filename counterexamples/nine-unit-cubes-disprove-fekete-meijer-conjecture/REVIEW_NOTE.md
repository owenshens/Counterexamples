# A Unit-Cube Box-Visibility Representation of K_9

`nine-unit-cubes-disprove-fekete-meijer-conjecture`

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
passes. The recorded run reports **24 checks, all passing**:

    VERDICT: ALL 24 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    4b405c18a5df23b10a5ed43db068db2677bc84e5f976c5189dc4a00c0d31aba4

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the cited nonexistence result for K_10 and the cited K_8 construction are external and are not recomputed here; only the K_9 representation is verified.
