# Cubic Graphs of Girth 11 on 138 and 142 Vertices

`cubic-girth-eleven-graphs-orders-138-142`

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
passes. The recorded run reports **51 checks, all passing**:

    VERDICT: ALL 51 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    ea140b634871765fc414f02a851edc299bb37ef19c8be1fa7341812ed7895476

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: the exhaustive searches of the cited works are not repeated. Specifically, that no connected cubic girth-11 graph exists on 128, 130, 132, 134 or 136 vertices is neither claimed nor tested; that every even order from 144 upward is realisable is taken as an input; and the SHA-256 of the source archive stated in the paper is not recomputed, since this program reads no external file and instead carries the four graph6 strings (the three order-140 lines and the order-152 line) as literals.
