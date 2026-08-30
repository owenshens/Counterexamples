# A nine-vertex P_4-sparse minimal (2,3)-polar obstruction\\ that is not a cograph

`a-nine-vertex-p4-sparse-minimal-2-3-polar-obstruction-that-is-not-a-cograph`

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
passes. The recorded run reports **71 checks, all passing**:

    VERDICT: ALL 71 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    119369d3571511ce67a47773edcfc2ff3f2d9656087d8bf55fb0478c71a99e16

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE
> NOT RE-RUN: the census of ALL P_4-sparse graphs on at most 12 vertices at the cell (2,3). Group E exhausts ONLY the thin spiders with |S| = |K| = 2 whose head is a disjoint union of exactly three cliques (16 shapes). Nothing here shows that 9 is the LEAST order of a P_4-sparse non-cograph minimal (2,3)-polar obstruction, and the paper claims no such thing.
> NOT RE-RUN: any cell other than (2,3) and its complement-dual (3,2). No claim is made or checked for general s, k, for k >= 4, for thick spiders with |K| >= 3, or for the finiteness or classification of the family of such obstructions.
> NOT RE-RUN: the list of the 50 published P_4-sparse minimal 2-polar obstructions is not reproduced here; what is checked is the weaker, self-contained fact that H* is not a minimal (2,2)-polar obstruction at all, which suffices to place it outside that list.
> NOT RE-RUN: nothing in this program is evidence about the literature. Priority, and in particular whether any earlier source exhibits this graph, is not a computation.
