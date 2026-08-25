# Mesh-pattern pairs 113–116 are not equidistributed with the X^(1)/Y^(1) family: a counterexample to an open problem of Lv and Zhang

`mesh-pair-113-refutes-lv-zhang-cross-pair-equidistribution`

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
passes. The recorded run reports **30 checks, all passing**:

    VERDICT: ALL 30 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    8932406ee1a809f15d99afc8e024085f2f3e7e5f4bb1e64f07e735ba443ac6a1

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the enumeration above covers S_n for every n <= 8, for the two pairs the paper exhibits in full -- the pair labelled 113 and the pair X_9^(1) -- and, over S_n for n <= 6, for the inverse of the 113 pair, which the paper's own argument identifies with a second member of the four-pair family. The two other members of that four-pair family, labelled 115 and 116, and the remaining nineteen members of the X^(1)/Y^(1) family, are named only by label in the source paper, whose shading tables are not reproduced here, so their distributions were not recomputed: two of the four P-pairs are not enumerated here at all and rest entirely on the citation to Lv and Zhang's Theorems 6.1 and 6.2, and for the nineteen unenumerated Q-pairs the refutation rests on the closed form quoted for that family, which is checked here cell by cell against the directly enumerated distribution of X_9^(1).
