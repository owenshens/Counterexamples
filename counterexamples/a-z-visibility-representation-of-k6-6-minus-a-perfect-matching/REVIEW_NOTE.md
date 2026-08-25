# A z-Visibility Representation of K_6,6K6,6 Minus a Perfect Matching

`a-z-visibility-representation-of-k6-6-minus-a-perfect-matching`

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
passes. The recorded run reports **66 checks, all passing**:

    VERDICT: ALL 66 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    22f074564e707cf84eeea3433d8e2568b1d7ee381eb3330f8a1ef264662ba167

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> INFO Not reproduced, because they are not computations: the bibliographic claim that Conjecture 6.5 of Bose et al. (1994) is the statement quoted, and the general proof of the paper's unit-cell Lemma (whose conclusion is nevertheless confirmed on this instance by the independent exact-geometry method in section 4); and, in section 7, the step from the ONE weakening of 'joins' computed there -- cylinder axis on both rectangles, end disks not contained -- to every weaker reading, which is the one-line monotonicity remark printed there and not a search over cylinders whose axis touches neither rectangle.
