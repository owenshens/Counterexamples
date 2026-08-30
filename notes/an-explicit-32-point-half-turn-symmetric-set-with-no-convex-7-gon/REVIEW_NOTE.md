# An explicit 32-point set with a free half-turn and no convex 7-gon

`an-explicit-32-point-half-turn-symmetric-set-with-no-convex-7-gon`

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

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    ce68b660416ca793f0b3119daf6077a3ae54157d7f8b4c26a168aa94b8aa28ce

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> Scope
> NOT RE-RUN: (1) the SEARCH that produced W* is not re-run here and is not reproducible from a seed -- each restart of the original grid descent was cut by a 90-second wall-clock budget, so W* is pinned as the integer list printed in the paper and not as a random seed; nothing above depends on how it was found. (2) The claim that 31 is MINIMAL for the largest absolute coordinate of such a set is NOT made and NOT tested: 31 is where a descending ladder of grid bounds stopped on budget, and no exhaustive search of any smaller grid was performed. (3) NO statement about g(7) is re-derived, because none is made: g(7) is unknown, g(7) >= 33 is classical and independent of this object, and the 4-fold negative result at n = 32 is quoted from the source paper and NOT recomputed here. (4) The source paper's own published 16-point 4-fold configuration is NOT re-counted here (its coordinates are not reproduced in this note); the calibration of the convex-position criterion against that published object was done elsewhere, and the criterion is instead pinned here by the four controls of Step 6. (5) The nine further 32-point witnesses at larger coordinate bounds mentioned in the note are NOT verified here; only W* is. (6) The REFLECTION (mirror) symmetry of order 2 at n = 32 is untouched: this program says nothing about it. (7) In Step 8 the exclusion of every rotation order s >= 7 is checked on one concrete rounded regular octagon, not on all such s; the general argument (an orbit of size s is a regular s-gon, and any 7 of its vertices are in convex position) is a hand proof and is given in the note.
