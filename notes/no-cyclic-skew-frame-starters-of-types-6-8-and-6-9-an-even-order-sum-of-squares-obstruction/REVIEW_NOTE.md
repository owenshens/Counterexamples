# No Cyclic Skew Frame Starters of Types 6^8 and 6^9:\\ An Even-Order Sum-of-Squares Obstruction

`no-cyclic-skew-frame-starters-of-types-6-8-and-6-9-an-even-order-sum-of-squares-obstruction`

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
passes. The recorded run reports **70 checks, all passing**:

    VERDICT: ALL 70 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    91134b47271d98deb0331b6a0c58310a57e6d5a5b1f17f220d0b0db01829ea36

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the exhaustive census of the 4^11 exact-cover instance.  The paper
>   does not claim the 4^11 cell and this program does not decide it.  Only the MODEL
>   arithmetic above (80 columns, 640 rows, 869 sub-shards) is checked; no search is
>   performed here, and the engine that ran that census is not part of this folder.
> NOT RE-RUN: nothing about NONCYCLIC groups of order 44, 48 or 54, about skew Room
>   FRAMES as opposed to starters, or about STRONG (non-skew) frame starters.  Those
>   boundaries are asserted in the paper on the cited literature, not recomputed here.
> NOT RE-RUN: the transcription itself.  The 35 verdicts, the nine "exhaustive search"
>   authorities and the three printed starters are read off arXiv:2211.12367v1 by hand;
>   this program checks that they are internally consistent with the criterion, which
>   is a strong test of the transcription but not a substitute for reading the source.
