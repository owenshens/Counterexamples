# A Counterexample to the Class 69 Equidistribution Conjecture on Involutions

`size-three-involutions-disprove-fang-et-al-class-69-equidistribution`

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
passes. The recorded run reports **21 checks, all passing**:

    VERDICT: ALL 21 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    aa6a3b7c54507f52b71bff155e1110295ca078e02d1046a18e6f627877074abf

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the equidistribution over the full symmetric groups is a result quoted from the literature for all n; it is confirmed here only for n <= 9, and the involution split only for n <= 11. That the four shadings above are the ones the conjecture displays, and that Class 69 is the last open length-two class, are bibliographic facts about the source and are not machine-checkable.
