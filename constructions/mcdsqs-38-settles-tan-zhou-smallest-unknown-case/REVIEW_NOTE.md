# A Steiner Quadruple System of Order 38 with Minimum-Colorable Derived Designs

`mcdsqs-38-settles-tan-zhou-smallest-unknown-case`

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
passes. The recorded run reports **32 checks, all passing**:

    VERDICT: ALL 32 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    59734838646641e7b7742407ec4d3c78eb242281460d243add3ed8517801b936

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE: the search that found this system and these two colourings is not re-run; the exhibited object and both colourings are verified exactly, and the lower bound 19 is an exact pigeonhole valid for every STS(37), so the theorem is fully established. The identity q'_0(3,4,38) = 1 + max_x chi'(D_x) is quoted from the cited literature and is applied, not proved, here.
