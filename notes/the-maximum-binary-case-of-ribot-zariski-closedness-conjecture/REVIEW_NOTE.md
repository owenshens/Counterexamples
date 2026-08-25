# The Maximum Binary Case of Ribot's Zariski-Closedness Conjecture

`the-maximum-binary-case-of-ribot-zariski-closedness-conjecture`

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
passes. The recorded run reports **35 checks, all passing**:

    VERDICT: ALL 35 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    863ff75154ee185efab0e86db0e961778c69c049bcc684f7ff508836135153c9

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the independent-set census for d >= 7 (search tree grows ~10^3 per dimension; d=6 already costs ~9*10^7 nodes).  Every other check below runs for the full stated range.
