# Two Residual Cases of Simmons's Polygonal Bigraph Conjecture

`residual-polygonal-bigraphs-close-simmons-finite-range`

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
passes. The recorded run reports **43 checks, all passing**:

    VERDICT: ALL 43 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    8827318baeb3c4eb0e072fc227810d49cd1d6c1efcf0d205467b96b98df20021

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: what the originally printed residual table did or did not display cannot be re-derived from arithmetic; the checks above establish which parameters are residual and which of them reduce, not which lines a printed table contained.
> NOT RE-RUN: no claim is made or tested here for m > 34; the statement verified is the finite one.
