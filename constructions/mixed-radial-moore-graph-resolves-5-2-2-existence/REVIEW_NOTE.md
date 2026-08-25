# A Mixed Radial Moore Graph with Parameters (5,2,2)

`mixed-radial-moore-graph-resolves-5-2-2-existence`

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
passes. The recorded run reports **18 checks, all passing**:

    VERDICT: ALL 18 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    bd4e40705697c43e0f176128e8e067805cdf64384eeea10000ce1a80613accc7

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: that the (5,2,2) case is listed as open in the two cited tables is bibliographic and is not recomputed here; the previously settled families (z=1 any r, r=1 any z) and the six pairs (2,2),(2,3),(3,2),(3,3),(4,2),(2,4) are taken as documented input to the least-Moore-bound comparison; no exhaustive search over all mixed graphs of order 52 is attempted, and no minimality of N_1 is tested (the paper makes no such claim). The falsification self-test is seven named local corruptions, five structural and two degree-preserving; it is not an exhaustive perturbation search, so it shows each condition family CAN fail, not that no corruption whatsoever escapes. The second computation of the distances (the bitmask closures) agrees with the breadth-first matrix on every one of the 2704 ordered pairs, but it reads the SAME decoded adjacency, so it cross-checks the traversal and not the transcription: an error in decoding the paper's table would corrupt both arms identically and is caught only by the printed first and last rows and the well-formedness, degree and count checks.
