# The Coprime Buratti–Horak–Rosa Conjecture Holds for the Support \1,6,18\

`hamiltonian-paths-settle-coprime-bhr-support-1618`

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
passes. The recorded run reports **26 checks, all passing**:

    VERDICT: ALL 26 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    ce10e45eecea4b84b226105c583f7bf9623c484c9df1fce7e2f8e909702363dd

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: (i) the paper's node counts (1.51e7 for the largest target, 4.11e9 in total) and its 3792/43/288 split of targets by reduction: the cap schedule in force here is [500, 5000, 100000] nodes per starting vertex, whose last rung is exactly the paper's own budget, so on each reduction this search is allotted neither more nor less than the paper allots itself, the cheaper rungs being tried first, and a target is credited to whichever reduction settles it first at the rung that settles it, so both the node counts and the split printed above are this program's and not the paper's; the SPLIT NOT COMPARABLE line above says which of those obstacles applies to this run, and running with --paper-protocol and no --stride= removes them all; (ii) the cited prior work, which covers the orders outside {47,49,59} and which is what makes these three the only open ones; (iii) the claim that (13,6,39) is one of ten triples that prior work leaves open at v=59, a fact about that reference.
