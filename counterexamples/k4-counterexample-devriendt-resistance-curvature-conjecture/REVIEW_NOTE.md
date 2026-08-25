# A four-vertex counterexample to Devriendt's Conjecture 7.6

`k4-counterexample-devriendt-resistance-curvature-conjecture`

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

    f2e03baaab73f6b991b6c94776462c90ad6b12b75ce12262845af76edd7637b1

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the cited article's own theorems are not reproved here; only the implication p>=0 ==> submodularity is spot-checked, on the 4031 primitive integer weightings of K_4 with conductances at most 4. Searches over conductances are restricted to the integer ranges named in each check, and stability is tested on a finite grid of perturbations rather than an open neighbourhood. The paper's two remarks are proved there analytically and are only spot-checked here: minimality of order four on the 258 integer weightings of K_2, P_3 and K_3 with conductances at most 6, together with 1029 triangle instances against the series-parallel closed form; and persistence under perturbation only at the 729 corners of the box w +/- 1/64, not on any open neighbourhood.
