# The Clow–Zaguia Cop Bound is Tight at Independence Number Four: f(4)=3

`the-clow-zaguia-cop-bound-is-tight-at-independence-number-four-so-f-4-equals-3`

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
passes. The recorded run reports **37 checks, all passing**:

    VERDICT: ALL 37 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    6291dc37a520a7075e9e7d37a1a33504e305b0d14fa241e957f0829bb9440758

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: no 3-cop solver was run. c(G) <= 3 is established from the printed dominating set (c <= gamma) and, independently, from the target paper's own theorem; the exact solver above decides only the 2-cop game.
> NOT RE-RUN: perfection of K_m x K_m for m >= 5 was NOT brute-forced -- the strong-PGT scan is exponential in the order, so only the 9-, 12- and 16-vertex graphs were scanned. For m >= 5 perfection is by theorem (line graph of a bipartite graph, plus the weak PGT), not by test.
> NOT RE-RUN: no census of connected perfect graphs was performed, so the MINIMUM ORDER of a connected perfect graph with alpha = 4 and cop number 3 is not determined here; it is only bracketed 10 <= min <= 16, the lower bound being a published value and the upper bound this witness.
> NOT RE-RUN: nothing here computes f(k) for any k >= 5. The family gives 3 <= f(k) <= k-1 for k >= 4 and this program checks only the lower half at k = 5 and k = 6. Problem 5.2 itself is NOT solved.
> NOT RE-RUN: the family check above stops at m = 6; m = 7 and beyond rest on the case analysis in the paper, which uses only m >= 4, and were not replayed.
> NOT RE-RUN: no literature was fetched. Whether the value c(K_n x K_n) = 3 is already proved in print -- in particular in Neufeld and Nowakowski, Discrete Math. 186 (1998) 253-268, which we could not read -- is a bibliographic question no program can settle.
