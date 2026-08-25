# The Hessian Arrangement Requires Four Added Lines for Supersolvability

`hessian-arrangement-needs-four-added-lines-disproves-dimca-kuhne-pokora-question`

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
passes. The recorded run reports **41 checks, all passing**:

    VERDICT: ALL 41 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    6bc3d6d6b28bab1f01a12438f7ccdfbd266d9ee7626e57cdc16e280ac44e9c4f

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: this program re-derives every number it prints from the twelve printed linear forms, the four printed added lines and the point (1:-1:0), but it does NOT reprove the external inputs it applies -- the rank-three equivalence between having a modular point and lattice supersolvability (G2), the agreement of Kabat's resolution definition of extSS with min |B \ H| (G3), the du Plessis-Wall freeness criterion (G4), the statement of Dimca Theorem 1.12(2) (G5), and the definitional fact that at most one line of H passes through a point off Sing(H), on which the off-Sing(H) bound rests (G1); the resolution of the labels Theorem 1.12(2), Example 6.5(i), Remark 2.3 and Definition 1.3 to numbered items INSIDE arXiv:2503.01624v7, arXiv:2505.01733v4 and arXiv:2201.04856v1 was done by hand against those preprints, since they are not shipped here, and step 10 verifies only that the paper cites those four labels and that its own bibliography sends each to the arXiv version used here (G6); the three figures |Sing(B)| = 45 and the forced-added-line counts 6 for a modular point in V and 6 for one off Sing(H) are derived here and are NOT printed by the paper, whose lower-bound argument claims only k >= 4 in each case; the six falsification probes above are computed in this run, but no separate mutation transcript accompanies the paper; and any transcript shipped beside the paper that reports 40 checks, or a step 10 gap reading 'no paper .tex was found', was produced by an earlier build of this file and is NOT a run of this one -- this version embeds its bibliographic input and therefore records step 10 in every run, whether or not a paper source sits beside it, for 41 checks in total.
