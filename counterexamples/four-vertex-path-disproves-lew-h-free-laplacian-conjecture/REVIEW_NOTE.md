# Counterexamples to Lew's H-Free Laplacian Conjecture

`four-vertex-path-disproves-lew-h-free-laplacian-conjecture`

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
passes. The recorded run reports **38 checks, all passing**:

    VERDICT: ALL 38 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    0a3c429c34286bf0cd2699b8314167bd979711f13c08649bd8182c7dfa535325

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: no exhaustive census beyond the ranges printed above -- the k=1 census covers all labelled graphs on at most 6 vertices, the order censuses all graphs on at most 5 vertices (all triangle-free graphs on exactly 4 vertices at k=2), Turan numbers are enumerated for n <= 6 and the increment formula is checked for r <= 7, k <= 40; the general theorems for arbitrary H and k are verified only on the listed instances. Remark 6 is only partly evaluated: its consistency assertion eps_k <= ex(k;H)+(4k-2)sqrt(k) is tested on the six instances listed above and not for all r and k, and its asymptotic comparison (excess O(k) against ex(k+1;K_r) = Theta(k^2)) is verified only through the per-instance excess bound, not in general. NOT CHECKED AGAINST ANY SOURCE: the wording and hypotheses of Conjecture 7.1, and the approximate bound, are transcribed from the paper's own quotation of arXiv:2601.17575v1; this program has no access to that preprint, so a transcription error in the conjecture -- for instance a dropped hypothesis on H, which Remark 8 shows would matter -- would not be detected here.
