# Vu's Cover-Ideal Linearity Question Through Twelve Nonisolated Vertices

`vu-cover-ideal-linearity-question-through-twelve-nonisolated-vertices`

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
passes. The recorded run reports **53 checks, all passing**:

    VERDICT: ALL 53 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    89ef4a2c61628cbe6094ee36144188c2a26073b32fb086bdf39e4e68316c0daf

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE: what is re-derived for EVERY bipartite graph with at most 12 nonisolated vertices is Proposition 5 -- every edge maximising alpha(G-N[e]) has an endpoint in a maximum independent set. Theorem 1's own identity v(J^t) = v(J) + (t-1)tau is exhibited by direct computation only for n <= 9 (t = 1..20); for n = 10, 11, 12 it is not computed here, but follows from Proposition 5 through Proposition 4, which is proved analytically in the paper and machine-checked here only for n <= 9, and through Vu's cited local formula (1)-(3).
> NOTE The n <= 10 layer is scanned encoding by encoding: all 955674 retained encodings of Table 1, one at a time, by the literal route of equations (4)-(5).
> NOTE NOT RE-RUN one encoding at a time: the 464811101 retained encodings with n = 11 and n = 12. They are covered by the 110249245 orbit representatives the transposed scan has at those two values of n (108373735 of them at n = 12; the whole census up to n = 12 takes 110570825 representatives): every retained encoding of a cell (p,q) is a relabelling of the p-side of exactly one representative, and alpha(G), the union of the maximum independent sets and M(G) are isomorphism invariants, so no graph is left untested. The equality of the two routes is checked cell by cell wherever both are affordable, including two n = 12 cells.
> NOTE Table 1's own counts, including the n = 12 row (610303679 generated, 448777436 retained), are nonetheless recomputed three independent ways: the paper's closed forms, an unbounded-knapsack dynamic program, and prefix enumeration.
> NOTE NOT RE-DERIVED (cited algebra, not computation): Lemma 2 on isolated variables, the standard gradedness J^t = J^(t) of Herzog-Hibi-Trung, and Vu's local formula (1)-(3). Everything downstream of them is recomputed here, and the local formula is used only in the affineness check.
> NOTE NOT RE-RUN: the program that produced Table 1, which the paper keeps in an auxiliary archive available on request and does not reproduce -- this program is an independent reimplementation in Python, written from the paper's printed data alone. NOT VERIFIED BY MACHINE ANYWHERE ABOVE: that the encoding of Section 3 reaches every labelled bipartite graph without isolated vertices on 8 to 12 nonisolated vertices. That completeness is checked by exhaustive construction only for n <= 7 (77340 distinct labelled graphs, reached through 169656 (bipartition, edge set) pairs); for 8 <= n <= 12 it is argued from the encoding, not computed.
