# A Counterexample to Cartwright–Chan's Principal 6 6 Criterion for Tree Rank

`a-counterexample-to-the-cartwright-chan-principal-6x6-tree-rank-criterion`

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
passes. The recorded run reports **110 checks, all passing**:

    VERDICT: ALL 110 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    868422cc0363c2e1683529417986b3011151373341baade7c126286c985add06

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> Scope of this re-run
> NOTE NOT RE-RUN: the uniform theorem is proved in the paper for every n >= 7; this program confirms it only for n = 7,...,14, so n >= 15 rests on the hand proof alone.
> NOTE NOT RE-RUN: whether some local size k >= 7 does certify tree rank 2. The paper leaves that open and this family cannot settle it, because for k >= 7 the principal k x k submatrices of the uniform family include M itself.
> NOTE NOT RE-RUN: the Pachter--Sturmfels notion of tree rank named in the same bullet of the source. It is a different notion of mixtures, nothing here touches it, and no computation could -- their text was not read.
> NOTE NOT RE-RUN: the wording of the conjecture in the printed Combinatorica text. The statement is quoted from arXiv:0912.1411v1, and no computation can check a transcription.
> NOTE NOT RE-RUN: dissimilarity matrices that are not 0/1. The source proposition used here as the decider applies only to the 0/1 sub-family, and the census of Steps 10-11 is a census of that sub-family alone; a real counterexample of smaller order is not excluded.
> NOTE NOT RE-RUN: anything about the tropical Grassmannian directly. Every tree rank above is obtained through the source paper's own combinatorial proposition, which this program takes as given rather than reproving.
