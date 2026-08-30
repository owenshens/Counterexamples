# The First r+1 Columns of the Cubic-Dimension Spectral Sequence of a Hypercube Vietoris–Rips Complex

`the-first-r-columns-of-the-hypercube-rips-spectral-sequence-are-the-cube-chain-complex`

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
passes. The recorded run reports **144 checks, all passing**:

    VERDICT: ALL 144 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    41c0edbec3e7ab512c4b31a0b23bf22f5f3dbb322b1eeb0ac60c504efaaf64ec

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the two results quoted from Galetto-Montano-Wellner and not reproved in the paper -- the n-freeness of the p-th column (source lines 930-939, label pro:1) and the induced description of F_i with its determinant character (lines 1029-1047, lem:1). Every check above is therefore made at n = p, which is exactly the case those results reduce to; nothing here verifies the reduction itself.
> NOT RE-RUN: the standard topology the proof cites -- the long exact sequence of a pair, the nerve theorem, the Mayer-Vietoris spectral sequence of a cover, and the Kunneth formula for a join. What is checked is that their CONCLUSIONS agree with direct computation at p <= 4 (check two-routes-agree-p*), not the theorems.
> NOT RE-RUN: Lemma 4.1 itself. The marginal map argument is a homotopy of topological pairs and is not modelled here; what is checked is its numerical consequence, degree by degree, at p <= 4, plus the H_p-character of the answer.
> NOT RE-RUN: any cell with p = r >= 5. That would need 2^32 subsets of V(Q_5). The proof of the paper is uniform in n, r and p; the machine evidence here stops at p = 4, and the equivariant character check stops at p = 5 (nerve) and p = 4 (chain level).
> NOT RE-RUN: the flag-transitivity and unit-coefficient steps of D4 at general n. They are checked only at n = 6 and p <= 4.
> NOT RE-RUN: the characteristic-zero labelling {lambda;mu} of H_m-irreducibles. This program computes characters as integer class functions and compares them with det and with the trivial character; it does not decompose any representation into the source's labelled irreducibles.
> NOT RE-RUN: the two off-range facts of Section 6(3) as statements about H_m-isomorphism TYPES. Their dimensions and Lefschetz characters are reproduced; the identifications {1^3;0}+{0;1^3} at (3,2) and {4;0} at (4,3) are quoted from the source.
> NOT RE-RUN: the literature. This program performs no search, fetches nothing, and checks no citation.
> [24.3s] 144 checks ran, 0 failed
