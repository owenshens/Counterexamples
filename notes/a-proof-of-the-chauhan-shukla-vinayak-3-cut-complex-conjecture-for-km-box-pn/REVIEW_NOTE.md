# A Proof of the Chauhan–Shukla–Vinayak Conjecture on the 3-Cut Complex of K_m P_n

`a-proof-of-the-chauhan-shukla-vinayak-3-cut-complex-conjecture-for-km-box-pn`

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
passes. The recorded run reports **155 checks, all passing**:

    VERDICT: ALL 155 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    fc678cb04c8cb7f27b53b1f89a7a2805e8ea239b3f5cd4b81697d09ee3bd23a4

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE: what this program does NOT cover
> NOT RE-RUN: the theorem itself.  The paper's claim is a statement about ALL m, n >= 3 and it is
> settled by the hand proof of Section 4, not by any census.  Every cell this program computes is a
> finite check of that proof's intermediate claims and of its two closed forms.  A program cannot
> verify an infinite family, and nothing below should be read as if it had.
> NOT RE-RUN: the large census.  The discovery run computed the reduced homology of D* on 66
> distinct cells, out to (10,10) at 100 vertices, on AWS instances.  This program re-derives the
> same quantities on the much smaller list printed above, chosen so that a referee can run it in
> under a minute on a laptop with nothing installed.  The cells (10,9) and (10,10) appear here only
> in the f-vector and Euler-characteristic check, not in a homology computation.
> NOT RE-RUN: Delta_3 computed directly, with no Alexander duality, at any cell beyond (3,3),
> (3,4) and (4,3).  Those three are the three smallest cells of the source's table and they are
> where the duality step is pinned to an independently computed answer; larger cells are checked
> through the dual only.
> NOT RE-RUN: the exhaustive-subset enumeration of D* beyond mn <= 16.  For larger cells the dual
> is enumerated by the incremental rule instead, and the two enumerations are shown to agree
> wherever both are affordable.
> NOT RE-RUN: anything bibliographic.  The wording of Conjecture 4.4, its position at lines 871-880
> of the source's 80,146-byte main .tex file, the numbering that follows from the shared theorem
> counter, the contents of the source's table, the two boundary lines quoted from Bayer et al.
> (2024) and Cut Complexes II (2025), and the prior-art search, were established by fetching and
> reading those sources.  This program checks mathematics, not provenance, and fetches nothing.
> NOT RE-RUN: the total k-cut complex Delta_3^t, and with it the source's sibling Conjectures 4.3,
> 4.5 and 4.6.  Nothing in this program touches the total operator.
