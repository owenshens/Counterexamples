# Two Zonotopes with the Same Arithmetic Matroid and Different Graded Ehrhart Series

`two-zonotopes-with-equal-arithmetic-matroid-and-different-graded-ehrhart-series`

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

Python 3.9 or later, standard library only (`fractions`, `itertools`, `math`): no third-party
package and no external data file. Exact integer and rational arithmetic throughout, no floating
point, no randomness and therefore no seed. It runs in a few seconds on one core. The program
prints one line per check and a closing verdict, and exits 0 only if every check passes. The
recorded run reports **48 checks, all passing**:

    VERDICT: ALL 48 CHECKS PASS

Everything it consumes is typed into the block marked `THE PAPER'S CLAIMS` at the top of the
file: the two generator matrices, the two lattice-point sets (3) and (4), the eight cubics of
Table 1 with their printed values at the omitted points, the coefficient lists of the two graded
series for `m <= 4`, and the control values quoted from the literature. Everything else it
derives. In particular it obtains the multiplicity function by integer Smith normal form rather
than by the gcd/determinant shortcut; the lattice-point sets twice, by an all-integer
support-function test and by an exact-`Fraction` inverse-matrix test, cross-checked against
Pick's theorem; and each graded series twice, once by the rank differencing of Lemma 2 and once
from an explicit basis of the vanishing ideal in each degree and its top forms — the second route
also verifying that `x*W_{j-1} + y*W_{j-1}` sits inside `W_j`, so that Lemma 2's identity is
checked rather than assumed.

It also runs the controls in both directions: it reproduces the four published graded values for
planar zonotopes that exist (two from the target paper, two from the ancillary file of Reiner and
Rhoades), returns *equal* at `q = 1` on the same pair for every `m <= 4`, stays silent on all 32
lattice equivalences of the first zonotope, and confirms the two anti-collapse predictions
`(1,0),(3,5) -> Z_2` and `(1,0),(4,5) -> Z_1`.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    af52f9ec93cba098aaa0c4dcdd06c9ec9cb316544e28b60f3ccb58028258301b

## Scope

The program **re-derives** the claims it checks; it does not merely confirm the exhibited object.
The only things it is given are the two generator matrices and the assertions listed above. The
lattice-point sets, the multiplicity function, the arithmetic Tutte polynomial, every graded series,
the determinant census and the family of pairs are all computed from the generators, and the listed
assertions enter only as comparison targets — so corrupting a generator, a point set, a coefficient
list, a cubic or a cubic's stated value makes the program report the failed check and exit non-zero
rather than pass.

The program's own statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the SECOND clause of Question 7.3. Nothing above decides whether some invariant
> strictly between the arithmetic matroid and the full lattice class supports a generalisation of
> the paper's Proposition on Ehrhart polynomials; that stays open.

> NOT RE-RUN: minimality outside the cell d = 2, n = 2. The census is exhaustive over GL_2(Z)
> classes of PAIRS only; no search over three or more generators is performed here, so nothing
> above excludes a smaller or differently shaped counterexample at n >= 3.

> NOT RE-RUN: the independent top-form route to dim Orb_j (checks `series-second-route` and
> `gr-ideal-equals-top-forms`) covers m = 1 and m = 2 only; the m = 3 and m = 4 series rest on
> rank differencing alone.

> NOT RE-RUN: bibliographic facts. The line numbers quoted for arXiv:2603.07873, the numbering of
> Question 7.3, the ancillary-file line numbers of arXiv:2407.06511 and every claim about the
> literature are read from sources, not computed here; only the NUMERICAL values taken from those
> sources are recomputed, as the `control-*` checks.

> NOT RE-RUN: the infinite family is checked for the listed finite ranges of k, not proved for all
> k; the proof for all k is the argument in Section 4 of the paper.

The program checks more than the paper states: the paper confines itself to the pair at `m = 1`, and
the census over determinants, the family indexed by odd `k`, and the series for `m = 2, 3, 4` are
left in the program alone. The paper makes no novelty claim.

Finally, the decisive claim needs no program at all: Section 3 of the paper proves the separation
by hand from the two eight-point sets printed in Section 2, and a referee who prefers not to run
code can check it there.
