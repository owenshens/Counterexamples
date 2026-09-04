# Referee note: the Liu--Li predicted degree sequence at (n, Delta) = (8, 5)

`liu-li-2008-degree-sequence-conjecture-fails-at-n-8-maximum-degree-5`

Besides this note the folder holds four files: `paper.tex` and `paper.pdf` (the paper), `verify.py`
(the program) and `verify.output.txt` (the program's recorded run). No file outside the folder is
named below, and the paper's argument rests on nothing outside it.

## 1. What the paper claims

Let `G(n, Delta)` be the set of graphs attaining the maximum spectral radius among all connected
nonregular graphs on `n` vertices with maximum degree `Delta`. The target, quoted verbatim in
Section 2 ("The conjecture, verbatim"), is Conjecture 2.12 of B. Liu and G. Li,
*Electron. J. Linear Algebra* **17** (2008), 54--61, in the forward-only restatement the paper takes
from two later sources: for `3 <= Delta <= n-2`, every `G` in `G(n, Delta)` has degree sequence
`(Delta, ..., Delta, delta)` with `delta = Delta-1` when `n*Delta` is odd and `delta = Delta-2` when
`n*Delta` is even.

**Theorem 6** settles that statement in one cell, `(n, Delta) = (8, 5)`: it is false there, no member
of `G(8,5)` carrying the predicted degree sequence `(5,5,5,5,5,5,5,3)`, written `(5^7,3)`. Three exact
facts carry it, none of them floating point (Section 5, "The refutation", says "Nothing floating point
enters it"). **Lemma 2**: the exhibited 18-edge graph `H` on 8 vertices, degree sequence `(5^7,1)`, is
connected, nonregular and of maximum degree exactly 5, so it lies in the family over which `G(8,5)` is
defined. **Lemma 3**: `lambda_1(H) > 49/10`, from an explicit integer vector of exact Rayleigh quotient
`61272483059723/12499997522989`, and independently because the characteristic polynomial (1) satisfies
`p_H(49/10) = -27065510199/100000000 < 0`. **Lemma 4**: every connected graph on 8 vertices with degree
sequence `(5^7,3)` has `lambda_1 < 97/20`, by Sylvester's criterion on the integer matrix `97I - 20A`
applied to the three isomorphism-class representatives `C1`, `C2`, `C3` of Section 4.2, cross-checked
against the polynomials (2) via `det(97I - 20A) = 20^8 * p(97/20)`. Since `97/20 < 49/10`, a maximiser
cannot carry the predicted sequence.

Two disclosures the paper makes in its own voice should be read first. **No priority is claimed for
refuting the conjecture** (Section 1, "Prior work"): Huang, Liu and Yang, arXiv:2411.17371, determine
`G(n, n-2)` (their Theorem 7), yielding counterexamples at every even `n >= 8` with `Delta = n-2`, and
`G(n, n-3)` for `n >= 59` (their Theorem 8), yielding counterexamples at every even `n >= 60`; the cell
`(8,5)` has `Delta = n-3` and `n = 8`, so it falls in the gap between those theorems, while L. Liu,
*J. Combin. Theory Ser. B* **169** (2024), proves the conjecture for `Delta = 3` and `Delta = 4`. The
stated increment is this one cell with an exact, solver-free integer certificate over the complete
predicted class. **The witness is not a discovery** (Remark 1): the `(5^7,1)` class on 8 vertices is
2,520 labeled graphs forming a single isomorphism class with `|Aut| = 16`, so `H` is forced by its
degree sequence and is the `n = 8` member `H_1(8)` of arXiv:2411.17371; what is proved is a comparison,
not a property of `H`.

## 2. What the program checks

`verify.output.txt` records **82 checks, all PASS**, closing `VERDICT: ALL 82 CHECKS PASS` with exit
status 0. Counting the `PASS` lines under each labelled step of that transcript:

* **Step 1, 6 checks** -- Section 3, the cell: `3 <= 5 <= 6`; `n*Delta = 40` is even, so `delta = 3`
  and the predicted sequence is `(5^7,3)`; 19 edges; the complement has 9 edges with degrees
  `(4,2^7)`; parity leaves only `delta` in `{1,3}`.
* **Step 2, 14 checks** -- Lemma 2 and the structure of Section 4.1: 18 distinct edges, degrees
  `(1,5^7)`, maximum degree exactly 5, nonregular, connected; both graph6 strings, `GI[z~g` and
  `GNz~s?`, decode to and re-encode from the printed edge lists; the core is `K_7` minus exactly the
  four edges `15, 23, 47, 67`, which span `P_3 + 2K_2` with the pendant on the `P_3` centre; the two
  labellings are cospectral; the characteristic polynomial is `[1,0,-18,-32,5,32,12,0,0]`, i.e. (1).
* **Step 3, 9 checks** -- Lemma 3 by both routes: the Rayleigh vector is a nonzero integer vector, its
  exact quotient is the fraction printed in the paper and exceeds `49/10`; `p_H` is monic and
  `p_H(49/10) = -27065510199/100000000 < 0`. Also recorded: `p_H(5) = 16800` and
  `p_H(99/20) = 197988807248001/25600000000 > 0`.
* **Step 4, 26 checks** -- Section 4.2 and Lemma 4: for each of `C1`, `C2`, `C3`, 19 edges with degree
  sequence `(5^7,3)` and connected, graph6 round-trip, characteristic polynomial as in (2), all eight
  leading principal minors of `97I - 20A` positive and equal to the lists printed in the proof of
  Lemma 4, and the cross-check `det(97I - 20A) = 20^8 * p(97/20)`; then exact rational brackets for
  `lambda_1(H)` and the three `lambda_1(C_i)`, giving `lambda_1(C1) < lambda_1(C2) < lambda_1(C3)`
  among the three representatives.
* **Step 5, 9 checks** -- Remark 5, the whole-class fact by the program's own enumeration: generating
  the `(5^7,3)` graphs through their 9-edge complements of degrees `(4,2^7)` gives exactly 21,000
  labeled members, all 21,000 connected, with Sylvester's criterion on `97I - 20A` holding for each
  one individually, **0 violations** (0 also for the weaker `49I - 10A`, which certifies only
  `lambda_1 < 49/10`). It recovers exactly three complement types, of orbit sizes 840, 10080, 10080
  summing to 21,000, reproduces the hand orbit formulas of Section 4.2, and records which of `C1`,
  `C2`, `C3` realises which type.
* **Step 6, 6 checks** -- Remark 1: exactly 2,520 labeled `(5^7,1)` graphs, all connected, a single
  isomorphism type, `2520 = 8!/16`, every one of the 2,520 with `det(49I - 10A) < 0` hence
  `lambda_1 > 49/10`, and `H` in that class.
* **Step 7, 8 checks** -- controls in both polarities: `K_8` and the disconnected `K_6 + K_2` must be
  refused at `97/20` and are (the latter has maximum degree 5, is nonregular and has `lambda_1 = 5`,
  so it is the anti-control for a broken connectivity filter); `P_8`, `C_8`, `K_{1,7}` must be
  certified and are; the certifier must refuse `H` itself and does; the lower-bound test must stay
  silent on `C3` and does.
* **Step 8, 4 checks** -- the assembly of Theorem 6: `97/20 < 49/10`; `H` beats every graph of the
  predicted degree sequence; hence no maximiser at `(8,5)` carries it; plus a robustness check that a
  graph of maximum degree at most 4 cannot compete (`lambda_1 <= 4 < 4.85`).

The program's header lines state "exact integer / exact rational arithmetic only; standard library
only"; it imports only `itertools`, `sys` and `fractions.Fraction`, opens no file and uses no network.
Characteristic polynomials are by Newton's identities on integer traces, minors and determinants by
fraction-free (Bareiss) elimination; both raise rather than round on an inexact division. The recorded
run reports Python 3.9.25; the paper's "Verification" section says Python 3.9+.

## 3. What the program does not check

The transcript's closing `NOTE SCOPE` block states most of this itself; it is carried over here.

* **The program does not parse the paper.** It re-derives quantities about objects transcribed into
  itself -- the two edge lists of `H`, those of `C1`, `C2`, `C3`, the graph6 strings, the integer
  Rayleigh vector, the polynomials, the minors of `97I - 20A`, the evaluations at `49/10`, `5` and
  `99/20`, and the four spectral brackets. `NOTE SCOPE` marks exactly this `NOT RE-RUN`: completeness
  of that transcribed list against the paper is not verified, and the paper's "Verification" section
  says the same. So every check named `..._matches_the_paper` compares against a constant held in the
  program, and for five of them the quantity is not printed in the paper at all (`p_H(5) = 16800`,
  `p_H(99/20)`, and the three values of `p_{C_i}(49/10)`), as are the four brackets' endpoints; for
  these the check is internal consistency, not agreement with the paper.
* **The statement of the conjecture is not re-derived.** The program's own opening lines say the
  Liu--Li statement "is NOT re-derived here and is taken from the paper, Section 2", the predicted
  sequence at `(8,5)` being assumed to be `(5^7,3)`; fidelity of the Section 2 quotations to the
  cited sources is a hand matter for the referee. The Step 4 ordering check repeats in its own text
  that "the conjecture's own statement is not verified here".
* **Lemma 4 as written rests on a hand classification.** The completeness argument of Section 4.2 --
  that the class has exactly the three isomorphism types -- is by hand. Remark 5 says so and says the
  Step 5 enumeration removes that dependency: the machine route replaces the hand completeness
  argument rather than checking it. Both routes are available to the referee, and Theorem 6 needs
  only one of them.
* **Which graph attains `lambda(8,5)` is not certified, and is not needed.** `NOTE SCOPE` and
  Section 6 ("Scope") both say so: optimality of `H` in `G(8,5)` is not proved (Remark 1), the paper
  does not assert as a certified statement that `H` lies in `G(8,5)`, and the replacement conjecture
  mentioned in Remark 7 is not proved either.
* **No cell other than `(8,5)`.** `NOTE SCOPE` states that nothing here concerns `Delta = n-2`,
  which belongs to Huang--Liu--Yang; Section 6 adds that at `Delta = 5` the cells `n = 10` and
  `n = 12` are untouched.
* **Automorphism orders are not computed as group orders.** Section 4.2 prints orders 48, 4, 4. The
  run records orbit sizes 840, 10080, 10080, and for `C1` only that orbit 840 is consistent with
  `|Aut| = 8!/840 = 48`; the two orders 4 are not checked. `|Aut| = 16` for the witness class is
  likewise recorded as `2520 = 8!/16`, from the orbit size.
* **One residual priority risk, stated in Section 6 and unresolved.** Table 1 of B. Liu, J. Shen and
  X. Wang, *J. Combin. Theory Ser. B* **97** (2007), 1010--1018 -- the table of `lambda_1`-extremal
  graphs the 2008 source explicitly consults ("by checking the Table 1 of [7]") -- could not be
  obtained; if it reaches `n = 8`, this cell was settled in 2007. The evidence against is indirect:
  Liu, Huang and You, *Electron. J. Linear Algebra* **18** (2009), Remark 2.11, quoted in Section 3,
  says the conjecture "need[s] to be verified for `n = 5,6,7`" and lists exactly six cells stopping at
  `n = 7`, and both 2024 sources treat `Delta >= 5` as open. Nothing in the program bears on this.

## 4. How to check it

```sh
shasum -a 256 verify.py
python3 verify.py
```

The digest of the shipped `verify.py`, computed from the file in this folder, is

    b001f1e24101a7c6c7c7d33f9b9380700317d3b9850ae0aa5185f6a90299d0db

and the header of `verify.output.txt` carries that same SHA-256 beside the program name, so the
transcript and the program can be paired before reading either. The run takes no argument and no input
file; it should reproduce the 82 `PASS` lines, close with `VERDICT: ALL 82 CHECKS PASS`, and exit 0.
