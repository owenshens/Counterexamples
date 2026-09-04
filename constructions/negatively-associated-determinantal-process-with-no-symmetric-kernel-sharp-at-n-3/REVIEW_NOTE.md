# Review note

**Paper.** *An Explicit Negatively Associated Determinantal Point Process with No Symmetric Kernel,
Minimal at n = 3.*

**Files a referee has here.** `paper.tex` and `paper.pdf` (the paper, 8 sections and 10 references),
`verify.py` (the verification program), `verify.output.txt` (the recorded run of that program).
Nothing else is needed to read or to check the paper.

## 1. What the paper claims

Poinas, *On determinantal point processes with nonsymmetric kernels* (Electron. J. Probab. **31**
(2026), no. 94; arXiv:2406.03360), ends the subsubsection *Negative association* of his section
*Discussion and open problems* with an unnumbered prose question. The paper quotes that paragraph in
full in its Section 1 and settles its **final sentence** — *does a negatively associated DPP exist
that does not have the same distribution as a DPP with a symmetric kernel?* — in the affirmative,
already on three points, for negative association in the standard finite-ground-set form the paper
displays as (3).

The witness is the single rational matrix
`W = [[1/2, 1/10, 1/25], [1/10, 1/2, 1/10], [1/4, 1/10, 1/2]]`.
**Theorem 1** (Section 2) asserts (i) `W` is a DPP kernel — the eight atoms are `1071/10000`,
`1329/10000`, `1271/10000`, `1129/10000`, all strictly positive; (ii) `DPP(W)` is negatively
associated, and (3) is strict whenever `A, B` are nonempty and `f, g` are both nonconstant, all 30
ordered inequalities holding with minimum slack `71/10000`; (iii) no real symmetric and no Hermitian matrix
has the principal minors of `W`, hence `DPP(W)` is the law of no symmetric-kernel and no
Hermitian-kernel DPP.

Part (iii) runs through one invariant of the *law*. With `u = K_ij K_ji`, `v = K_ik K_ki`,
`w = K_jk K_kj`, `p = K_ij K_jk K_ki`, `q = K_ik K_kj K_ji`, `sigma = p + q` (display (4)),
**Lemma 2** gives `pq = uvw` with no hypothesis, so `sigma^2 - 4uvw = (p - q)^2`, and exhibits
`u` and `sigma` as functions of the principal minors alone. **Theorem 3** then says that for fixed
`u, v, w > 0` a real matrix with those invariants exists iff `sigma^2 >= 4uvw`, a real symmetric one
iff `sigma^2 = 4uvw`, a Hermitian one iff `sigma^2 <= 4uvw`. For `W` one has `u = v = w = 1/100`,
`p = 1/400`, `q = 1/2500`, `sigma = 29/10000` and `sigma^2 - 4uvw = 441/10^8 > 0`, which excludes both
classes. **Remark 4** re-derives the real symmetric half without Theorem 3: the off-diagonal moduli are
forced to `1/10`, and the eight sign patterns realise only `det S = 14/125` or `27/250`, never
`det W = 1129/10000`.

**Theorem 5** (Section 4) makes `n = 3` sharp: such a process exists for exactly the `n >= 3`, by
padding `W` with independent Bernoullis upward, and by a case check at `n <= 2` downward.
**Theorem 6** (Section 5) adds that in the 7 coordinates `(a_1, a_2, a_3; u, v, w; sigma)` the
negatively associated 3-point laws have nonempty interior while real-symmetric-kernel laws are confined
to the hypersurface `sigma^2 = 4uvw`, so Lebesgue-almost every negatively associated 3-point law is the
law of no real-symmetric-kernel DPP.

Two limitations are stated by the paper itself, in the abstract and in Section 7, and they matter for
how the result may be quoted. `DPP(W)` is **not** strongly Rayleigh, so the same existence question
with negative association replaced by strong Rayleigh is untouched; and non-symmetrisable determinantal
laws are not new — Section 6 records two published constructions, one in the asking paper itself, whose
three-point laws already have positive invariant, for neither of which negative association is
established.

## 2. What the program checks

`verify.py` is exact-rational and self-contained; `verify.output.txt` prints one `PASS` line per check
in ten numbered sections and closes with `VERDICT: ALL 71 CHECKS PASS`. Counts below are the `PASS`
lines of those sections; the pairing with the paper's statements is the paper's own, given in its
Section 8.

| transcript section | checks | statement checked |
|---|---|---|
| 1. the primary witness `W` (n = 3) | 8 | Theorem 1(i): the minors, exchangeability, the eight atoms with sum exactly 1 and minimum `1071/10000`, the atoms re-derived by a second independent `det(K - I_{A^c})` route, and `P(X = empty) = det(I - W)` so the `L`-ensemble exists |
| 2. negative association, brute-forced from the definition | 11 | Theorem 1(ii): exactly 30 ordered inequalities, 0 violations, minimum slack `71/10000`, the three slack classes `1/100`, `71/10000`, `129/10000`, and the source's necessary condition `K_ij K_ji >= 0`. The up-set enumerator is pinned to the Dedekind numbers `2, 3, 6, 20, 168, 7581` (OEIS A000372) and the test counts to `30, 348, 4560, 140058`; all 7780 generated families are independently re-tested for up-closure by a tester shown non-vacuous |
| 3. symmetric and Hermitian exclusion | 8 | Lemma 2 (`pq = uvw`; `u` and `sigma` from the principal minors), Theorem 1(iii) (`gap = 441/100000000`), and Remark 4 (forced squares `1/100`; the complete 8-pattern enumeration yields only `14/125` and `27/250`). The Hermitian identity (5) is confirmed on 1728 Gaussian-rational candidates with the forced moduli |
| 4. controls, in both polarities | 9 | none of the paper's statements; the program's own controls, including a sign-flipped kernel that is still a DPP kernel but fails NA, and a nonsymmetric kernel of gap exactly 0 carrying the symmetric control's law |
| 5. "Theorem B" | 6 | none; a closed-form slack formula compared against the brute force. The label is the program's |
| 6. the second `n = 3` witness `W'` | 3 | none; an object the paper does not discuss |
| 7. "Theorem C" | 5 | none; a `sigma`-interval. The label is the program's |
| 8. "Theorem D" | 10 | the padding and `n <= 2` parts of Theorem 5: `W (+) diag(1/3)` is a DPP kernel and passes all 348 NA checks, the atoms factorise as an independent superposition, the obstruction on `{1,2,3}` is inherited, a 450-kernel `n = 2` sweep, and vacuity at `n = 1`. The `n = 4` witness `K4` in the same section is not a statement of the paper |
| 9. the box certificate | 7 | none; the program's own object |
| 10. not strongly Rayleigh | 4 | the second item of Section 7: at `x = (-4,-4,-2)` and `(i,j) = (1,2)`, `D_12 f = (2735^2 - 7675*987)/10^8 = -19/20000 < 0`, with the symmetric control giving `1/2500 >= 0` at the same point. The conditional-NA checks beside it are the program's own |

Of sections 4, 5, 6, 7, 9 and the conditional-association checks the paper says in its Section 8:
"Those labels are the program's own and do not correspond to statements of this paper." Nothing above
them depends on them.

## 3. What the program does not check

* **The general Hermitian exclusion is a hand proof; the program is a bounded control.** The
  transcript's own section-3 heading records this — "the general Hermitian case is proved in the paper,
  not here" — and its last check repeats it. What runs is the identity `sigma^2 - 4uvw = -4 Im(p)^2` on
  **1728** Gaussian-rational Hermitian candidates with the forced moduli `|S_ij|^2 = 1/100`, plus the
  fact that none of them attains `W`'s gap or `sigma`. The quantifier "no Hermitian matrix" in
  Theorem 1(iii), and both realisability (converse) directions of Theorem 3, are proved by hand only.
  The **real symmetric** exclusion, by contrast, is a complete enumeration: the moduli are forced and
  `1/100` is a rational square, so the eight sign patterns exhaust the case.
* **The `n <= 2` impossibility of Theorem 5 is proved by hand.** The run's closing scope paragraph says
  so in those words: the sweep of 450 rational NA `2x2` DPP kernels whose pair invariant is a rational
  square is "a corroboration, not the proof". (300 further kernels with `u < 0` were correctly rejected
  by the NA test.)
* **Theorem 6 is not computed.** The program certifies only that the eight atom inequalities, the 30 NA
  inequalities, `u, v, w > 0` and `sigma^2 - 4uvw > 0` hold *strictly* at the law of `W`. The step from
  that to nonempty interior in the 7-dimensional space, and the Lebesgue-nullity of the hypersurface
  `sigma^2 = 4uvw`, are the hand argument of Section 5. The box certificate of transcript section 9 is a
  statement about a box of matrix *entries* and is not one of the paper's claims.
* **Two groups of the paper's numbers are outside the 71 checks**, and the paper lists them itself at the
  end of Section 8 rather than leaving a reader to hunt for them: (a) all arithmetic of Section 6 — the
  triple invariants of Poinas' coupling (`u = w = b^2`, `v = K_11(K_11 - p_1)`, gap `(b^2 p_1)^2`) and
  the Borodin–Diaconis–Fulman numbers (`det K = 0` against a symmetric competitor's `-d^3`, the external
  checksum `P(0,0,0) = 1 - 3d + d^2` from their Theorem 4.1 and Example 6, and `u = w = 1/12`, `v = 0`
  for the descent process); the program's "source's own kernel" control is the uniform-on-even-subsets
  kernel, a *different* object. (b) The Hermitian counterexample in the first item of Section 7 (atoms
  `11/100` and `13/100`, NA slack `1/100`), which the program never evaluates. Both groups were checked
  by hand; neither is load-bearing — the first is attribution, the second only weakens Theorem 6.
* **Transcribed from cited sources, not recomputed:** the quoted question and the block-diagonal
  sentence of Poinas, the closure of negative association under independent superposition (Joag-Dev and
  Proschan, property P7) used in Theorem 5, Brändén's criterion, the Borodin–Diaconis–Fulman checksum
  above, and the Dedekind numbers against which the up-set enumerator is pinned.
* **Further scope, carried over from the run's own closing `NOT RE-RUN / SCOPE` paragraph:** the failure
  of strong Rayleigh is established at **one** explicit point, which suffices to refute the criterion,
  and no sweep over `R^3` is performed; conditional NA is tested with one coordinate fixed at a time and
  the full CNA+ hierarchy is not enumerated; "Theorem C" is verified on a finite rational sample of its
  `sigma`-interval, the interval statement itself being proved in the paper; the box certificate is a
  *sound lower bound*, so the negative bounds it reports at radii `1/50` and `1/20` are not evidence of
  a bad matrix, and `1/100` is a certified lower bound on the true radius, not the radius; and Poinas'
  first question at `n >= 4`, the continuous-space version of the question, and the facet description of
  the 7-dimensional NA region are open, with nothing here bearing on them.
* **The paper does not rest on the program at all.** Section 8 opens by stating that every object is
  printed in full in exact rationals, so the verification can be redone with a pen: the eight atoms by
  (2), the three slack classes, the six-term determinant, Lemma 2, and Theorem 1(iii) in four
  multiplications and one subtraction, with Remark 4 replacing even Theorem 3.

## 4. How to check it

```sh
shasum -a 256 verify.py
python3 verify.py
```

The first command prints

    14121aa907789e10f8610ea32ecdefd759f003187f8ad5d7e33e2afea87f714c

which is exactly the value on the `sha256:` line in the header of `verify.output.txt`, so the shipped
transcript can be paired with the shipped program. That header, which accompanies the output rather than
being printed by the program, also names `verify.py` and `Python 3.9.25`.

`verify.py` needs Python 3.9+ and the standard library only (`sys`, `fractions`, `itertools`); it takes
no arguments, reads nothing from disk — the matrices are typed into it as printed in the paper — and
uses `Fraction` and integer arithmetic throughout, with no float deciding anything. It emits one
`PASS <name> [detail]` line per check, then its closing `NOT RE-RUN / SCOPE` paragraph, then the
verdict, and exits 0 iff every check passed; the recorded run ends `program exited with status 0`
under `VERDICT: ALL 71 CHECKS PASS`. `paper.pdf` is the compiled `paper.tex`.
