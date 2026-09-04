# Review note

Paper: *A counterexample, with a hyperbolic summand, to a printed subadditivity question for
topological volume* (`paper.tex`, `paper.pdf`; 5 sections, 4 pages).

The folder contains exactly five files: `paper.tex` and `paper.pdf` (the paper), `verify.py`
(an arithmetic program), `verify.output.txt` (the recorded run of that program), and this
note. Nothing else is needed to read or to check the paper.

## 1. What the paper claims

Kegel, Ray, Spreer, E. Thompson and Tillmann, *On a volume invariant of 3-manifolds*
(reference [4] of the paper), define the topological volume `vol_t(M)` of a closed orientable
3-manifold as the infimum of `vol_h(M \ L)` over hyperbolic links `L ⊂ M`, the empty link
included — equation (1) of §1 — and in the subsection *Connected sums* ask whether
`vol_t(M # N) ≤ vol_t(M) + vol_t(N)` for all closed orientable `M`, `N`. §1 quotes that
Question verbatim and observes that as printed it carries exactly two hypotheses, closed and
orientable.

The paper answers the Question **no**, at one pair. **Theorem 3** (§3, *The counterexample*):

> `vol_t(W # W) ≥ 2v_0 > 1.886 ≥ vol_t(W) + vol_t(W)`,

`W` the Weeks manifold and `v_0` the volume of the regular ideal hyperbolic tetrahedron; so
the Question of §1 is false at `(W, W)`. The proof is three lines from the two statements of
§2, *Two lemmas*: **Lemma 1**, a closed orientable non-hyperbolic 3-manifold has
`vol_t ≥ 2v_0`, from Cao–Meyerhoff's minimum volume for *orientable cusped* hyperbolic
3-manifolds ([2]); and **Lemma 2**, `M # N` with neither summand `S^3` is reducible and hence
not hyperbolic, by the Sphere Theorem. The two numerical inputs are equation (2),
`vol_t = vol_h` on hyperbolic manifolds, quoted from [4], and equation (3),
`vol_h(W) ≤ 0.943`, quoted from the abstract of Gabai–Meyerhoff–Milley ([3]).

Both summands are hyperbolic, and §4, *What is not settled*, says in the paper's own words
that the reading of the question in which the summands are **not** hyperbolic "is untouched".

## 2. What the program checks

`verify.output.txt` records **38 checks, all passing** — closing line `VERDICT: ALL 38 CHECKS
PASS`, program exit status 0. The header records the program name, its SHA-256 and the Python
version (3.9.25); the run reports decimal working precision 80 with "all decisions in exact
rationals". Blocks, as the transcript labels them:

* **Step 1, 8 checks** — the two constants from their definitions: π by two Machin-like
  formulas, negligible Clausen-series truncation bounds, the Bernoulli recurrence against
  `B_2..B_10`, then `2v_0 = 2.029883212819307250042405` and
  `v_0 = 1.014941606409653625021203` to the 25 digits §1 prints, cross-checked by the Clausen
  duplication identity. Two of the eight concern values the paper does not print: the source's
  own rounding `2.02988321281931` and `4v_0 = 4.0597664256386145`.
* **Step 2, 5 checks** — transcribed rows and level sizes of Table 1 of [4]
  (`11+13+13+6+14+4+10+5 = 76`; level-one round-up `2.02989`; level 2 volume
  `2.56897060093671 > 3.07 - 0.943 = 2.127`). The paper states none of this.
* **Step 3, 3 checks** — `2 × 0.943 = 1.886`, and `2v_0` exceeding that sum with
  `gap ≥ 0.143883212819307`: precisely the two arithmetic steps in the proof of Theorem 3. The
  third check is a second route through a "census reach" of `3.07`, which the paper does not
  use.
* **Step 4, 11 checks** — ten pairs `W # N` on that second route
  (`0.943 + 2.02989 = 2.97289 < 3.07`, `gap ≥ 0.09711`), plus a count of eleven pairs. Not in
  the paper.
* **Step 5, 4 checks** — one of them re-decides
  `v_0 = 1.014941606409653625021203 < 1.886`, which is the content of the Remark after
  Lemma 1 (the Gieseking manifold). Another records that subadditivity *holds* at the source's
  own worked example, `vol_t(L(2,1) # L(3,1)) = 2.66674478344906 < 5.23571538438577` — the
  strict inequality displayed under *Adjacent work* in §4, for which the paper prints no
  decimals. The other two are the transcribed sum
  `2.66674478344906 + 2.56897060093671 = 5.23571538438577` and the observation that it exceeds
  `3.07`, so the second route "never fires there".
* **Step 6, 4 checks** — boundary cases: `4v_0 = 4.0597664256386145 > 3.07`,
  `3 × 0.943 = 2.829 > 2v_0`, one check on the floating-point value `1.014941606409654`, and
  one recording that doubling that value lands on `2v_0`, an equality the Question's `≤`
  satisfies.
* **Step 7, 3 checks** — by-product ratios (`1.07`, `1.62`, `0.510`), which the transcript
  itself heads as "not themselves printed in the paper".

So **9 of the 38 checks bear on numbers the paper prints**: six of Step 1 on the two 25-digit
expansions of §1, two of Step 3 on the proof of Theorem 3, one of Step 5 on the Remark. §5,
*Verification*, lists the remainder as material the paper "does not state" and for which the
checks "are not evidence".

One naming caution: the transcript's check names read `X0_route_A_...`, `X0_route_B_...` and
`..._under_theorem_B`, and one detail line says "theorem B never fires there". Those are the
program's own labels; the paper numbers a single Theorem 3, which is the `route_A` check, and
states no "theorem B".

## 3. What the program does **not** check

**Theorem 3 is a hand proof and the program is a control.** §5 puts it in those terms: "The
result above is a hand proof … Nothing in §§1–3 needs a computer." The program re-derives
numbers and re-decides comparisons; it decides no topology. Carrying over the transcript's
closing `NOTE SCOPE` block:

1. The topological content is quoted from the literature and is "not machine-checkable here":
   `vol_t = vol_h` on closed hyperbolic manifolds (equation (2)), Cao–Meyerhoff's `2v_0` floor
   for **orientable** cusped manifolds, Gabai–Meyerhoff–Milley's `vol_h(Weeks) ≤ 0.943`, the
   Sphere Theorem, and Kneser–Milnor uniqueness.
2. It "verifies no irreducibility of any manifold named". This is the unverified premise of
   the second route: each of the ten Step 4 checks is followed by a `NOTE` saying that an
   imported census theorem with reach `3.07` plus irreducibility of the row — "neither stated
   nor checked here" — would be needed to make the inequality topological. §5 of the paper
   says the same and adds that Theorem 3 uses none of it.
3. The values `0.943`, `3.07`, `2.02988321281931`, `2.66674478344906` and `2.56897060093671`
   are hand-typed literals, transcribed rather than recomputed.
4. One value, `1.014941606409654`, is a hand-typed floating-point literal that is **not
   interval-certified**; it is used only in Step 6. §5 calls it a floating-point SnapPy
   measurement and flags it the same way; no such measurement is part of this folder, and
   Theorem 3 does not rest on it.
5. No volume of the Weeks manifold is used anywhere except through the bound `0.943`; §3 says
   the paper never uses a decimal expansion of `vol_h(W)`.

Limits the paper sets for itself in §4, which no program addresses: the case where **both**
summands are non-hyperbolic is untouched, and the paper explains why this route cannot reach
it (Lemma 1 already puts both terms at `≥ 2v_0`, so a counterexample there needs a lower
bound on `vol_t(M # N)` above a sum of at least `4v_0`, and nothing in the paper supplies
one); nothing in the paper supports or refutes `vol_t(M # N) ≤ vol_t(M) + vol_t(N) + 2v_0` or
any other repaired form; the paper exhibits "the one pair `(W, W)`" and counts no others; and
it settles the Question as printed and makes no claim about what was meant. Three
bibliographic limits are stated too: [1] is not open access and was consulted only through
works restating its definition, including [4]; [3] is cited at `arXiv:0705.4325` only, with no
journal coordinates because none were checked; and [5] carries a caution that the coordinates
[4] gives for it do not resolve.

## 4. How to check it

```sh
python3 verify.py          # one line per check; exit status 0 iff all pass
shasum -a 256 verify.py
```

The second command should print the digest of the shipped `verify.py`:

    f2a30cc92abe4c19aced27467f9ee5a49491dcc840cec2f5c4a360bc62f21ed7

The header of `verify.output.txt` carries this same SHA-256 beside the program name, so the
transcript and the program can be paired. Standard library only: no third-party package, no
external data file, and no input beyond the literals typed at the top of `verify.py`.
