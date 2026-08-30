# A Counterexample to the Singh--Singh Width-2 ROABP Modular-Stability Proposition

`a-counterexample-to-the-singh-singh-width-2-roabp-modular-stability-proposition`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

What the paper settles is that two printed statements of Shalender Singh and Vishnupriya Singh,
*An Algebraic Rigidity Framework for Order-Oblivious Deterministic Black-Box PIT of ROABPs*,
arXiv:2602.13449v1 (13 February 2026, 48 pp, PDF-only), are false: the numbered, unconditional
**Proposition 4.7.3** on printed p. 28 ("Conjecture 4.7.2 holds for width-2 diagonal ROABPs",
carrying a two-sentence proof), and the printed quantifier form of **Conjecture 4.7.2** on printed
pp. 26--27. The witness is `C = x_1 - x_r` on `n >= r` variables, a width-2 diagonal ROABP of
individual degree 1 whose bad set is all of `Z_r^*`; the mechanism is Fermat's little theorem.

What the paper does **not** settle is the conjecture the authors *describe* one page earlier, in
their own Remark on printed p. 27: `C` of width at least 3 fixed in advance with non-degenerate
coefficient support density. Nothing in this folder bears on that statement, and it remains open.
The witness object itself is not new either --- it is the source's own Example 4.7.3 (printed
p. 27) at `S = {1}`, `S' = {r}` --- and Section 3 of the paper says so; what is new is the reversed
order of quantifiers and the resulting count `|B_r(C)| = r - 1`. See `## Scope` below.

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

Python 3.9 or later, standard library only: no third-party package, no network and no external
data file. Every object it consumes is a literal in its own source, and each of those literals is
printed in the paper --- the substitution `Gamma_g : x_i -> lambda^(g^i)` of the source's
Definition 4.7.1, the two witness families `x_1 - x_r` and `x_1 - x_((r+1)/2)`, the primes
31, 131, 1039, 4099 and 10007, the class parameters `w = 2`, `d = 1`, the source's own illustrative
threshold `(wd)^7 = 128`, `epsilon = 1/10`, and the six primes 1009, 1013, 1019, 1021, 1031, 1033.
Arithmetic is exact integer arithmetic; there is no floating-point decision anywhere in the
program, and the two decimal figures the paper prints (`131^(9/10) = 80.45...` and
`1039^(9/10) = 518.74...`) are checked in the equivalent exact integer form `(r-1)^10 > r^9`. The
program prints one line per check and a closing verdict, and exits 0 only if every check passes.
The recorded run reports **40 checks, all passing**:

    VERDICT: ALL 40 CHECKS PASS

The run takes about one second on a laptop and needs no special hardware.

The design point worth naming is that the program does **not** trust the derived criterion the
paper's proof uses. `gamma_vector()` implements Definition 4.7.1 literally, building `C_g` as a
full length-`r` coefficient vector over the basis `lambda^0, ..., lambda^(r-1)` and testing it for
being identically zero; the exponent-congruence criterion `g^i = g^j mod r` is implemented
separately in `bad_set_exp()`; and the check `two-routes-agree` compares the two bad sets at every
one of the nine (prime, characteristic) pairs where both were run. The width-2 diagonal ROABP of
the paper's Proposition 2 is likewise not asserted but *executed*: `roabp_width2_diagonal()`
multiplies out the explicit diagonal matrices symbolically over `Z[x_1..x_n]` and the check
`roabp-width-2-diagonal-realisation` confirms the product is exactly `x_i - x_j` with every matrix
entry of degree at most 1 in its own variable.

Four controls guard against a decider tuned to manufacture large bad sets:

- `control-forced-positive-x1-plus-x2` --- `C = x_1 + x_2` at `r = 1039` must have `B = {1}`, since
  `g = g^2 mod r` only at `g = 1`; the engine returned `B = {1}`.
- `control-proved-silent-x1` --- `C = x_1` maps to a single basis vector and can never vanish; the
  engine returned `|B| = 0` in characteristics 0, 2 and 3.
- `control-characteristic-sensitivity` --- `C = x_1 + x_2 + x_3` at `r = 31` sends all three
  monomials to `lambda^1` at `g = 1` with coefficient sum 3, so `g = 1` must be bad over `F_3` and
  good over `F_2` and `Q`; the engine returned exactly that, so it is not blind to characteristic.
- `control-forced-negative-difference-1` --- `C = x_1 - x_2` has index difference 1, so the lemma
  forces `|B| = gcd(1, 30) = 1`; the engine returned a singleton.

A fifth check, `reading-is-nested-not-product`, records that the whole result depends on one
character of the source: under the alternative reading `x_i -> lambda^(g*i)` the bad set of
`x_1 - x_r` is **empty** at every prime tested. The nested exponent of Definition 4.7.1 is
load-bearing, and the paper says so in Section 1.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    d5fd3bc8150b2bcd155241ca7585cea815f84c7bde326641f37eb5f12944cdf8

The run recorded there was made on Python 3.9.6.

## Scope

**Does the program re-derive the paper's claim, or only confirm the exhibited object? Both, at
finite parameters --- and neither establishes Theorem 1 in general.** Taken check by check:

- It **confirms the exhibited object**: `C = x_1 - x_r` has `B_r(C) = Z_r^*` and `|B| = r - 1` at
  `r = 31, 131, 1039` by the full length-`r` ring engine in three characteristics, and at
  `r = 31, 131, 1039, 4099, 10007` by the independent congruence route.
- It **re-derives the lemma the paper's proof turns on**, rather than assuming it:
  `|B_r(x_i - x_j)| = gcd(j - i, r - 1)` is checked at *every* index difference `1..r-1` at
  `r = 31, 131, 1039`, the bad set is checked to be a subgroup at every index difference at
  `r = 131`, and `ring-basis-independence` establishes over all 465 exponent pairs at `r = 31`, in
  each of the characteristics 0, 2 and 3, that for a two-term `C` the condition `C_g = 0` really is
  an exponent congruence. It also *executes* rather than asserts the class membership of
  Proposition 2, multiplying out the explicit diagonal matrices symbolically. The exact-integer
  inequalities `130^10 > 131^9`, the least-`r` result `r >= 7`, and the halved witness's least prime
  `1039` are likewise recomputed, not quoted.
- It does **not** establish Theorem 1 for all `r`, and cannot: the theorem is quantified over every
  prime and every field, and its proof is the two-line Fermat/Euler argument in Section 1 of the
  paper. A referee who distrusts the program still has the theorem; a referee who distrusts the
  theorem does not get it from the program. The program is corroboration, and the paper says so.

The program's own closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the full length-r ring engine was run only at r = 31, 131 and 1039. The two largest
> primes the paper prints, r = 4099 and r = 10007, are checked by the exponent-congruence route
> alone; the two routes are shown to agree at the three smaller primes, but no coefficient vector
> of length 4099 or 10007 is built here.

> NOT RE-RUN: nothing here tests any ROABP of width >= 3, and nothing here tests any object of
> non-degenerate coefficient support density. The regime the source DESCRIBES in its p.27 Remark
> is therefore untouched by this program as it is by the paper, and remains open. This program
> exhibits no counterexample there and refutes nothing there.

> NOT RE-RUN: no search over a family of ROABPs is performed, and in particular the small-parameter
> census at r = 31, n = 5 that appears in this result's internal record is NOT reproduced here and
> is NOT offered as evidence anywhere. A sweep at n << r is structurally blind to the family of
> Theorem 1, so it could not bear on the surviving regime in either direction; Section 3 of the
> paper says so.

> NOT RE-RUN: the source PDF (arXiv:2602.13449v1, 1,086,923 bytes, PDF-only) is not fetched or
> parsed. Conjecture 4.7.2, Proposition 4.7.3, the p.27 Remark, Example 4.7.3, Definition 4.7.1,
> Algorithm 1 and the Appendix B table are TRANSCRIPTIONS from the printed pages named in the
> paper; this program cannot re-check a transcription, and a reader who doubts one should compare
> it against the PDF by eye.

> NOT RE-RUN: no prior-art search is performed here. The paper states that the witness OBJECT is
> the source's own Example 4.7.3 at S = {1}, S' = {r} and claims novelty only for the reversed
> quantifier order and the count |B_r(C)| = r - 1; this program checks the count, not the
> literature.

Three further limits belong to the paper rather than to the program, and it states all three in
Section 3.

The first is the one a referee should weigh first. The object `x_1 - x_r` is a difference of two
monomials of width 2 and support size 2, and the source's own p. 27 Remark removes objects of
"very small effective width (constant or logarithmic)" and "extremely low algebraic density" from
the intended scope of Conjecture 4.7.2 --- naming "differences of two monomials" as its own
example. So this is a **corpus correction** rather than the refutation of a live open problem: what
dies is the printed quantifier form of the conjecture and, decisively, the numbered unconditional
Proposition 4.7.3 on p. 28, which asserts the conjecture *precisely on width 2* and carries no
degeneracy hedge. The sharpest statement of the finding needs no reading-dependent judgement at
all: p. 27 excludes the class and p. 28 asserts the conjecture on it, so the source contradicts
itself. Section 3 of the paper quotes, in full, all **three** hedges the source prints between the
conjecture and Proposition 4.7.3 --- the p. 27 Remark, the separately headed p. 27 **Clarification**,
and a second, shorter Remark on p. 28 --- and answers each. Two are quantified `g`-first
("adversarial constructions that tailor `C` to a specific choice of `g` do not contradict the
conjecture"; "for a fixed `g`, there exist tailored nonzero polynomials `C` that are annihilated by
the map") and so do not reach an object fixed before any `g` that annihilates all `r-1` parameters
at once. The third is the strongest, and a referee will raise it first: the Clarification says the
conjecture "nor does it preclude the existence of isolated structured counterexamples". That does
not rescue either printed statement. It cannot be read as a hypothesis of Conjecture 4.7.2, whose
printed inequality is quantified over every nonzero width-`w`, degree-`d` ROABP with no exception
clause attached; and it cannot touch Proposition 4.7.3 at all, which asserts *unconditionally* that
the conjecture holds for width-2 diagonal ROABPs, so one counterexample in that class falsifies it
however isolated or structured it is. The Clarification does show that the authors foresaw objects
of this shape --- it exhibits none, bounds nothing, and withdraws neither printed statement. The
conjecture as
described remains open; this folder offers no repair of Proposition 4.7.3, and --- deliberately ---
offers no evidence about the surviving regime either, since a small-parameter sweep at `n << r` is
structurally blind to the failing family and would only look like evidence.

The second is attribution. `x_1 - x_r` is an instance of the source's own Example 4.7.3, so no
novelty is claimed for the object. The new content is the count `|B_r(C)| = r - 1`, which is
maximal and is exactly the quantity Conjecture 4.7.2 bounds; that Kronecker-style substitutions
collide once the exponent range meets the modulus is folklore, and the paper records that it has
**no fetched locator** for a prior statement of the count for a fixed `C` at a fixed `r`. A referee
who supplies such a locator reduces Theorem 1 to a restatement. Two literature channels were also
holed at search time and are recorded as open risk: keyword search on OpenAlex and on Semantic
Scholar both returned HTTP 429 (their by-identifier record endpoints did answer, with zero
citations), MathSciNet was not accessible, and Google Scholar was not consulted.

The third is the algorithmic reading. The paper's Remark on Algorithm 1 (printed p. 40) says that
the procedure needs only one good parameter and that `C = x_1 - x_r` leaves none, so on that input
it reports a nonzero polynomial as zero. That is a statement about the printed algorithm on this
input, taken from the printed text; no claim is made about whether the algorithm can be repaired,
and Theorem 7.5 is in any case stated conditionally on Conjecture 4.7.2.
