# A 71-Point Arc Refutes the Printed Finite-$p$ Form of the Entropic Cauchy--Davenport Candidate of Gavalakis, Goh and Kontoyiannis

`a-71-point-arc-refutes-the-printed-finite-p-form-of-the-gavalakis-goh-kontoyiannis-candidate`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

The source paper (Gavalakis, Goh and Kontoyiannis, arXiv:2604.20233v3) proposes, as a candidate
entropic Cauchy--Davenport theorem, the numbered display it labels `eqiidCD`:

 H(X+X') >= min{ log p, log( (2^{H(X)+1} - 1) / sqrt(2) ) }

for every prime `p`, every `F_p`-valued `X`, and `X'` an independent copy of `X`, with `log` to
base 2. The claim of this note is that the first line of that display, read literally at finite
`p`, is **false**, and false in the branch where the minimum is *not* the cap `log p`. The witness is `p = 101` with `X` uniform on the arc
`S = {0,1,...,70}`.

Two things about that scope are load-bearing and are stated in the paper rather than buried:

* **`eqiidCD` is a hedged, author-labelled candidate display, and the paper does not call it a
 conjecture.** The authors introduce it with "a more appropriate candidate ... might be that".
 It is nevertheless numbered, back-referenced by number by the authors, and a definite
 assertion, which is why it can be refuted. The paper says all of this explicitly in Section 1.
* **The failure is in the *unsaturated* branch.** A failure at `min = log p` would be a
 degenerate artefact: there the display degenerates into the demand that `X+X'` be exactly
 uniform. That is *not* what this note reports; the paper places the failure strictly below the
 saturation threshold, in the branch that carries the content `H(X) + 1/2 - o(1)`, precisely so
 a reader does not mistake the degenerate case for the content of Theorem 1.

## What was checked, and how

The decisive step is an inequality between two explicit positive integers, so the whole check is
integer arithmetic and no logarithm is ever evaluated. Two facts do the work:

1. For `X` uniform on a set of size `m`, `2^{H(X)+1} - 1 = 2m - 1` **exactly**, so which branch
 of the `min` is active is the integer comparison `(2m-1)^2` against `2p^2`. At `(p,m) =
 (101,71)`: `19881 < 20402`, so the active branch is `log(141/sqrt(2))`, strictly below
 `log 101`.
2. With `C_k = #{(i,j) in S x S : i+j = k mod p}` and `P = prod_k C_k^{C_k}`,
 `H(X+X') < log((2m-1)/sqrt(2))` is *equivalent* to `(2m^4)^{m^2} < (2m-1)^{2m^2} P^2`
 (Lemma 2 of the paper, two lines). At `(101,71)` those integers have 38847 and 38874 decimal
 digits and the inequality holds.

`verify.py` re-derives all of it from the two numbers printed in the paper --- the prime 101 and
the arc `{0,...,70}` --- and in particular:

* computes the convolution counts twice, once by an `O(m^2)` brute-force double loop and once by
 the closed form `C_k = f(k) + f(k+p)` of equation (5), and checks the two agree on all 101
 residues;
* checks the printed multiset `41^(42) 42^(2) ... 70^(2) 71^(1)` and its sum `5041 = 71^2`;
* checks the branch test, and that the saturation floor at `p = 101` is `m = 72`, one step above
 the witness;
* checks `L < R` and both digit counts, and brackets the deficit **in integers** as
 `2^91 L < R < 2^92 L`, i.e. `91/10082 < deficit < 92/10082` bits --- so the violation is
 bounded away from zero without evaluating a logarithm;
* re-runs the same two tests on the six further witnesses of the table in Section 4 and on the
 complete per-prime out-of-band census of Section 4, and brackets **every** witness's deficit
 in integers (the `k` column of that table, deficit in `(k/2m^2, (k+1)/2m^2)`);
* certifies, again by integer cross-multiplication, that those seven deficits are **not
 monotone** in `p` and that the deficit at `p = 401` is more than 35x smaller than at
 `p = 151` --- the arithmetic behind the Section 5 statement that the census is no evidence
 either way about the asymptotic form;
* certifies the two scope statements of Section 5 (`0.50 < log p - H(X) < 0.53`, and
 `H(X+X') - H(X) > 0.477`) by integer comparisons of the form `x^b > 2^a`.

**Both control polarities are exercised, and one of them is a disclosed misfire.** Seven
out-of-band arcs that do *not* violate the display are checked to return "no violation", so the
decider can say NO. And the criterion of Lemma 2 is *branch-conditional*: applied in the
saturated branch it reports a violation where there is none, because it compares against the
inactive branch. `verify.py` demonstrates that misfire at `p = 7, 11, 13, 61, 101` with
`S = F_p`, checks that the branch test classifies every one of those as in band --- so the guard
is what stops the criterion being applied there --- and checks that all seven reported witnesses
are out of band, the regime where the criterion is valid. The paper states this as a Remark in
Section 3. It is disclosed rather than suppressed because in the row this folder comes from it
was originally a **failed control**.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
Runs in about 4 seconds. Every decision is an integer comparison; the only floating-point
numbers printed are `NOTE` lines labelled "orientation only" or "corroboration (NOT a
decision)", and no `PASS` depends on one. The program prints one line per check and a closing
verdict, and exits 0 only if every check passes. The recorded run reports **51 checks, all
passing**:

 VERDICT: ALL 51 CHECKS PASS

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 b708e0b3f14bdbd42b2e401e1c93e2c8f650bd506257fdd8177c4a4c9902e39b

`verify.py` is **new to this folder**: it was written for this hand-over, on the objects printed
in the paper, and it was run locally under CPython 3.9.6 on the machine that assembled the
folder. It is not the program that originally found the witnesses.

The original computation is recorded in this wave's artifacts manifest, at
the run record's `MANIFEST.json` in the producing repository. **That manifest and the
scratch files it indexes are not shipped inside this folder**; the four files in this folder are
the complete contents. Nothing in the paper depends on the manifest, and everything about it that
bears on the claims is stated below rather than referred out. Its record has limits which must not
be laundered here:

* The witnesses were found by `arcs3.py` (a complete arc census, `m = 1..p`, over every prime
 `p <= 509`) and re-decided independently by `intcert4.py` (pure-integer certificates plus the
 controls). Both ran as detached jobs on the fleet, `arcs3.py` on one slot
 (the recorded dispatch id `(a dispatch id, redacted)`) and `intcert4.py` on
 a second slot (the recorded dispatch id `(a dispatch id, redacted)`), both fleet-runner Status
 Success with job `RC=0`. A first attempt at `intcert4.py` (the recorded dispatch id
 `(a dispatch id, redacted)`) exited `RC=1` on CPython's `int`-to-`str` 4300-digit
 limit and was read as INCOMPLETE, not as a negative.
* **No digest was filed at dispatch time.** The dispatch harness printed
 `ARTIFACT_NOT_FILED=no wave` for all four jobs, so the files in the run record are the
 attack agent's scratch copies from a scratch directory, copied byte-for-byte and digested where
 they now sit. The invocation lines in that manifest are the runnable form against the filed
 copy, not a transcript of the line that actually ran.
* **Two of the four recorded outputs are truncated and there is no archived copy.**
 `arcs3.out` and `census.out` are each exactly 20045 bytes and begin mid-line, because the
 dispatch harness prints only the last 20000 bytes of a job log; the planned object-store keys were never
 written. A referee who wants the full arc-census log must re-run rather than fetch. Nothing in
 the present paper depends on the missing head of those logs: every number the paper prints is
 either in the intact `intcert4.out` or in the intact tail of `arcs3.out`, and all of them are
 recomputed from scratch by `verify.py`.
* **Wall-clock timings were not captured by any stage** of the original run and are therefore
 not stated anywhere.

**One correction to the original record, found by this folder's re-run.** `intcert4.py`
estimated decimal digit counts as `1 + int(bit_length * log10 2)` without correcting the
estimate, and that is off by one for the three largest cases. The exact counts for the right-hand
integer are `436214` at `p = 307`, `778821` at `p = 401` and `1307509` at `p = 509`, not the
`436215`, `778822`, `1307510` recorded in the run record's `certificate.txt`. The paper prints
the exact values, `verify.py` computes them by integer comparison against powers of ten, and
`verify.output.txt` states the discrepancy in a `NOTE`. No verdict anywhere depends on a digit
count.

**One correction to the paper itself, found by the bundle review.** An earlier draft of Section 5
asserted that "the deficit does not shrink over `61 <= p <= 509`". That is **false**: the deficit
at `p = 401` is more than 35 times smaller than at `p = 151`. It was the single sentence in the
paper that editorialised past its own integers, and it was doing so in the direction that
flattered the result. Section 5 now states the seven deficits, states that they are non-monotone
and shrinking, and claims nothing about `p -> infinity`; the table in Section 4 carries the
integer bracket `k` for each witness so a referee sees the trend directly; and `verify.py`
certifies both the brackets and the shrinkage by integer cross-multiplication (Step 7b). No other
claim in the paper moved, and Theorem 1 is untouched --- it refutes the first line of
`eqiidCD`, read literally at finite `p`, that line carrying no `o(1)`.

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN / SCOPE. This program re-derives the claims of the paper and nothing more. It does
> NOT search primes p > 509; it does NOT examine any set other than arcs {0,...,m-1} (no
> two-arc, no general subset, no non-uniform X); it does NOT establish anything about the
> o_{p->infinity}(1) term of the candidate, which the paper explicitly leaves open; it does NOT
> touch the bulk regime where 2^{H(X)} is far below p/sqrt2; it does NOT re-derive the
> saturated-branch Fourier remark of Section 5, which the paper states as folklore and does not
> claim; and it does NOT verify the constant K_epsilon of the source paper's windowed theorem,
> which is used only qualitatively in Section 5.

Beyond that, the following are open or unchecked and are stated as such in Section 5 of the
paper:

* **The asymptotic form is not refuted, and the census is no evidence either way.** Only the
 first line of `eqiidCD` is refuted; the second carries an `o_{p->infinity}(1)` and is
 untouched. The deficits at the seven witnesses are *not* monotone in `p` and they do **shrink**
 substantially over the census --- `5.7e-3, 9.0e-3, 1.07e-2, 3.70e-3, 4.1e-4, 3.0e-4, 1.27e-3`
 at `p = 61,101,151,211,307,401,509`, so `p = 401` is more than 35x below `p = 151`. On this
 data the deficit of the smallest out-of-band arc may well tend to 0, in which case the
 asymptotic line survives intact; the paper says so and conjectures nothing. (Pulling the other
 way, and equally unresolved: each row is the *smallest* out-of-band arc at its prime, not the
 worst, so these are not maximal deficits; and only arcs are censused.) Note also that
 out-of-band failure is *intermittent* at small primes: there is no out-of-band violating arc
 at all at `p = 31, 37, 43, 47, 53`.
* **No theorem of the source paper, or of the literature, is contradicted.** Every witness has
 `log p - H(X)` between 0.50 and 0.53 bits, i.e. sits at the near-uniform end of the entropy
 range. The paper does not determine the constant `K_epsilon` of the window of the authors' own
 Theorem `thmFpEPI` and claims nothing about its size, so it does not assert that the window
 excludes these witnesses; what is checked directly is that all seven witnesses satisfy that
 theorem for every `eps > 0.023`.
* **The bulk regime is untouched**, and it is not addressed by the paper at all.
* **No minimality claim.** Nothing here asserts that `p = 101` or `m = 71` is extremal.
* **Prior-art channels that were not readable.** The prior-art search behind this row reached
 the arXiv API, Semantic Scholar (both the arXiv and the DOI records, two citers, both read in
 full as LaTeX) and Crossref. **OpenAlex, zbMATH and MathSciNet were NOT read** --- OpenAlex
 returned HTTP 429 and then timed out, zbMATH's API returned HTTP 404 and its front end 403,
 and MathSciNet is paywalled here. So a prior observation recorded only in a review database,
 or in correspondence with the authors, would be invisible to us. The strongest evidence
 against that is that version 3 of the source (28 August 2026) was a revision made to patch
 errors pointed out by a reader and it left `eqiidCD` untouched.
