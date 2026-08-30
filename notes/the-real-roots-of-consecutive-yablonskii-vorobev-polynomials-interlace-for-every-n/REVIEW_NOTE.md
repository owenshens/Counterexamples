# Gap-Only Interlacing of the Real Roots of Consecutive Yablonskii--Vorob'ev Polynomials, Given the Published Real-Root Count

`the-real-roots-of-consecutive-yablonskii-vorobev-polynomials-interlace-for-every-n`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run. The theorem itself is proved by hand in the paper and needs
no machine; the program is a control on the interior of that proof.

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
Exact integers and exact rationals throughout -- there is no floating-point number anywhere in
the file and no tolerance of any kind; every decision is an integer sign or an integer
equality. It takes about fifteen seconds on a laptop. The program prints one line per check
and a closing verdict, and exits 0 only if every check passes. The recorded run reports
**24 checks, all passing**:

    VERDICT: ALL 24 CHECKS PASS

It reads the seed `Y_0 = 1`, `Y_1 = z` and the recurrence (1) of the paper, and derives from
them, and from nothing else, every quantity it compares against the paper's statements. It
generates `Y_0, ..., Y_43` by exact division in `Z[z]` rather than reading a coefficient table
in; it isolates real roots exactly with a Fujiwara bound and rational bisection, where the
count `floor((n+1)/2)` expected from (F3) is what certifies the isolation complete, and an
independent Sturm computation that uses (F3) nowhere confirms that count only for `n <= 12`;
and it re-derives the occupancy statements (3) and their forced counterparts for `Y_{n+1}` for
`1 <= n <= 40`. Three checks in the recorded run (`regression-coefficients`,
`regression-intervals`, `regression-words`) compare the generated objects against coefficient
lists, isolating brackets and merged root words held as literals inside `verify.py`. Those
literals are not displayed in the paper, and the checks are labelled in the run as the internal
regression checks on the generator, the isolator and the merge that they are; no check asserts
agreement with a display, a table or a numbered result the paper does not contain.

Two of the 24 checks are deliberately negative, because a control that can only pass is not a
control. The **anti-control** exhibits `M = (x-1)(x-2)`, `A = 2x-3` and
`B = (5x-6)(5x-7)(5x-8)`: `A` and `B` satisfy both the hypothesis and the conclusion of the
paper's Lemma 4 with respect to `M`, yet `(M,A)` interlace and `(M,B)` do not -- so the parity
lemma alone is strictly insufficient and the published root count (F3) is load-bearing rather
than decorative. The **forced-negative control** perturbs the constant term of `Y_3` and
confirms that the same isolation routine which returns the passes above raises an error instead
of returning when the count it is asked for is unattainable.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

## Scope

The program's own closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the published inputs are used as inputs. The program re-derives (F1), (F2) and
> the COUNT asserted by (F3) for n <= 41 by exact algebra, but it does not verify the PROOFS of
> Fukutani-Okamoto-Umemura (Nagoya Math. J. 159 (2000) 179-200), of Clarkson (Semin. Congr. 14
> (2006) 21-52) or of Roffelsen (SIGMA 8 (2012) 099). The theorem of the paper is conditional
> on those, exactly as it says.

> NOT RE-RUN: the theorem holds for EVERY n and is proved by hand in the paper. The root-order
> data above is a control on the interior of that induction for 1 <= n <= 40 only. n > 40 is
> covered by the proof and by nothing in this file.

> NOT RE-RUN: interlacing means throughout the GAP-ONLY definition quoted in section 1 of the
> paper. Nothing in this file establishes the statement under the stronger variant that also
> demands |m_{n+1} - m_n| = 1; the stronger-definition-variant check above records instead that
> under that variant the statement FAILS for every odd n.

> NOT RE-RUN: no claim of minimality, extremality or uniqueness is tested; no polynomial family
> other than the Yablonskii-Vorob'ev family is examined; and the generalised Okamoto
> polynomials that are the subject of the source paper are not touched.

> NOT RE-RUN: the bibliographic material of section 1 -- the quoted sentence, the definition of
> interlacing and the numbered results attributed to arXiv:2402.15887v1 and to the other sources
> -- is not machine-checked here. It was read off those sources by hand, and the journal version
> of record was not reachable, so no preprint-to-journal diff of the quoted sentence was
> performed.

The first of those limits is the one that shapes the headline: the result is stated under the
definition of interlacing used by Roffelsen and Stokes, a two-sided condition on gaps, and the
title, the abstract and the theorem all say so. Under a stronger variant that additionally
demands that the two root counts differ by one, the same proved word makes the statement true for
every even `n` and **false for every odd `n`**. One further limit is recorded in the paper rather
than in the program, because it is not computational: the index-2 interlacing theorem of Clarkson
for the pair `Y_{n-1}, Y_{n+1}` is neither implied by nor implies the result here, and no share of
it is claimed.
