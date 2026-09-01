# A Primitive Weird Number with Two Square Odd Prime Factors and Omega = 8

`an-omega-8-primitive-weird-number-with-two-square-odd-prime-factors-settles-omega-2`

Supporting material for this paper: the program that re-derives its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

That

 m = 2^2 * 13^2 * 19^2 * 46219 * 1108619 = 12504224434300196

is a primitive weird number with exactly two square odd prime factors (13 and 19) and
`Omega(m) = 8`; and that, combined with the lower bound `Omega_2 >= 8` proved by Amato, Hasler,
Melfi and Parton, this gives `Omega_2 = 8` exactly, against the `8 <= Omega_2 <= 12` recorded in
Open Question 1 of their paper.

The decisive part needs no computer and no search. `sigma(m) - 2m = 8`, and the divisors of `m`
that are at most `8` are `1, 2, 4`, whose total is `7`; so `Delta(m) = 8` has no representation
whatever as a sum of distinct divisors of `m`, and `m` is weird by the criterion proved as
Lemma 2 of the paper. Primitivity is five exact integer comparisons, printed as a table in the
paper (Section 2, item (iv)). A referee who wants to check nothing else should check those two
lines.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **48 checks, all passing**:

 VERDICT: ALL 48 CHECKS PASS

Its only inputs are the strings printed in the paper -- the factorisation in the paper's own
`2^2 * 13^2 * ...` notation, the value of `m`, the product form of `sigma(m)`, and the five-row
table of maximal divisors, together with an index sequence for `m` that the paper does not print.
Everything else is re-derived from the
factorisation and compared against what the paper prints: `Omega` and `omega`, the two square odd
prime factors under both readings of "square odd prime factor", `sigma(m)` by three independent
routes (the multiplicative formula, the paper's printed product form, and explicit summation over
all 108 divisors), `Delta(m) = 8`, weirdness by exhaustive subset-sum over the eligible divisors,
and primitive abundance by checking **all 107 proper divisors** rather than only the five maximal
ones. It also re-derives that index sequence, `[1^2, 2^2, 1^2, 167, -1]`, from its own
sieve, and re-derives `Delta` a second time through the terminal-prime recursion
`sigma(n) - r(2n - sigma(n))`, a route that never forms `sigma(m)`.

Both polarities of the weirdness decider are exercised on published integers before it is
believed: `70`, `836`, `4030`, `5830` (OEIS `A002975`) are decided WEIRD, and `12`, `20`, `945`
are decided semiperfect, so the decider can say both yes and no.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 97326964b8647292e97035b4db95063839f2ac6d36cd1e282161f6045477211b

`verify.py` was written for this folder and run locally (Python 3.9.6, exit status 0). It is not
one of the programs that produced the result; it is an independent re-derivation of the paper's
printed quantities.

The result itself came out of two recorded jobs, whose provenance is transcribed here from the
run's artifacts manifest and not reconstructed:

- `run7.py`, sha256 `17c630a25d892f900d76e53bc27bb15b550e8bb4b32eb2cb1868f9933d600784`,
 27,454 B, standard library only -- the exhaustive census of the cell that found `m`. Dispatched
 as `aws/slot_run.sh --detach AUTO .../run7.py 3300 t8202-final`, a separate cloud instance, 24 worker processes,
 target phase 203 s, `RC=0`. No randomness and no seed; the workers shard a fixed prefix tree.
- `idx.py`, sha256 `af9c6c1d5a5e67f3b3b58cc9f18447c40046eb6a1b879b69191c579668681c07`, 1,383 B --
 the index-sequence job. a separate cloud instance, `RC=0`. It imports `sympy`, so it is *not* the program
 shipped here; `verify.py` redoes the same derivation with a standard-library sieve.

One limit of that record, which the manifest states and which is repeated here rather than
smoothed over: both jobs were dispatched from scratch paths outside the run directory, so the
harness printed `ARTIFACT_NOT_FILED=no wave` and never filed the code at dispatch time. The two
scripts were copied byte-for-byte into the run's artifacts directory afterwards, and the two
invocation lines above are the **reproducing form** rather than a transcribed
`ARTIFACT_INVOCATION=` line, which was never printed. Neither script is shipped in this folder,
and nothing in the paper depends on either of them: the paper's proof is items (i)--(v) of
Section 2, and `verify.py` re-derives every number in it from the factorisation alone.

## Scope

- **The case `n = 2` only.** Open Question 1 quantifies over every `n`. `Omega_1 = 7` was already
 known; `Omega_3` and `Omega_n` for `n >= 3` are untouched here.
- **The equality is conditional on a published theorem.** `Omega_2 <= 8` is what the witness
 proves. `Omega_2 >= 8` is `thm:patterns` of Amato-Hasler-Melfi-Parton, transcribed and *not*
 re-proved here; `verify.py` says so in its closing scope lines and treats the bound as an input.
- **The census is not claimed.** `m` was found by an exhaustive census of the cell, which also
 reports that cell fully decided with exactly one PWN in it -- which would give `Omega_3 >= 9`
 and exclude an exponent-`>= 3` companion at `Omega = 8`. Those by-products are deliberately
 *not* claimed by the paper and are *not* checked by `verify.py`; only the witness has been
 independently re-derived.
- **The integer is not new.** `12504224434300196` is the tenth term of OEIS `A063788`
 ("numbers `k` such that `sigma(k) = 2k + Omega(k)`"), deposited 27 May 2025 under that
 unrelated definition. The paper makes no priority claim for the observation that `m` is a PWN of
 this class, nor for the closure: it has not been established that either is unrecorded.
 `verify.py` checks the `A063788` relation as a consistency check on that third-party record.
- **The quotations are from the arXiv v2 source, not the journal galley.** All verbatim text in
 Section 1 comes from `many-factors.tex` inside the arXiv:1802.07178v2
 e-print. The published galley was obtained but was not compared clause by clause with v2, so a
 proof-stage change to the Open Question's wording, or to the `Omega` column of the table that
 supplies `Omega_2 <= 12`, is unchecked. Theorem 1 does not depend on either.
- **Prior-art search is not exhaustive.** arXiv, zbMATH, OEIS and OpenAlex work-record queries
 answered and are reported in the run record; two channels did not answer (Semantic Scholar
 keyword search and the OpenAlex citer query, HTTP 429 on independent retries), and MathSciNet
 was not consulted. Nothing scooping the witness was found, but those gaps are real.
- **Floating point decides nothing.** Every comparison in `verify.py` is over `int` or
 `fractions.Fraction`; the index-sequence ratio is handled in exact form (`488061/11`).
