# The Iterated Beatty Range Criterion Fails Below 5/4 and Holds Above 2

`the-iterated-beatty-range-criterion-fails-below-5-4-and-holds-above-2`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

The paper answers Question 1 of Khani, Valizadeh and Zarei (arXiv:2607.12817), which asks
whether their range criterion for iterated Beatty sequences extends from irrational
`alpha > (3+sqrt5)/2 = 2.618...` down to `alpha` in `(1, 2.618...)`. The answer is mixed:
**no** for every irrational `alpha` in `(1, 5/4]`, where an explicit fixed point of
`f_alpha` is a counterexample already at level `n = 2`; **yes** for every irrational
`alpha > 2` and every `n`, and at level `n` as soon as `alpha >= a_n`, the `n`-bonacci
constant. Nothing here contradicts the theorem of arXiv:2607.12817, whose hypothesis is
untouched; what is settled is the question.

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

Python 3.9 or later, standard library only: no third-party package and no external data
file. It runs in about 45 seconds on one core. The program prints one line per check and a
closing verdict, and exits 0 only if every check passes. The recorded run reports
**31 checks, all passing**:

    VERDICT: ALL 31 CHECKS PASS

It reads the objects exhibited in the paper --- the two moduli as (minimal polynomial,
isolating interval), the membership chain `3 -> 4 -> 6 -> 9 -> 14 -> 22`, the window
endpoints `(1055+896*sqrt2)/223` and `8+8*sqrt2`, the closed forms `13phi-21`, `8phi-13`,
`5phi-8`, `3phi-4` and `(9x+38-8(x-1)sqrt2)/47` --- and re-derives every quantity the
paper claims from them. There is no floating-point arithmetic and no floating-point
comparison anywhere: each modulus is an exact element `a + b*sqrt(d)` of a real quadratic
field with `a, b` rational, every inequality between algebraic numbers is decided by
integer arithmetic, and every decimal expansion printed in the paper (23 of them) is
checked by exact rational bracketing of its truncation rather than by re-printing a float.
Membership in the range of an iterate is settled by a forward sieve that is exact and
complete on `[1, X]`, and is checked against the printed chain independently of it.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by
an exit status, both written by the run harness. The header records the SHA-256 of the
program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    31adfec127c6eae0c7bf65d9f67ac4a01e0dab9982c1f0d9db2472dd6b1347ad

## Scope

The program's own statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the following are recorded in the paper but are NOT re-derived here.
>
> (a) the machine census reporting 38,417 disagreements at alpha=(7+sqrt2)/8 for x <= 10^5.
> No claim of the paper depends on it; this program checks x <= 24 exhaustively and nothing
> beyond.
>
> (b) any statement about which inequalities the proof printed in arXiv:2607.12817 uses
> (the threshold remark of Section 3). That is a reading of a text, not a computation; the
> line numbers, byte count and sha256 given in the paper let a referee repeat it, and this
> program does not fetch the source.
>
> (c) the literature search behind the novelty remark of Section 5 and behind Section 6.
> No database is contacted.
>
> (d) the open cells are not resolved: whether the equivalence fails anywhere on
> [phi,2), and whether it fails at alpha=phi for n=3, are left open. The searches here are
> BOUNDED (x <= 20000 at n=3, x <= 3000 otherwise) and a bounded negative is not a proof.
>
> (e) the family sweep behind Theorem 7 (the negative half) enumerates quadratic
> irrationals only, of the shape (p+q sqrt d)/r with d <= 43, q <= 3, r <= 25, p < 60;
> the theorem itself is proved for all irrational alpha in (1,phi) and its proof, not this
> sweep, is what carries it.

Two further limits are worth naming explicitly, since they bound what the paper claims
rather than what the program checks.

* **Neither theorem needs the program.** Theorem 6, Theorem 7 and Corollary 8 are proved
  in the paper by hand, and every witness is exhibited there in closed form, so a referee
  who ignores `verify.py` entirely loses nothing but corroboration. The program exists to
  catch transcription error, and it did: an earlier draft of this result carried the window
  as `(10.41352, 19.313626)`; the exact endpoints are `(1055+896*sqrt2)/223 =
  10.4131630129...` and `8+8*sqrt2 = 19.3137084989...`, both wrong in the fifth
  significant digit in that draft, with the integer window `11..19` unaffected either way.
* **Prior art caps the novelty of one witness, not of the result.** The golden-ratio
  witness of Section 5 is one step from the identity `f_phi^2 = f_phi + id - 1`, which is
  in print in Connell 1959, in Khani--Zarei 2022 and as an axiom of the target authors'
  own earlier preprint arXiv:2508.02303 --- the first two of which arXiv:2607.12817 itself
  cites. The paper says so in place. The `(1, 5/4]` family of Section 4 (Theorem 7 and
  Corollary 8) carries no Wythoff or Fibonacci structure and is free of that caveat. Separately, the `n = 1` clause of the
  criterion, which arXiv:2607.12817 calls folklore without a reference, is published
  (Fraenkel--Levitt--Shimshoni 1972; Szuesz 1985); the Szuesz citation is made from its
  Zentralblatt review, not from its three printed pages, and the paper records that.
