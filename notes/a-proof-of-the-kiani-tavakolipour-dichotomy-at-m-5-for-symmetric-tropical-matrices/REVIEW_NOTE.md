# A Proof of the Kiani–Tavakolipour Dichotomy at m = 5 for Symmetric Tropical Matrices

`a-proof-of-the-kiani-tavakolipour-dichotomy-at-m-5-for-symmetric-tropical-matrices`

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

Python 3.9 or later, standard library only: no third-party package, no external data file and no
network access. All arithmetic is exact — Python integers and `fractions.Fraction` — and
`-infinity` is carried as the sentinel `None`, never as a float, so no decision in the program is
taken on a floating-point comparison. The program prints one line per check and a closing verdict,
and exits 0 only if every check passes. The recorded run reports **47 checks, all passing**:

    VERDICT: ALL 47 CHECKS PASS

It reads the objects exhibited in the paper — the sample certificates and the census
specification — and re-derives every quantity the paper claims about them. It checks a number of
further objects that the paper, in its revised shorter form, no longer prints: a third sample
certificate at m = 4, a nonsymmetric 4×4 both-fail example and its all-finite variants, a symmetric C₅ witness that −∞ does not walk down the
cycle types, two single-branch control matrices C_I and C_II, and the δ-sequences of the source's own
4×4 and 9×9 examples, whose entries are transcribed into the program from the source e-print. In
particular it re-runs the m = 5 census itself
(378,000 triples, falling into 53,900 classes with identical constraint systems under
order-preserving relabelling of the index set), re-substitutes each of the 53,900 multiplier vectors
into its own forms, and re-runs the m = 3 and m = 4 censuses and the independent unpooled n = 5
census of 7,000 linear programs with the full 135-constraint system. It also carries anti-controls:
on each of the two symmetric 5×5 matrices that fail exactly one branch, the corresponding
single-branch system must come back **live**, and does, so the engine is not certifying everything
put to it; and the nonsymmetric m = 4 census must be non-silent, and is (132 live classes of 670).

The run takes about six minutes single-threaded; the m = 5 census is about two minutes of that and
the unpooled n = 5 control census about three.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    f453c0728e71858b94ba3d2d7883560f215ca8d745c20cd4588294ebda061563

## Scope

The program's own statement of what it does not cover, quoted from its output:

> NOT RE-RUN: m >= 6.  This program checks m = 5 (and re-proves the published m = 3 and m = 4
> cases); the conjecture ranges over 3 <= m <= n and m >= 6 is untouched here.  For m = 6 the same
> reduction gives the window V0 = [13] and a raw triple count of about 1.9 x 10^8, which this
> program makes no attempt at.
>
> NOT RE-RUN: the 378,000 individual multiplier vectors are not shipped as data.  They are
> regenerated and re-verified by this run, but no file in this folder lists them, so a referee
> wanting the full list must run this program.
>
> NOT RE-RUN: the converse of the reduction.  Only the forward direction is checked and only it is
> used -- a pooled certificate is a certificate at every n.  The claim that a LIVE pooled triple
> always lifts to a real counterexample is not verified here and no conclusion of the paper rests
> on it.
>
> NOT RE-RUN: no prior-art or literature search is performed by this program; it settles
> mathematics, not novelty.

One further limit belongs to the paper rather than to the program: the literature paragraph rests on
the arXiv API, zbMATH Open and the Semantic Scholar citation and reference graphs; MathSciNet was
not available and OpenAlex answered HTTP 429 on every attempt, so a prior observation recorded only
in those places would not have been seen.
