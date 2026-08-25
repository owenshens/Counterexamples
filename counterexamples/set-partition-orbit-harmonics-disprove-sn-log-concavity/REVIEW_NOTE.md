# Counterexamples to Zhu's S_n-Log-Concavity Conjecture

`set-partition-orbit-harmonics-disprove-sn-log-concavity`

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

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **32 checks, all passing**:

    VERDICT: ALL 32 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    6e7c509e04cf62217c1ccc54357f33509a79598d8a38a7c4ec9f1d93ea62f382

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE not re-run here: the isomorphism R(Pi_{n,m})_d = C[X_{n,m,d}] (Zhu, Prop. 3.29) is taken as given, so the orbit-harmonics quotient is never built; the multiplicities at n=20 and beyond come from the complement census, which is cross-validated against brute-force character theory only for 4 <= n <= 8.  The dependence on k is settled for all odd k (source: the surviving 4-edge classes are the same in every even ambient grid; target k >= 7: every 5-edge core uses at most 5 vertices per side, so the isolated-vertex veto always fires), with censuses run explicitly up to k = 61.
> NOTE scope: this program does not re-derive the paper's full claim.  What it re-proves by exhaustion is Lemma 3 of the paper (the four-edge / five-edge classification) together with the sign multiplicities and the family arithmetic that follow from it; the paper has four further load-bearing steps this program does not verify on its own terms.  (a) Lemma 2, the orbit criterion via Frobenius reciprocity, is assumed by the census rather than proved.  (b) So is the vanishing of NON-TRANSVERSE orbits (a block pair with |B and C| >= 2 admits an odd transposition in the stabiliser): the census enumerates transverse pairs only.  Steps (a) and (b) are the exact content of the census-vs-character-theory agreement, since direct_sign_mult counts every orbit straight from the definition -- so they are corroborated, but only for 4 <= n <= 8, not at n = 20.  (c) The closing inference -- complex S_n-representations are semisimple, so an injection would force source multiplicities to be dominated by target ones -- is pure representation theory and is not checked numerically; the checks establish only the strict drop 4 > 1 and 4 > 0.  (d) Over the m-range k+1 <= m <= n, every m is swept for k = 5,...,13, but for k up to 61 only m = k+1 is run; that is the binding case because the admissibility cap 'every degree of G is at most m' is monotone in m, so larger m can only admit more.  No claim is made for even k or for k = 3, neither of which the paper asserts.
