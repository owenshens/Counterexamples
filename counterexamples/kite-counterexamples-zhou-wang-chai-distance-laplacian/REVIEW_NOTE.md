# Kite Families with (n^3) Failure of the Zhou–Wang–Chai Distance-Laplacian Bound

`kite-counterexamples-zhou-wang-chai-distance-laplacian`

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
passes. The recorded run reports **22 checks, all passing**:

    VERDICT: ALL 22 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    46726f496ca934a99f6cb11f89eb3a4688de40b171ef3528dc975f865d774b88

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the theorem is an infinite family; the spectral verification above covers H_n for n = 5..40 and the cubic index r = 2m for m = 4..12 only.
> NOT RE-RUN: the remark's path statements are checked on finite ranges too: the excess of P_n at r = n-2 spectrally for n = 4..40, the selected-sum excess of P_{3m} in closed form for m = 2..40, and the true U_{2m}(P_{3m}) against the bound only for m = 3..6.
> NOT RE-RUN: the inequality under attack is transcribed from its published source and is never collated with it, because this program has no network access; every check tests U_r(G) <= W(G) + C(r+2,3) for 2 <= r <= n exactly as the paper quotes it, so a hypothesis dropped in transcription would not be detected here.  The only offline corroboration is the check diameter_two_violators_are_exactly_the_two_stars, which reproduces for this same inequality the diameter-two classification quoted from earlier work.
> NOT RE-RUN: the closing algebra is verified as a polynomial identity (hence for all n and m), but the two variational principles used by the proof, Rayleigh-Ritz and Ky Fan, are classical inputs and are exercised numerically, not proved here.
> NOT RE-RUN: the exhaustive census covers connected graphs of order at most 6; the diameter-two classification quoted from earlier work is confirmed only within that range.
