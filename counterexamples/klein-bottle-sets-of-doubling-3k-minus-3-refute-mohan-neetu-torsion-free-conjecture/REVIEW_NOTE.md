# Klein Bottle Sets of Doubling 3|S|-3 Refute a Torsion-Free Small-Doubling Conjecture of Mohan and Neetu as Stated

`klein-bottle-sets-of-doubling-3k-minus-3-refute-mohan-neetu-torsion-free-conjecture`

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
passes. The recorded run reports **304 checks, all passing**:

    VERDICT: ALL 304 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    e46199f25989cd6723c886df4f61dcab1abd74613aa3e14dbb803b9c864073ba

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: torsion-freeness of K is confirmed structurally and over a bounded box (62 elements, exponents 1..12) only; the unbounded statement is Lemma 1, a hand proof.
> NOT RE-RUN: no search for a smaller or different counterexample, in K or in any other group; no minimality of 3|S|-3 is examined, and k > 30 is not enumerated.
> NOT RE-RUN: the byte locator of section 1 (source member sha256, tex lines 906-909, bytes 90,881-91,090) is not re-fetched -- this program has no network and ships no copy of the e-print.
> NOT RE-RUN: the BS(1,q) argument that (x g^i)^2 = x^2 forces q = -1 is algebra over Z[1/q], not a finite computation.
> NOT RE-RUN: no prior-art, attribution or bibliographic claim of the paper is verified.
