# Two-Circulant Butson Hadamard Matrices of Orders 46 and 74

`two-circulant-butson-hadamard-matrices-of-orders-46-and-74`

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
passes. The recorded run reports **63 checks, all passing**:

    VERDICT: ALL 63 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    8ac467d190f63b4b4bea549567f6a24c14c9659e96e011b2988e0bae3dfa8ca3

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the two lane censuses. This program does NOT re-derive the search programs' complete-count claims -- 10,120 ordered / 5,060 unordered complementary pairs over all 6^11 = 362,797,056 first rows per side at (p,k) = (23,11), and 5,904 ordered / 2,952 unordered over 6^12 = 2,176,782,336 at (37,12) -- nor the shard-cover accounting of the 128-thread p = 37 run, nor the searchers' own forced-positive controls at p = 3,5,7,11,13,17 and p = 31 (k=6), nor the checksum of the C engine. Those rest on the search programs' counters; their standard output was not preserved, and reproducing them needs 456 s single-threaded at p = 23 and 903 s on 128 threads at p = 37. Nothing above depends on them: existence needs one object, and the objects are the printed exponent vectors.
> NOT RE-RUN: minimality, or any classification of BH(46,6)/BH(74,6). Nothing here says these are the smallest or the only such matrices. The only separation established between W1 and W2 is that they sit on different Loeschian branches, so no symmetry of the two-circulant ansatz carries one to the other; whether H(W1) and H(W2) are inequivalent as Butson Hadamard matrices, under the far larger group of row/column permutations and sixth-root scalings, is NOT tested. The k = 22 lane at p = 23 (6^22 vectors) and the k = 23, 46 lanes at p = 47 were not searched, and BH(86,6) at p = 43 was not searched at all.
> NOT RE-RUN: the literature search. This program checks mathematics, not novelty; the unread bibliographic channels are named in the paper.
