# An Order-14 Counterexample to Conjecture 4 of Pioge et al.

`order-14-gram-matrix-pioge-conjecture-4`

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
passes. The recorded run reports **40 checks, all passing**:

    VERDICT: ALL 40 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    06483ebc3042d86926e8807bf2a93736559745a87154699d9dd42435fed159bf

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the original authors' unsuccessful numerical search in order 15; orders at most 13, about which the paper claims nothing; the padding corollary beyond k = 1, i.e. orders 16 and above are not each recomputed, only the order-15 instance of the same block identity; the explicit unitary completion U in C^{14x14} of the two orthonormal rows of Y = eta^{-1/2} X W^{1/2}, which is not constructed here and is not exactly representable, eta^{-1/2} being irrational -- what is checked instead is X W X^* = eta I_2, i.e. that those two rows are orthonormal, together with A W A = eta A and the trace-2 identity, so that P = eta^{-1} D A D is an exact Hermitian rank-two idempotent and the completion follows by Gram-Schmidt; and the transcendental Gaussian time-delay matrix S(d) itself, whose d^2 coefficient is instead obtained exactly from F_A by the two independent routes checked above.
