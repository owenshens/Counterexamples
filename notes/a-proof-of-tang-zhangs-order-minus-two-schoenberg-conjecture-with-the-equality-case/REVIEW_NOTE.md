# A Proof of Tang and Zhang's Order -2 Schoenberg Conjecture, with the Equality Case

`a-proof-of-tang-zhangs-order-minus-two-schoenberg-conjecture-with-the-equality-case`

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
passes. The recorded run reports **67 checks, all passing**:

    VERDICT: ALL 67 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    f5dc0380ca031ec8f51d4a096c3f0e5f6ddf780deaa62bdbd264f9d570fa8edd

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: Schur's inequality itself (Lemma 2 of the paper) is NOT proved here. It is quoted from Horn-Johnson and is also the source paper's own lem:weyl. Every check above is an identity or a finite instance; the passage from 'K is normal' to 'equality in eq:C1' is the LEMMA, not the program.
> NOT RE-RUN: the 'only if' half for n >= 3 is NOT machine-checked. Step 8 decides eq:C1 exactly, including strictness, only at n = 2, where the critical points are the roots of a quadratic; for n >= 3 the program verifies the normality criterion (Step 9) but not the strict inequality, because sum_k 1/|w_k|^2 is not a rational function of the z_j when the w_k are complex.
> NOT RE-RUN: no statement about the other orders. The paper's order -1, -4, -2m and -p inequalities, and the dual order -2 corollary of the source's later section, are outside both the paper and this program.
> NOT RE-RUN: nothing is verified about arbitrary real or complex z_j -- every instance above has Gaussian-rational or rational entries, chosen from a seeded generator over [-5,5]^2 and [-6,6]. The theorem is proved for all z_j by the argument of Section 3; the batteries corroborate the identities the argument uses.
> NOT RE-RUN: no search for a counterexample was performed and none could be, the statement being proved. The anti-controls (constant 3 -> 2, and the unsandwiched Frobenius norm) are the only negative-direction evidence, and both fire.
