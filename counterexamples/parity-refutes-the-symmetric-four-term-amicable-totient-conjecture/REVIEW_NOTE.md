# Parity Refutes the Printed Four-Term Amicable-Totient Equations

`parity-refutes-the-symmetric-four-term-amicable-totient-conjecture`

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
passes. The recorded run reports **66 checks, all passing**:

    VERDICT: ALL 66 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    aa92eeb641d80ec17404736f38f6bc784cbeea44991f4d82e10c3ac51566c89d

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE
> NOT RE-RUN: the 415,523-term census of OEIS A002025. It needs the b-file, which is an external data file this program is forbidden to read, so the figures 14,484 cell members, 8,342,784 labellings, the n-distribution and the 1,073 type-(4,2) members are checked here only for internal arithmetic consistency (section H). No statement of the paper rests on them: Theorem 2 is a parity argument and Theorem 3 is one explicit pair.
> NOT RE-RUN: the two minimality claims that do rest on the census -- that A = 32642324 is the SMALLEST type-(4,4) amicable pair with gcd a power of two, and that no member of the cell has n = 1 or n > 8 below 10^17. Section H checks only that the witness is smallest among the eight pairs printed in the paper.
> NOT RE-RUN: the provenance of the target. The wording of the conjecture, the line numbers 576-582, 565, 222-230, 567-573, 585, 21, 33, the file size 25,122 B and the 664 lines of Amicable27122025.tex were established by fetching and reading the arXiv:2512.22319v1 e-print source; nothing here re-fetches it. This program checks mathematics, not bibliography.
> NOT RE-RUN: the three theorems of the target paper are checked here only for PARITY SHAPE and on the one (2,2) example the paper prints (section F). Their proofs are not reproduced, and the paper's (3,3) example is not recomputed.
> NOT RE-RUN: the offset-free equations for type-(4,4) pairs OTHER than the eight tabulated here. Parity does not reach them and no general argument is offered; they are open, as the paper says in its Status list.
> NOT RE-RUN: the prior-art question. Whether Proposition 4 is already in print was not settled -- the full texts of Borho-Hoffmann, Costello and Garcia-Pedersen-te Riele were not obtained, and MathSciNet was never consulted.
> NOT COVERED: the assertion at line 585 of the target, that x is in general a function of the (k-1)-th and (k-2)-th elementary symmetric polynomials of its OWN prime factors. Nothing here refutes it.
