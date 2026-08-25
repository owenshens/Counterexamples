# A Counterexample to the Weak Drury Permanent Inequality

`rank-two-counterexample-zhang-weak-drury-inequality`

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
passes. The recorded run reports **52 checks, all passing**:

    VERDICT: ALL 52 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    338e6456518e249d9a3b79c0b2c7b3342e7e9795741ab2b25b1615a6bfb829fa

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: what is refuted here is the inequality Phi_1(A) <= per(A) exactly as printed in this paper, with Phi_1 defined by the paper's own displayed formula Phi_1(A) = sum_{k=2}^n |a_{1k}|^2 per A[hat{1,k}].  The paper's ATTRIBUTION of that inequality -- that it is Zhang's Eq. (21) at p. 315, that it is the weak consequence of Eq. (20) which survived Hutchinson's disproof of the stronger conjecture, and that under Hutchinson's indexing B_{k-1,k-1} = A[hat{1,k}] -- is transcribed from the cited articles and is NOT checked here: this program has no network access and no copy of those articles, so it cannot confirm the equation numbers, the page, or the index convention.  If that identification were wrong, the arithmetic below would still hold but it would refute a differently-attributed statement.
> the identities the proofs rest on
