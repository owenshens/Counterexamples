# A 32-Digit Prime x with x^2+x+1 Square-Full: A Question of Ochem and Zelinsky Answered in the Affirmative

`a-32-digit-prime-x-with-x2-plus-x-plus-1-square-full-answers-ochem-zelinsky`

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
passes. The recorded run reports **74 checks, all passing**:

    VERDICT: ALL 74 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    289a9b2f3da65350b94ccf993d60a9b0d55d9a10953702fc159d936c78cfc338

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE -- what this program does NOT cover. (i) Nothing bibliographic is re-derived: the locator of the question (restrained.tex line 852 of the arXiv:2607.19746v1 e-print), the provenance of the witness integer as a(27) of OEIS A296376 and line 27 of its b-file b296376.txt, and the 1902 Majol reference are transcribed from those sources and are not fetched or checked here, and neither is the (negative) fact that x=18 is NOT an example recorded in that e-print.  (ii) NOT RE-RUN: any x > 10^8. The census in Step 7 decides every x with 1 <= x <= 10^8 and finds no prime witness there, so no claim of minimality is made or supported for the exhibited 32-digit x; a smaller prime witness with 9 to 31 digits is not excluded by anything here.  (iii) NOT RE-RUN: whether infinitely many primes x have x^2+x+1 square-full, and whether the two branches exhibited in Step 5 exhaust the solutions of u^2 - 1372 a^2 = -3 -- Step 5 checks 20 members of one branch and 4 of the other, which supports the remark about why the object is hard to find and proves no structure theorem.  (iv) The Lucas certificates prove primality of x and of the factors of n; they say nothing about whether the factorisation of n printed in the paper is the only one, beyond the re-multiplication check in Step 1.
