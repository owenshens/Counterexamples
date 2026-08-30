# A Proof of Calderon's mod-p Bailey Congruence at Every Odd Inert Prime

`a-proof-of-calderons-mod-p-bailey-congruence-at-every-inert-prime`

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
passes. The recorded run reports **82 checks, all passing**:

    VERDICT: ALL 82 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    35d79341b188b3ab68e88caaba3f71b67b4d311c5f25130bbe2e728effa65889

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the paper's Theorem 1 is a hand proof and is NOT machine-verified here. No proof assistant was used. Sections 3-4 above check the finite combinatorial content of its Steps 1, 2 and 4 and of its integrality Lemma over stated finite ranges, and Sections 6-11 corroborate the conclusion on finitely many instances; neither is a proof.
> NOT RE-RUN: the split cells of Table 5 at their FULL digit range. Those would be 78*78*9 = 54756 instances at p = 13 and 55*55*9 = 27225 at p = 11. Table 5 is bounded to alpha,beta <= 4 and alpha,beta <= 6 and the bound is printed in the table.
> NOT RE-RUN: the ramified prime p = 2 is only formally exercised: 1 <= gamma <= alpha <= p-1 = 1 forces alpha = gamma = beta = delta = 1, so its 9 instances are degenerate and establish nothing about ramified primes.
> NOT RE-RUN: no minimality or optimality claim is tested. Nothing here searches for a counterexample outside the inert regime beyond the three families of Table 5, and nothing here bears on Conjectures 6.3 and 6.4 of the source (the reciprocal-moment and mod-p^{3k} Ljunggren statements), which are untouched and remain open.
> NOT RE-RUN: no prior-art search. This program makes no network access and reads no file; the standing of the result against the literature is argued in the paper, not here.
> NOT RE-RUN: the earlier working programs of this project are not executed by this file. verify.py was written independently against the paper and reproduces their counts.
> checks run: 82 passed, 0 failed, elapsed 28.3s
