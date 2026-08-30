# A 24-cell weak 4-resolving set of K_11 K_11 refutes a surmised exact value

`a-24-cell-weak-4-resolving-set-of-k11-times-k11-refutes-a-surmised-exact-value`

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
passes. The recorded run reports **109 checks, all passing**:

    VERDICT: ALL 109 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    c3f438ddeb1a37bb559cafd3a9ab36c1533e68cf77f5adf176ac2f8b4c018bfd

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: n = 9 and n = 10 are NOT decided here.  This program computes the lower bounds 20 and 22 and confirms the published upper bounds 21 and 23, and performs no search in between; the surmised equality may well be correct at those two values of n.
> NOT RE-RUN: no exhaustive or solver search over candidate sets is performed at any n.  Minimality comes only from Theorem 1.2, whose combinatorial steps are hand proofs; this program checks that theorem's ARITHMETIC (the chain 2n <= 4h+3E <= 11E for n = 4..300) and not its Steps 1, 4 and 5.
> NOT RE-RUN: the block-diagonal family is verified against the raw distance rule only at m = 1, 2, 3, and against the Lemma 2.2 criterion at m = 1..6.  The statement for all m >= 1 is the hand proof of Section 5.
> NOT RE-RUN: no text of arXiv:2605.22307v1 is fetched or hashed here, so the quotations, line and byte locators, statement numbers and file digests of Section 1 are NOT machine-checked by this program.  The one numerical value taken from that source on trust is wdim_4(K_3 x K_3) = 6.
> NOT RE-RUN: the constraint-programming runs that originally found S_12, S_13, S_14 and a second 24-cell set at n = 11 are not repeated; they are needed only as a source of upper-bound witnesses, and the witnesses shipped are re-verified above from scratch.  No solver is imported.
> NOT RE-RUN: exact values for n outside {11, 12, 13, 14} and outside the multiples of 11 are neither computed nor bounded better than by 2n+ceil(2n/11) <= wdim_4 <= 2n+1+floor(n/4).
