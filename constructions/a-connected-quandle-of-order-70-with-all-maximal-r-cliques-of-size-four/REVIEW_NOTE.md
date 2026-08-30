# A Connected Quandle of Order 70 All of Whose Maximal R-Cliques Have Size Four

`a-connected-quandle-of-order-70-with-all-maximal-r-cliques-of-size-four`

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
passes. The recorded run reports **85 checks, all passing**:

    VERDICT: ALL 85 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    d5fbe36a7d2595a791d31c31e4b3c00e69d4f9ea4d0e4cbea08edaa31ccfb6d0

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the minimality question.  This program does NOT enumerate connected quandles by
> NOT RE-RUN: order, so it neither confirms nor contradicts the claim that no connected quandle of
> NOT RE-RUN: order at most 24 is a witness; that census is reported in the paper as unreproduced,
> NOT RE-RUN: orders 25-69 were never examined, and connected RACKS THAT ARE NOT QUANDLES were
> NOT RE-RUN: never examined at any order.  Nothing here asserts that 70 is the least order.
> NOT RE-RUN: the version of record.  The problem statement and its locator in S1 are quoted from
> NOT RE-RUN: arXiv:2504.09368v1; the journal text (J. Algebra 698 (2026) 493-532) was not
> NOT RE-RUN: accessible to us and its problem numbering may differ.  No program can check that.
> NOT RE-RUN: Part C sweeps the commuting graph for 5 <= n <= 12 and the full rack identity only for
> NOT RE-RUN: n <= 8; for n = 13 and beyond only the closed-form arithmetic (C101-C108) is checked.
> NOT RE-RUN: no claim about Problem 5.22 (do all maximal R-cliques of a connected rack have the
> NOT RE-RUN: same size?) is made or tested here; in Q every maximal R-clique has size 4.
