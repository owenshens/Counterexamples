# An explicit (6,4) and an explicit (8,4)

`explicit-of-6-4-and-of-8-4-close-the-even-g-4-cells-of-distance-three-one-factorizations`

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
passes. The recorded run reports **56 checks, all passing**:

    VERDICT: ALL 56 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    326c6211683fb6d98d0e81a649b9a3953b8b847b2deeb9092b23f60c7f520f62

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the source's doubling lemma itself.  Checks F1-F2 verify its stated hypotheses and side conditions at (n,g) = (6,4); the lemma is QUOTED, not reproved, and NO OF(12,4) is constructed here.  Only OF(6,4) and OF(8,4) are exhibited and verified edge by edge by this program.
> NOT RE-RUN: the clauses labelled "source:" in F4/F5.  They are the target paper's own results (OF(g+1,g) for even g, odd n >= 2g-1, the doubling lemma, and the n = 0 mod 4 proposition) and are taken on trust; F4 checks only that the six clauses TOGETHER (four of the source, two of this paper) leave no n > 4 uncovered, and it checks it over the finite range 5 <= n <= 4000; beyond it the clauses are residue conditions and the coverage argument is the paper's, not the program's.
> NOT RE-RUN: MINIMALITY and UNIQUENESS.  No census of the prescribed-automorphism lane was run here, neither base factor is claimed unique or smallest, and no OF is claimed to be the only one on its parameters.
> NOT RE-RUN: g >= 5, the range n <= g (the source's ODAR objects), and Hamming distances other than 3.  Nothing in this program touches them.
> NOT RE-RUN: the wording, numbering and line offsets quoted from arXiv:2602.16319v1, and the prior-art channels (arXiv API, Crossref, Semantic Scholar, zbMATH) that bound the novelty.  MathSciNet was never consulted and OpenAlex never answered.  This program checks mathematics, not provenance.
