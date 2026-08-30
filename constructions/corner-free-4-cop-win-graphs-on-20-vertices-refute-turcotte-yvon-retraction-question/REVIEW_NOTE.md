# Corner-free 4-cop-win graphs on 20 vertices,\\ and a retraction question of Turcotte and Yvon

`corner-free-4-cop-win-graphs-on-20-vertices-refute-turcotte-yvon-retraction-question`

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
passes. The recorded run reports **51 checks, all passing**:

    VERDICT: ALL 51 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    2e4b4832507e03107ad7a95ec6fcec9edd03fa9ed412da9ce8d45ad850f29a76

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the completeness of the census lane.  This program generates no graphs, so it does not re-derive that the connected 4-regular graphs of girth at least 5 on 20 vertices are exactly two.  That count is Meringer's published census (file 20_4_5.asc) and OEIS A058343; what is checked here is that the two graphs named in the paper have the claimed automorphism orders, 96 and 20.
> NOT RE-RUN: M_4 = 19.  It is the target paper's own corollary and is used, not verified; the paper's claim is conditional on it in exactly the way the source states it.
> NOT RE-RUN: girth 3 and girth 4 graphs on 20 vertices.  Nothing here bounds how many 20-vertex graphs with cop number 4 and no corner exist, so no minimality or uniqueness is claimed.
> NOT RE-RUN: k = 2, k = 3 and k >= 5 of the question.  k = 2 is refuted in the source itself, k = 3 is TRUE there by its own classification, and k >= 5 is untouched because M_5 is unknown.  Only k = 4 is settled here.
> NOT RE-RUN: the authors' published table of 25148 cop numbers.  The solver above is controlled on five published values, not regression-tested against that table.
