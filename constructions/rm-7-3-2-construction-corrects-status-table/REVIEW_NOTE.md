# A Certified (7,3,2) Mixed Radial Moore Graph

`rm-7-3-2-construction-corrects-status-table`

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
passes. The recorded run reports **45 checks, all passing**:

    VERDICT: ALL 45 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    d5b8d73e8974bee16e1f115cb86dbccc2208fc1507f01977cacc0aa85561e811

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE, stated explicitly:
> (a) the bibliographic claim that a particular printed table
> carries the entry 18 in the cell (r,z) = (7,3) is a fact
> about published pages and cannot be checked by this
> program.  Concretely, the row and column headers of the
> two disputed tables are NOT reproduced anywhere in this
> program or in the paper, so nothing here excludes a table
> whose axes are indexed the other way round, i.e. by
> (z,r); on that reading the 'transposed cell' correction
> would itself be a misreading of a correctly printed
> table, and no check above would detect it.  The same
> holds for the paper's supporting quotation that the
> source text names the vector s_{3,7,2}.  What IS checked
> is the parameter arithmetic, which uses no table layout
> at all: order 108 with baseline 204 forces undirected
> degree 3 and directed degree 7, whereas (r,z) = (7,3) has
> order 104 and baseline 196.
> (b) no census over an external catalogue of mixed graphs is
> attempted, and no minimality of N_1 = 3378 is tested;
> the paper makes no optimality claim either.
> (c) the exhaustive parameter search above ranges over the
> pairs (r,z) with 0 <= r,z <= 40, which covers every pair
> with M(r,z,2) <= 1681 and so both orders in question.
