# Reduce-Max Exceeds Backlog 2 on Ten Bamboos

`ten-bamboo-counterexample-reduce-max-backlog-bound`

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
passes. The recorded run reports **21 checks, all passing**:

    VERDICT: ALL 21 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    72a5002681ab6626a046001af94bab6fff1bd1be361efa2025d46094a98be7ad

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: (a) the 2000- and 2702-bamboo instances of the earlier disproof, and their reported bounds 2.0004 and 2.076, are not reconstructed, since those rate vectors are not exhibited in the text being checked; (b) the corollary is simulated exactly for n = 10..22 only, its general growth bound t/(1000m) <= 60/1000 being verified in closed form for all m >= 1; (c) nothing about fewer than ten bamboos is claimed or tested; (d) in the census of H <= 35 the backlog is exact only for the partitions whose state recurred within the 4000-round horizon, is a horizon-limited lower bound for the remaining 95, and is computed under lowest-index tie-breaking alone -- that census only corroborates the remark about the earlier search and carries no part of the main theorem, whose tie-independence is established by exhaustive expansion above.
