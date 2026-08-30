# A Disconnected Graceful Bipartite Graph with No Near alpha-Valuation

`the-smallest-graceful-bipartite-graph-with-no-near-alpha-valuation`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run. The decisive claim -- that the $7$-vertex, $6$-edge graph
$C_4 \cup P_3$ is bipartite and graceful and that no graceful labelling of it is a near
alpha-valuation -- is proved in the paper by hand, in a table of fifteen cases, and needs no
program at all.

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
Every quantity it handles is an integer; there is no floating-point arithmetic anywhere in it.
It reads the objects printed in the paper -- the edge list and labelling of $C_4 \cup P_3$, the
fifteen rows of the case table, the larger witnesses with their labellings, the census table
and the named controls -- and re-derives each of them. The program prints one line per check and a
closing verdict, and exits 0 only if every check passes. The recorded run reports **86 checks, all
passing**:

    VERDICT: ALL 86 CHECKS PASS

It takes a little over three minutes, almost all of it in the census of every graph with at most
ten edges.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    a88c05df5e29a21a56dbf63e6d1642155a097162e30b2e7910e0184002c4441e

## Scope

The program's closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the census at m = 11. This program re-derives the classification for m <= 10 only,
> so the paper's statement that C_4 u P_3, C_6 u K_{1,4} and 2C_4 u K_{1,3} are the ONLY
> counterexamples with at most 11 edges is verified here only up to 10 edges. The 11-edge witness
> 2C_4 u K_{1,3} is itself fully verified above, over all 4608 of its beta-valuations.

> NOT RE-RUN: any m >= 12. Nothing in this program, and nothing in the paper, bounds the number of
> counterexamples with 12 or more edges.

> NOT RE-RUN: the connected form of the question. It is OPEN. What is checked here is only that
> every connected bipartite graceful graph with at most 10 edges has a near alpha-valuation; its
> tree sub-case is an open conjecture of El-Zanati, Kenig and Vanden Eynden.

> NOT RE-RUN: vertex-minimality. Graphs on at most 6 vertices with more than 6 edges were never
> examined, so no claim is made or checked that 7 is the least possible number of vertices.

> NOT RE-RUN: the literature. That the non-existence half is new, and that the gracefulness of
> C_4 u P_3 is already published, are bibliographic statements no program can settle.

The first of these bears on how the result should be read, and the paper states it in its own
text. The counterexamples are all **disconnected**, so the refutation is of the question exactly as
posed, under the source paper's standing convention that graphs are not assumed connected; the
strengthened question for connected bipartite graphs is not settled by the paper, and the census
finds no connected counterexample in its range. The paper claims no census beyond $m \le 10$, so
the program's $m = 11$ caveat applies to material the paper does not assert; the program does
verify the $11$-edge witness $2C_4 \cup K_{1,3}$, which the paper does not use.
