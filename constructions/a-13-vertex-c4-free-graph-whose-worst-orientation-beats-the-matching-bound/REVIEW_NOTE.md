# A 13-Vertex Graph With No 4-Cycle Whose Worst Orientation Beats the Matching Bound

`a-13-vertex-c4-free-graph-whose-worst-orientation-beats-the-matching-bound`

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

Python 3.9 or later, standard library only: no third-party package, no external data file, no
network and no randomness, and no floating-point value enters any decision. The recorded run
reports no timing. The program prints one line per check and a closing verdict, and
exits 0 only if every check passes. The recorded run reports **57 checks, all passing**:

    VERDICT: ALL 57 CHECKS PASS

Its inputs are the objects printed in the paper and nothing else: the graph6 string
`L??E@agXAp@wdC`, the 22-edge list and adjacency lists of Section 2, the six-edge matching and
the independent set `A = {0,...,5}` of Section 2, the seven traces of the orientation `D*` of
Section 3, and the published `mold` values of the source paper. Every
quantity the paper states is re-derived from those: that the graph6 string is **label-equal**
(not merely isomorphic) to the printed edge list, with 78 bits consumed and no padding bit set;
the adjacency lists, the degree sequence and connectedness; that all 78 vertex pairs have at
most one common neighbour and, by an independent route, that none of the 715 four-element vertex
subsets carries a 4-cycle; the two triangles and the girth; `alpha'(G) = 6` and `alpha(G) = 6`
by exhaustive search, hence `n - alpha'(G) = 7` and `alpha + alpha' != n`; that `A` is a set of
sources in `D*` for each of the 16 ways of orienting the four edges inside `V \ A`, that its
seven traces are exactly the printed ones, and that no set of size at most 5 is
locating-dominating in `D*`; that **all 2^22 = 4 194 304 orientations** of `G` admit a
locating-dominating 6-set; that exactly `3 192 596` of them admit none of size 5; and finally
`mold(G) = 6 < 7 = n - alpha'(G)`.

The load-bearing claim is `mold(G) <= 6`, and it is decided **exhaustively**, not sampled. A
naive sweep would be about 10^9 predicate evaluations, out of reach for CPython, so the program
builds for each candidate 6-set `S` the *exact* set of orientations in which `S` is
locating-dominating, as a 2^22-bit Python integer, and unions those sets; 980 of the 1716
six-sets already union to all of `{0, ..., 2^22 - 1}`. The reduction is stated in Section 4 of
the paper, so a referee can check the algorithm and not only run it: the predicate depends on
the orientation only through the traces of the seven vertices outside `S`, each "empty trace" and
each "two traces collide" event is a cylinder in the orientation cube, and the no-`C_4`
hypothesis is what makes at most one collision cylinder per pair.

That machinery is **calibrated in both polarities before any conclusion is drawn**, on other
people's integers. It reproduces every `mold` value we could locate in the source paper:
`mold(C_4) = 3`, `mold(P_4) = 2`, `mold(C_n) = ceil(n/2)` for `n = 3` and `5 <= n <= 10`,
`mold(K_n) = ceil(n/2)` for `3 <= n <= 7`, `mold(P_n) = ceil(n/2)` for
`n = 2,3,4,5,7,9`, and `mold(K_{1,4}) = 4`. Two of those discriminate conventions:
`mold(C_4) = 3 > 2 = n - alpha'(C_4)`, which alone rules out a minimum-over-orientations
misreading of `mold`; and the Petersen graph, which lies in the **same class** as the witness (no
`C_4` subgraph, maximum codegree 1 over all 45 pairs) and **attains** the bound,
`mold = 5 = n - alpha'`, so the decision procedure demonstrably answers "attains" where it must.
The bowtie, `C_5`, `C_7`, `C_9` and `K_7` also attain. Independently, on all 19 control graphs
with at most ten edges the bitset construction is checked against a **naive per-orientation
implementation** that shares no code with it, with no disagreement; and on the 13-vertex witness
the two agree on 2048 sampled orientations at `k = 6`, where a locating-dominating 6-set is
found, and on 256 sampled orientations at `k = 5`, where the naive scan over all 1287 five-sets
confirms that none exists. The informative non-zero `3 192 596` on the same code path that
returns `0` at `k = 6` is what shows the answer `0` is not a silently empty sweep.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    a905912f956fc52e4126745d4dc27be46b2614702cbb5080f13ccc0e571beb91

## Scope

The program's own statements of what it does not cover, quoted from its output:

> NOT RE-RUN: minimality of any kind. This program says nothing about whether 13 is the least
> order, or 22 the least size, of a graph with these properties, and it does not reproduce the
> census that suggested so; that census found 32 such graphs at n = 13, fourteen of them at
> m = 22, and its per-graph decisions were never re-run by a second implementation. The paper
> claims existence only, which needs none of it.

> NOT RE-RUN: the other thirteen graphs at (n, m) = (13, 22), and the two at m = 24. Only the
> single witness printed in the paper is decided here.

> NOT RE-RUN: n >= 14. No graph of order 14 or more is examined, so nothing here bounds how
> common such graphs are.

> NOT RE-RUN: the bibliographic locators. The line numbers, the page, the problem number "Open
> problem 37" and the DOI printed in the paper are not checked by this program; it checks
> mathematics only.

> NOT RE-RUN: prior art. Nothing here bears on whether the answer was known.

> NOT RE-RUN: ld(G) and the proof of the source's Lemma 34. The chain ld <= mold <= n - alpha'
> is not re-derived; only the two ends that the strict inequality needs, mold(G) = 6 and
> n - alpha'(G) = 7, are computed.

Two limits of the surrounding record, beyond what the program can see. First, the paper's
"Scope" section is load-bearing and should be read: the witness is **not** claimed
smallest or unique, and the search that reported 32 witnesses at order 13 is not reproduced
anywhere in this folder and is not mentioned in the paper. Second, the recorded run above was made on the machine that prepared
this folder (Python 3.9.6), not on a second independent machine; the program is deterministic
and depends on nothing outside the standard library, so any reader can produce the same bytes,
but this transcript is not itself an independent replication.
