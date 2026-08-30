# The Cyclic Antibandwidth of the 3x3x4 Mesh is 13, Not the Conjectured 14

`cyclic-antibandwidth-of-the-3x3x4-mesh-is-thirteen-not-the-conjectured-fourteen`

Supporting material for this paper: the program that checks its computational claims, a second
program that re-decides the paper's decisive negative with the one non-trivial step of its proof
removed, and a recorded run of each.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |
| `exhaust.py` | auxiliary program: the `k = 14` exhaustion with no symmetry reduction at all |
| `exhaust.output.txt` | recorded run of `exhaust.py` |

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package, no external data file, no
network and no randomness, and no floating-point value enters any decision. It is deterministic —
every search-tree node count it prints is reproducible exactly — and it runs in about three
minutes on one core, nearly all of that in the `k = 14` exhaustion of its section 10. The program
prints one line per check and a closing verdict, and exits 0 only if every check passes. The
recorded run reports **56 checks, all passing**:

    VERDICT: ALL 56 CHECKS PASS

Its inputs are hard-coded objects: two 36-label grids `W1` and `W2` — the first of which is the
labeling `W` printed in Section 2 of the paper — the degree-6 calculation at `(2,2,2)` printed
beside it, and the conjectured row as a closed form. The program checks more than the paper
claims: besides `W`, it also reads the second grid, decides `k = 15` separately, and does
arithmetic on six reported values `cab` and `ab` for `n3 = 3..8`; the paper states none of that,
and none of it is needed for the paper's theorem. What matters here is that **it re-derives every
layer of the paper's proof, with no solver**:

* the graph, by two independent constructions (axis successors, and `L1` distance exactly 1);
  `|V| = 36`, `|E| = 75`, `Delta = 6` at exactly `(2,2,2)` and `(2,2,3)`, the 18/18 bipartition,
  and — via a maximum matching — the independence number 18, which is what shows that no
  independence-number argument can decide `k = 14`;
* that `W1` and `W2` are bijections onto `1..36` with minimum cyclic edge distance exactly 13, on
  15 and 16 tight edges respectively, hence `cab >= 13`;
* the arc identity `|{l' : d_c(l,l') >= k}| = n - 2k + 1` for all 36 labels and all `k <= 18`, and
  the window reformulation (Lemma 3 of the paper) verified over **all** bijections of two small
  graphs — all 40,320 of `P_2x2x2` and all 720 of the 2x3 grid;
* the degree cap, which kills every `k >= 16`;
* `k = 15` INFEASIBLE by complete exhaustion, twice: 843,415 search-tree nodes with **no**
  symmetry reduction and 202,414 with it, 0 complete labelings both times;
* `k = 14` INFEASIBLE by complete exhaustion: 281,957,255 nodes, 0 complete labelings, ~156 s.

The decision procedure is calibrated in both polarities on integers the pipeline did not produce,
before any negative is read: `cab(P_2x2x3) = 4` and `cab(P_3x3x3) = 9`, each decided two-sidedly
(a witness found and re-verified at the published value, infeasibility one above), the latter
being the one cell of the target row that a published exact method certifies. There is also a
proved-silent control (`k = 12` on the 27-vertex mesh, above its cap) and a check that the degree
cap does *not* decide `k = 14` or `k = 15`, so the searches are not decorative.

The only step of the upper bound that needs an argument rather than a run is the 16-fold symmetry
reduction, and the program attacks it from two sides: it exhibits that **exactly one** of the 16
images of `W1` under (eight grid automorphisms) x (label reflection) satisfies the three
canonicity constraints — so the reduction discards no solution and keeps no duplicate — and it
replays the reduced search along that canonical image, confirming that the label survives
forward checking and all three constraints at each of the 36 placement steps. `exhaust.py` then
removes the reduction entirely.

## The auxiliary program

```sh
python3 exhaust.py
```

Same contract (Python 3.9+, standard library only, deterministic). It decides one question:
`cab(P_3x3x4) >= 14`, by complete exhaustion with **no symmetry reduction beyond the rotation
pin**, which is a proof rather than a heuristic. It therefore pays a much larger tree than
`verify.py` pays — 1,399 s against 156 s, a measured factor of 8.74 in nodes — and it agrees:
`INFEASIBLE`, 2,463,942,873 search-tree nodes, 0 complete labelings, the search dying 4 vertices
short of a complete labeling on every branch. Its recorded run closes with

    VERDICT: ALL 4 CHECKS PASS

If the two runs had disagreed, the symmetry argument would be wrong.

## Provenance

`verify.output.txt` and `exhaust.output.txt` each hold their program's output, preceded by a short
header and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the pairs can be matched:

```sh
shasum -a 256 verify.py exhaust.py
```

    7f02fa53909e089ef8ce504a6743f7ec916b5247ed42480e1259ad1f754d6f2f  verify.py
    9fe82e9d862fe2cb828f8aa5367c5524f60ac06d453ff8d5d745fa8a3d30578a  exhaust.py

## Scope

The programs' own closing statements of what they do not cover, quoted from `verify.py`'s output:

> NOT RE-RUN: a HAND-CHECKABLE proof that k = 14 is impossible. The k = 14 layer is decided above,
> exhaustively and without any solver, but it is decided by a search of a few hundred million
> nodes: a referee must re-run it, not read it. No DRAT, LRAT or other independently checkable
> UNSAT proof object exists for it anywhere in the record, and the obstruction is nearly global
> (its minimal infeasible induced subgraph was measured elsewhere at 30 of the 36 vertices), so it
> does not compress to a hand-sized core the way k = 15 does. The 16-fold symmetry reduction used
> here IS argued, and its soundness is checked in section 8; exhaust.py in this folder re-decides
> k = 14 with that reduction removed entirely.

> NOT RE-RUN: the cells n3 = 5, 6, 7, 8. The values 18, 22, 27, 31 quoted in the paper are
> reported, not re-derived here: no labeling attaining them survives in the record, so there is no
> object for this program to read. Only n3 = 4 is settled by what is printed.

> NOT RE-RUN: the linear antibandwidths ab(P_3x3xn3) = 9, 14, 18, 23, 27, 32. Their ARITHMETIC
> relation to the published row is checked above; their computation is not reproduced, again for
> want of a printed witness.

> NOT RE-RUN: anything at n3 >= 9. The uniform closed form floor(9(n3-1)/2) is a conjecture; 392 of
> the 398 cells the published row quantifies over are untouched.

> NOT RE-RUN: the value 1779 attributed to a published heuristic at n3 = 400, and the transposition
> of two cells in the restating paper's results table. Only the arithmetic around them is checked;
> neither is read off the source by this program.

> NOT RE-RUN: the bibliographic locators. The arXiv identifier, the DOIs, the table label and the
> source line numbers printed in the paper are not verified here; this program checks mathematics
> only.

> NOT RE-RUN: prior art. Nothing here bears on whether cab(P_3x3x4) = 13 was known.

Two further limits belong to the paper rather than to the programs. The body of the 2011 primary
source is closed access and was read by nobody here, so what is refuted **by name** is its 2026
restatement, quoted in Section 1 of the paper; because the refutation is downward, that gap cannot
rescue the row on any reading of the primary. And the row is not repaired: one of the 398 cells it
quantifies over is decided, and the paper claims nothing about the other 397.
