# Review note

Files in this folder: `paper.tex`, `paper.pdf`, `verify.py`, `verify.output.txt`, and this note.
Nothing outside the folder is named below.

## 1. What the paper claims

*A Matching Labelling for the Cyclic Antibandwidth of the 2x2xn Mesh* proves, as Theorem 1, that

    CAB(P_2 x P_2 x P_n) = 2(n-1)   for every integer n >= 2,

equivalently for the cylinder `C_4 box P_n`, and records that the identity fails at `n = 1`, where
`CAB(P_2 x P_2 x P_1) = CAB(C_4) = 1`.

The statement it settles is displayed as equation (1) of Section 1, "The statement": the first data
row of a table of conjectured three-dimensional mesh values attributed to Bansal and Srivastava,
asserting `CAB(P_2 x P_2 x P_{n_3}) = 2(n_3-1)` for `n_3 = 3,...,500`. Theorem 1 gives equality at
all 498 tabulated cells and for every `n_3 > 500` besides. Section 1 states that the 2011 original
was not read and that the row is quoted from the e-print arXiv:2601.04239v1, which reproduces and
attributes the table; that two of the 498 cells (`n_3 = 3` and `n_3 = 168`, values 4 and 334)
already carry certificates of optimality in that benchmark paper; and that 667 and 997 are the
entries recorded in the OPTSICOM repository at `n_3 = 335` and `n_3 = 500`.

The proof is two hand arguments, both short. Proposition 2 (Section 2) exhibits the closed-form
labelling, equation (2),

    f(u,v,i) = (1 + 2i + u) + 2n * ((u+v+i) mod 2),

shows it is a bijection onto `1..4n`, and computes its edge distances exactly: `2n` on every
`v`-edge, `2n-1` on every `u`-edge, `2n-2` on every `i`-edge, so the minimum is `2n-2` for `n >= 2`.
Lemma 3 (Section 3) is the counting bound: if `CAB(G) >= k` with `1 <= k <= N/2` then every vertex
has degree at most `N-2k+1`. Theorem 4 applies it at `k = 2n-1` to a degree-4 vertex of `G_n`,
`n >= 3`, obtaining `4 <= 3`. Proposition 5 settles the endpoints, `CAB(C_4) = 1` and
`CAB(Q_3) = 2`, the latter by showing that a labelling attaining 3 would make `Q_3` isomorphic to a
non-bipartite circulant. Section 4, "Two instances written out", prints the labelling at `n = 3` and
`n = 4` and evaluates the formula at `n = 335` and `n = 500`: 668 on 1340 vertices and 998 on 2000
vertices, one above each of the entries 667 and 997.

Section 5, "What is already in print, and what is not settled", limits the novelty claim itself: the
upper bound is also a consequence of the published two-dimensional formula
`CAB(P_{n_1} x P_{n_2}) = n_2(n_1-1)/2` at `n_2 = 4`, since `P_n x P_4` is a spanning subgraph of
`C_4 box P_n` with exactly `n` fewer edges; and the matching labelling of the cylinder is not claimed
as new either, the 2005 and 2009 references carrying that formula not having been read in full.

## 2. What the program checks

`verify.output.txt` is the recorded run of `verify.py`. Its closing line is
`VERDICT: ALL 42 CHECKS PASS` and the run exited with status 0. The 42 `PASS` lines fall in seven
labelled blocks:

* **A, the label tables printed in the paper, read literally (12 checks)** — Section 4, and the edge
  classification of Section 2. For each of the `n = 3` and `n = 4` tables: bijectivity onto `1..12` /
  `1..16`, agreement of every printed label with the closed form (2), the edge counts `20 = 8*3-4`
  and `28 = 8*4-4`, the distance multisets `[(4,8),(5,6),(6,6)]` and `[(6,12),(7,8),(8,8)]`, and the
  minima `4 = 2*3-2` and `6 = 2*4-2`. Two structural checks back the rest: mesh and cylinder have the
  same edge set for `n = 2..30`, and the three `O(n)` edge families used by the sweep coincide with
  the literal all-pairs product adjacency predicate for `n = 2..14`.

* **B, the sweep `n = 2..500` (8)** — the arithmetic of Proposition 2, cell by cell. 499 values of
  `n`, containing all 498 cells `n_3 = 3..500`: `|V| = 4n`, `|E| = 8n-4`, the degree multiset (eight
  3s and `4(n-2)` 4s for `n >= 3`, 3-regular at `n = 2`), bijectivity of (2), the exact distance
  multiset `{2n-2: 4(n-1), 2n-1: 2n, 2n: 2n}`, and hence the minimum `2n-2` everywhere. Plus the two
  largest cells written out: 1340 vertices, 2676 edges, minimum 668; 2000 vertices, 3996 edges,
  minimum 998.

* **C, the counting upper bound by brute force in `Z_{4n}` (9)** — Lemma 3 and Theorem 4. Exactly
  `3 = N-2k+1` residues lie at cyclic distance `>= 2n-1` from a given one at every swept `n`; the
  count is monotone in `k`, so ruling out `k = 2n-1` rules out every larger `k`; the side conditions
  `k <= N/2` and `2(k-1) < N` hold; `Delta(G_n) = 4` at `n = 3..30` and by the sweep degree multiset
  beyond; hence the contradiction `4 > 3` at every `n = 3..500`, with attained minimum and ceiling
  both reading `2(n-1)` there. Also the four benchmark cells (`3, 168, 335, 500` giving
  `4, 334, 668, 998`), the comparison against 667 and 997, and `P_3 x P_3 x P_3`: `N = 27`,
  `Delta = 6`, ceiling 11 — the witness Section 5 uses for the bound not being tight in general.

* **D, complete rotation-pinned search, the constraint as the only prune (5)**, so the upper bound is
  tested rather than assumed. Pinned solution counts at `k = 2n-1` are `0` for every `n = 3..8`;
  optimal-labelling counts `(n, pinned, total)` are `(3,8,96)`, `(4,8,128)`, `(5,8,160)`; the
  endpoints of Proposition 5 come from unpruned enumeration, all `3! = 6` pinned labellings of `C_4`
  giving 1 and all `7! = 5040` of `Q_3` giving 2; and the pruned engine is cross-checked against full
  enumeration on `Q_3` at `k = 2`.

* **E, negative controls (3)** — the checker and the search must be able to say NO. A naive
  layer-by-layer labelling scores minimum 1 at `n = 3, 4, 5, 50`; swapping labels 1 and 2 in the
  `n = 4` witness drops the minimum to 5; the closed form attains `2n-2` and *not* `2n-1` at
  `n = 3, 4, 5, 10, 100`.

* **F, the two-dimensional route of Section 5 (3).** `floor(n_2(n_1-1)/2) = 2(n_1-1)` at `n_2 = 4`
  for `n_1 = 4..500`, so the floor is vacuous there; `P_4 box P_n` is a proper subgraph of
  `C_4 box P_n` with exactly `n` fewer edges for `n = 2..40`; `(8n-4) - (7n-4) = n` at every
  `n = 2..500`.

* **G, two tori by complete search (2)**, bearing on the second priority caveat of Section 5, that a
  published toroidal value `2n-2` would supply Proposition 2. `C_4 box C_3` is 4-regular on 12
  vertices, the Lemma 3 ceiling is 4, and the search finds nothing at `k = 5`; `C_4 box C_4` has no
  labelling with every edge at cyclic distance `>= 6 = 2*4-2`.

One convention when reading the transcript beside the paper: its block headings and closing lines
call the closed-form labelling "Proposition 1" and the counting bound "Lemma 2", whereas the shipped
paper numbers them Proposition 2 and Lemma 3. Equation (2) and the Section 5 reference agree; no
mathematical content differs.

## 3. What the program does not check

**Theorem 1 is a hand proof and the program is a control.** No block establishes the theorem: the
exhaustive search of block D independently decides only `n = 3..8`, plus the endpoints `n = 1` and
`n = 2`. The claim on all 498 cells, and a fortiori every `n_3 > 500`, rests on Proposition 2 and
Lemma 3 as proved by hand; blocks B and C re-derive their arithmetic cell by cell, and the searches
are corroboration, not the load-bearing step. The transcript says exactly this in its closing
`NOT RE-RUN` lines, and the paper's "Scope" paragraph in Section 5 says that beyond `n = 500` the
statement of Theorem 1 rests on those two proofs alone.

The quantifiers are proved for all `n >= 2` and only sampled by the program: the sweep at
`n = 2..500`, the cylinder identification at `n = 2..30`, the literal adjacency predicate at
`n = 2..14`, `Delta = 4` directly at `n = 3..30`, the subgraph relation at `n = 2..40`, the complete
search at `n = 3..8`, the over-attainment control at `n = 3, 4, 5, 10, 100`.

Carried over from the transcript's closing `NOT RE-RUN` block:

* **The quoted source statement itself.** The conjectured row, the certified optima 4 and 334, and
  the best-known values 667 and 997 are facts about external documents. The program checks the
  mathematics *against* those integers; it does not fetch them. They are transcribed from the cited
  sources, not recomputed, and the paper says so in Section 1 and in the Section 5 paragraph
  "Sources not read, and what rests on them".
* **The priority question.** Nothing in the run can tell whether an extremal labelling of the
  `4 x n` grid, or a toroidal formula, already yields the closed-form labelling. The paper does not
  claim the lower bound as new.
* **The exhaustive census beyond `n = 8`.** Block D decides `n = 3..8` and counts the optima at
  `n = 3, 4, 5` only.
* **The count `32n` of optimal labellings.** It is computed at `n = 3, 4, 5` only. The paper states
  no such count.
* **`CAB(C_4 box C_n)` beyond `n = 4`.** Block G searches `n = 3` at `k = 5` and `n = 4` at `k = 6`
  only.
* **The other four rows of the conjectured table** (3x3, 4x4, 5x5, 6x6). Only the ceiling 11 at
  `P_3 x P_3 x P_3` is computed, as the witness that the bound is not tight in general; Section 5
  states that nothing here transfers to those rows.

## 4. How to check it

```
shasum -a 256 verify.py
python3 verify.py
```

`verify.py` is standard library only (`re`, `sys`, `itertools`, `collections`), reads no input file
and no network, uses no floating point, and exits 0 if and only if every check passes. The recorded
run used Python 3.9.25. The header of `verify.output.txt` carries the SHA-256 of the program it ran,
so transcript and program can be paired; that digest is

```
5420af9a710b77e9bcefb57b96722b40e52ffd792f5712cc8b84c0bf474f362b  verify.py
```

and it is the digest of the `verify.py` shipped in this folder.
