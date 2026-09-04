# Referee note

Paper: *The Hexagonal Prism under a Closed-Neighbourhood One-Round Rule: Fractional Eternal
Domination Number 46/13*.

Files in this folder: `paper.tex` and `paper.pdf` (the paper), `verify.py` (a verification
program), `verify.output.txt` (a recorded run of it), and this note. Nothing else is needed.

## 1. What the paper claims

**Theorem 1.** For the closed-neighbourhood one-round rule read off in Section 1,
`gamma_f^infty(C_6 [] K_2) = 46/13`, and the infimum is attained: an eternal strategy of total
weight exactly `46/13` exists.

`C_6 [] K_2` is the hexagonal prism on 12 vertices. Section 1 identifies it as the cell `k = 1`
of the first bullet of Problem 7.5 of Devvrit, Krim-Yee, Kumar, MacGillivray, Seamone, Virgile
and Xu — reference [1] of the paper, Discuss. Math. Graph Theory **44** (2024), no. 4, 1395–1428,
arXiv:2304.11795 — and as the smallest cell of that problem left open there, [1] proving for it
only `7/2 < gamma_f^infty(C_6 [] K_2) <= 4` and, through its own connectivity estimate
`(n + kappa)/(kappa + 1)` at `n = 12`, `kappa = 3`, also `<= 15/4`.

**The claim is explicitly conditional on a reading of the movement rule.** The budget constraint
of [1] is the open-neighbourhood sum `sum_{y in N(x)} m_{xy,i} <= w_i(x)`; the paper reads
unmoved weight as staying where it is, so that a legal round transports weight along *closed*
neighbourhoods with budget `w_i(x)` at `x`. The abstract, Section 1 and item (1) of Section 4 all
state that this reading is the paper's own and not a quotation, that no update equation deciding
between the two readings was found in [1] and none is decided here, and that under a stricter
rule nothing in the paper bounds `gamma_f^infty` and the cell of Problem 7.5 stays open.

The two halves:

* **Upper bound (Section 2).** Lemma 2: a family `{f_v : v in V}` of fractional dominating
  functions indexed by `V`, each of total weight `S`, with `f_v(v) = 1`, and with every ordered
  pair joined by a closed-neighbourhood transportation plan, gives `gamma_f^infty <= S` with the
  infimum attained. It is applied to the 12 translates of the seed built from the printed
  profiles `h_0, h_1`, of total weight `1 + 9/13 + 24/13 = 46/13`. The seed, its twelve closed
  neighbourhood sums, and one full plan `f_(0,0) -> f_(3,0)` on 13 arcs are printed.
* **Lower bound (Section 3).** A finite rooted attack tree yields the linear program (1), and
  Lemma 3 says its optimum is at most `gamma_f^infty`. For the printed 9-node tree, (1) has 493
  variables, 201 equalities and 116 inequalities, and Proposition 4 exhibits a dual feasible
  point — 31 multipliers on the `>= 1` rows and 123 nonzero equality potentials, all of
  denominator 13 — of objective `17/13 + 29/13 = 46/13`, so weak duality gives
  `gamma_f^infty >= 46/13`.

An unnumbered remark closes Section 3: under this reading the last step of the published proof of
`7/2 < gamma_f^infty` is not justified, because from the `x = 0` member of the source's own
family, of total weight `7/2`, the attack at `j` is answerable, leaving `h = l = 1/4` where the
published argument asserts `0`. The published inequality is nevertheless true, since
`7/2 < 46/13`.

## 2. What the program checks

`verify.output.txt` records **100 checks, all passing** (`VERDICT: ALL 100 CHECKS PASS`, exit
status 0), in ten steps. Every object is rebuilt from printed values, in exact integer or
`Fraction` arithmetic.

| step | checks | which claim |
|---|---|---|
| 1. the graph | 7 | the Notation of Section 1: order 12, 3-regular, 18 edges, bipartite, and the printed `N[0], ..., N[11]`; plus `gamma = 4`, and that `Cay(Z_12, {+-1, 6})` is not bipartite and so is a different graph |
| 2. controls, run before any claim of the paper is touched | 14 | the witness of [1] on `Cay(Z_8, {+-1, 4})` reproducing `8/3` with all 64 ordered pairs reconfigurable; the published `28/5` at `n = 10`; the published lower-bound formula shown *not* to apply at `n = 6` (its hypothesis is `n = 10 mod 12`; at `n = 6` it would give `40/11 > 46/13` and refute the value); uniform `1/4` totalling 3; the `15/4` family confirmed a genuine eternal strategy; and one forced negative, a demand of `11/10` at vertex 0 refused from uniform `1/4` |
| 3. the family of Section 2 | 8 | Lemma 2 (i) and (ii) for the denominator-13 family: `h_0, h_1` even so the 12 translates are well defined; seed and its twelve closed-neighbourhood sums exactly as printed; total `46/13` and fractional domination for all 12 members; weight 1 at each member's own centre; the 12 translates pairwise distinct |
| 4. a second family | 7 | a further 12-member family of total weight `46/13` in denominator 26, whose seed is supplied in `verify.py` and is not printed in the paper (see §3) |
| 5. transportation plans | 24 | four plans arc by arc: every arc an edge or a self-loop, every vertex shipping exactly its source weight and receiving exactly its target weight, total shipped `46/13`, attacked vertex receiving exactly 1. One is the 13-arc `f_(0,0) -> f_(3,0)` printed in Section 2; the other three belong to the second family |
| 6. all ordered pairs | 3 | Lemma 2 (iii): all 144 ordered pairs of each family decided feasible by exact integer maximum flow, hence `gamma_f^infty(C_6 [] K_2) <= 46/13` |
| 7. the LP and its dual | 16 | Proposition 4: the 9-node tree with the printed parent and attack maps; the second attack set at distances 3, 4, 3; the dimensions 493 / 201 / 116 re-derived; the 31 printed inequality multipliers (6 attack + 25 domination), attack part `17/13`, domination part `29/13`, sum `46/13`, all of the sign the dual requires and none on the `>= 1` rows of nodes 3 and 8; the 123 equality potentials, their denominators, and the vanishing of those the paper sets to zero; and all 493 dual column inequalities satisfied in exact rationals, with dual objective exactly `46/13` |
| 8. the bounds meet | 6 | Theorem 1 together with each published number it must respect: `7/2 < 46/13` with margin `1/26`; `46/13 < 4`; `gamma_f(C_6 [] K_2) = 3`; the corollary of [1], `46/13 - 3 = 7/13 < 1`; and `46/13 < gamma = 4` |
| 9. the remark | 8 | the remark closing Section 3: the `x = 0` configuration totals `7/2` and is fractionally dominating; the printed response uses only edges and self-loops, respects each vertex's one-hop budget, conserves `7/2`, puts exactly 1 on the attacked vertex `j`, and leaves all twelve closed-neighbourhood sums at least 1, with `h = l = 1/4` nonzero |
| 10. an ancillary bound | 7 | `gamma_f^infty(Cay(Z_16, {+-1, 8})) <= 9/2`, from a seed, its 16 translates and all 256 ordered pairs, with `gamma = 5`. The paper asserts no such bound (see §3) |

## 3. What the program does not check

Taken from Sections 4 and 5 of the paper and from the closing `NOTE SCOPE` line of
`verify.output.txt`; intended to be complete.

* **It cannot decide which game is being played.** Every check encodes the same
  closed-neighbourhood reading, in the same closed neighbourhoods, so no check can discriminate
  between that reading and a stricter one; for a stricter rule nothing in the run bounds
  `gamma_f^infty`. The `NOTE SCOPE` line says exactly this, as do Section 1 and Section 4 (1).
* **The main theorem is a hand proof and the program is a control on its certificates.** Lemma 2
  (a reconfigurable family gives an eternal strategy), Lemma 3 (a finite attack tree is a
  relaxation), their proofs, and the weak-duality step of Proposition 4 are proved by hand and
  are not checked at all. What the program checks is the finite data those reductions consume:
  conditions (i), (ii), (iii) of Lemma 2 for the 12 members and the 144 ordered pairs, and dual
  feasibility together with the dual objective for the LP of Proposition 4.
* **Only `>=` is certified for the LP optimum.** Proposition 4's assertion that the optimum of
  (1) for this tree is exactly `46/13` combines the dual certificate with the upper bound of
  Section 2 and Lemma 3; no primal optimum of (1) is exhibited or verified.
* **No search and no solver is re-run.** Section 4 (5) states that a floating-point linear
  programming solver located the witness family and the dual vector; the `NOTE SCOPE` line lists
  exactly those floating-point programs as `NOT RE-RUN`, the objects being re-checked only as
  exhibited rationals. That search is not part of this folder, and the paper does not rest on it:
  every object above is exhibited as exact rationals and every verification is exact.
* **Only the 9-node tree is certified.** The `NOTE SCOPE` line states that no attack tree other
  than the 9-node tree printed in the paper is built or certified. Section 4 (4) adds that other
  trees were solved, that no optimum of theirs is quoted or used and no certificate for any of
  them is printed, and that Theorem 1 rests on the printed 9-node dual and on nothing else.
* **Nothing bears on any other cell, or on general `k`.** The `NOTE SCOPE` line: nothing bears on
  any cell of Problem 7.5 other than `C_6 [] K_2`, and no general-`k` statement is tested.
  Section 4 (2) lists the cells left open under this reading, `{C_{4k+2} [] K_2 : k >= 3}`
  together with `{Cay(Z_{8k}, {+-1, 4k}) : k >= 2}`; Section 4 (3) adds that the method gives no
  construction for general `k` and does not address the source's separate questions of whether
  `gamma_f^infty` is always rational, or whether `n` fractional dominating functions always
  suffice.
* **What is attributed to [1] is transcribed, not recomputed.** The bounds
  `7/2 < gamma_f^infty <= 4`, the values `8/3` and `28/5`, the estimate
  `(n + kappa)/(kappa + 1)`, the corollary `gamma_f^infty - gamma_f < 1`, and in the closing
  remark both the quoted sentence and the vertex labelling `a, ..., f`, `g, ..., l` of the
  relevant figure of [1], are taken from the cited source. Step 2 does recompute `8/3` from the
  source's own printed witness and Step 8 checks arithmetic consistency with the published
  numbers, but no check can confirm that [1] says what the paper reports it says.
* **The run is wider than the paper, and the paper depends on none of the excess.** The closing
  paragraph of Section 5 says so and lists it: controls, a second family of total weight `46/13`,
  three further transportation plans, an ancillary upper bound for `Cay(Z_16, {+-1, 8})`, and a
  closing scope note whose phrase "the accompanying note" is not to be read as a scope statement
  about the paper's sections (in `verify.py` and `verify.output.txt`, "the note" means this
  paper). Concretely: Steps 2, 4 and 10, three of the four plans of Step 5, and the second-family
  half of Step 6 concern objects the paper never states, and `gamma(C_6 [] K_2) = 4` and
  `gamma(Cay(Z_16, {+-1, 8})) = 5` occur only in the run.

## 4. How to check it

By hand, from the paper alone: the seed's total weight and its twelve closed-neighbourhood sums;
the 13-arc plan of Section 2, whose data are integral after scaling by 13; the dual objective
`17/13 + 29/13`. The 493 dual column inequalities are, as Section 5 says, tedious rather than
difficult, and are what the program is for.

By program, in this folder:

```sh
python3 verify.py
```

Standard library only, no arguments and no input files. It prints one line per check and a
closing verdict, and exits 0 only if every check passes. The recorded run used Python 3.9.25 and
exited 0.

The first lines of `verify.output.txt` are a provenance header, which that file itself marks as
not written by the program; they carry the SHA-256 of `verify.py` beside its output, so program
and transcript can be paired:

```sh
shasum -a 256 verify.py
```

    d401441a5b805f1d727a5b9ff6bd7b7e5bd93a6e5f7fe2a97f27c30db307a145

That digest was computed from the `verify.py` in this folder and equals the one in the header of
`verify.output.txt`.
