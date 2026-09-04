# Review note — *A 25-Vertex Weak Tensor Product of α-Free Trees With an α-Valuation*

Besides this note the folder holds `paper.tex` and `paper.pdf` (the paper), `verify.py` (a checking
program) and `verify.output.txt` (the recorded run of that program). Nothing else is needed to read
the paper or to check it.

The paper has seven sections — §1 Introduction, §2 The factor, and its α-freeness, §3 The product,
and the α-valuation, §4 The product labelling σ, §5 Two distinct graphs, §6 Relation to known
results, §7 Scope — and numbers eight statements: Theorem 1, Fact 2, Proposition 3, Corollary 4,
Lemma 5, Corollary 6, Proposition 7, Theorem 8, plus one unnumbered Remark in §3.

## 1. What the paper claims

It answers **no** to the second of the four unnumbered items in the "Further Work" section of Anglin,
Huczynska and McCourt, *Graph labellings and external difference families* (arXiv:2603.05662v1),
which asks whether the weak tensor product of two graphs that each admit a near α-valuation but no
α-valuation must itself admit no α-valuation. §1 quotes that item verbatim and writes (∗) for the
universally quantified statement it invites.

**Theorem 1.** Let `S_{3,2}` be the 7-vertex spider (`K_{1,3}` with every edge subdivided once) with
the near α-valuation γ of equation (1), the one printed in Anglin–Huczynska–McCourt. Then `S_{3,2}`
has no α-valuation (Proposition 3), while `K = S_{3,2} ⊗̄ S_{3,2}`, taken with respect to the
bipartitions induced by γ — connected bipartite, 25 vertices, 36 edges, parts of sizes 16 and 9 —
admits the α-valuation printed in §3, with boundary λ = 27. Corollary 4 concludes that α-freeness is
not preserved by ⊗̄. **Theorem 8** repeats this on the non-isomorphic pair `S_{3,2}`, `S_{4,2}`
(Proposition 7): a connected product `L` on 32 vertices and 48 edges, parts 20 and 12, with an
α-valuation of boundary λ = 36, so the question is not rescued by reading "two graphs" as two
distinct graphs.

§4 records why a positive answer was plausible. **Lemma 5** gives, for connected `G ⊗̄ H`, the
criterion `m(h_C − h_D + 1) < g_B − g_A` for the specific product labelling σ of El-Zanati, Kenig and
Vanden Eynden — the labelling on which the question rests — to be an α-valuation; **Corollary 6**
concludes that σ is never one when neither factor labelling is an α-valuation, so that construction
cannot itself produce a counterexample, and that when both factor labellings are α-valuations σ is
one, recovering Snevily's theorem. §1 and §7 are explicit that nothing published is contradicted:
Theorem 7 of El-Zanati–Kenig–Vanden Eynden and Snevily's theorem both stand.

## 2. What the program checks

`verify.output.txt` closes with `VERDICT: ALL 44 CHECKS PASS` and `program exited with status 0`; its
header records Python 3.9.25. The 44 checks fall in five parts.

* **Part 1 — 8 checks — §2, γ, Proposition 3.** `S_{3,2}`'s structure (|V|=7, |E|=6, tree, connected,
  degrees 3,2,2,2,1,1,1, bipartition `{c,l1,l2,l3}|{m1,m2,m3}`); γ is a β-valuation, is near α, and is
  **not** α (max γ(V_small)=4 > 3=min γ(V_large)); both case sums of Proposition 3 (`24−2b(c)` and
  `2b(c)+12`, all even, never 21); the 288 candidates permitted by Fact 2(ii) examined with 0
  α-valuations; and a census of all 7! = 5040 injections giving 60 β, 24 near-α, 0 α.
* **Part 2 — 15 checks — `K` and the witness `b` of §3 (Theorem 1).** `K` rebuilt from the ⊗̄
  definition quoted in §1 and matched edge-for-edge onto §3's hub/grid/row/column/pendant description
  under `h→cc, R_i→l_ic, C_j→cl_j, P_ij→l_il_j, Q_ij→m_im_j`; |V(K)|=25, |P|=16, |Q|=9, |E(K)|=36;
  connected, hence with a unique bipartition; all 36 edges joining P to Q; degree sequence
  9^1 4^9 3^6 1^9; `b` given once per vertex and injective in [0,36], the 12 unused labels being those
  §3 lists; the 36 edge differences exactly {1,…,36}; §3's four printed difference blocks reproduced
  value-for-value in the printed order, and those 36 printed values sorting to 1,…,36; 27 the **only**
  straddling x in {0..36} (max_P=27 < 28=min_Q); the Remark's two consistency checks (weighted degrees
  666 = 1+2+⋯+36, and difference 36 carried by `cl_2 = 0 ~ m_1m_2 = 36`, which are indeed adjacent);
  §4's non-separability observation (3 against 8); and a refusal check — swapping `b(cl_2)` with
  `b(l_1l_2)` is not graceful, and a non-injective labelling is refused outright.
* **Part 3 — 5 checks — Lemma 5 and Corollary 6.** σ is graceful on the same 36 edges of `K` but is
  **not** α there (max_P σ = 28 > 15 = min_Q σ); Lemma 5's criterion evaluates to 12 and −1, so
  12 < −1 is false, agreeing with the previous check; `P_3` labelled 0,2,1 is α with boundary 1; and
  two α factors give an α σ on the product (|V|=5, |E|=4, boundary 3), which is Snevily's theorem.
* **Part 4 — 13 checks — §5, Proposition 7, Theorem 8.** `S_{4,2}`'s structure (|V|=9, |E|=8, tree,
  degrees 4,2,2,2,2,1,1,1,1, not isomorphic to `S_{3,2}`); δ of equation (2) is β, is near α, is not α
  (5 = δ(b_2) exceeds 3 = δ(a_3)); both cases of Proposition 7 (the sums force b(c)=2 and b(c)=6, the
  centre edges give {3,4,5,6}, the leaf edges owe {1,2,7,8}, and differences 1 and 2 cannot coexist);
  the 5760 candidates permitted by Fact 2(ii) examined with 0 α-valuations; |V(L)|=32, |P|=20, |Q|=12,
  |E(L)|=48, connected, every edge P→Q; the degree profile of Theorem 8 (12 at `cc`, 3 at each `cb_j`,
  4 at each `b_ic`, 1 at each `b_ib_j`, totalling 48 edges); 32 distinct labels in [0,48]; the 48
  differences exactly {1,…,48}; 36 the only straddling x (max_P=36 < 37=min_Q); weighted degrees
  1176 = 1+2+⋯+48.
* **Part 5 — 3 checks — published controls in both polarities.** The zigzag α-labelling of paths for
  m = 2,…,8; Rosa's labelling of `K_{p,q}` for all 1 ≤ p,q ≤ 4, each with boundary p−1; and, by
  exhaustion over every injection, "`C_m` has an α-valuation precisely if m ≡ 0 mod 4" on 4 ≤ m ≤ 7 —
  yes for `C_4`, no for `C_5`, `C_6`, `C_7`.

Per its own header `verify.py` uses the standard library only, with no external data file and no
network, and exact integer arithmetic throughout, so no decision depends on rounding. Both products
are rebuilt by evaluating the quoted adjacency predicate on every pair of vertices, not read from a
stored edge list.

## 3. What the program does not check

Its own closing block, verbatim:

> SCOPE.  This program checks the objects PRINTED in the paper and the
>         finite claims made about them.  It does not search for further
>         counterexamples, does not claim K is of least order, and
>         establishes nothing about factors on more than 9 vertices.

No claim in the paper rests on any computation outside this folder. In more detail:

* **Every statement is proved by hand and the program is a control.** Propositions 3 and 7 are parity
  and counting arguments; Fact 2, Lemma 5 and Corollary 6 are general arguments. The decisive
  verification of each witness is a finite computation the paper carries out itself: §3 does `K` in
  four blocks with a sorted, block-tagged union of the 36 differences, and the proof of Theorem 8
  calls the corresponding step for `L` "a finite check on the printed labels". The program redoes that
  arithmetic independently; it is corroboration, not the proof.
* **Fact 2 is not itself checked, and two enumerations lean on it.** Fact 2 is general, part (ii)
  about connected bipartite graphs with |V| = q+1. Both α-freeness enumerations examine only the
  labellings Fact 2(ii) permits — 288 for `S_{3,2}`, 5760 for `S_{4,2}` — so each is exhaustive only
  granted Fact 2(ii). For `S_{3,2}` that dependence is removed by the unrestricted 5040-injection
  census; for `S_{4,2}` there is no such census, and its full set of injections is never enumerated.
* **Lemma 5 and Corollary 6 are quantified and only instantiated.** Lemma 5 is an if-and-only-if for
  every connected `G ⊗̄ H`, and Corollary 6 asserts its conclusion "for every such choice of γ and δ
  and in either coordinate order". The program evaluates the criterion on `K` alone, plus the single
  positive instance built from `P_3`. Lemma 5's proof cites Theorem 7 of El-Zanati–Kenig–Vanden Eynden
  for σ being a near α-valuation with those classes; that theorem is not reproved. §4 states that the
  disconnected case is **not** settled, and that no claim is made that a counterexample must be
  non-separable.
* **Some inputs are transcribed from the literature, not recomputed.** γ in equation (1) is read off
  the figure printed in Anglin–Huczynska–McCourt. The statement that `S_{3,2}` is a tree of smallest
  size with no α-valuation is quoted from that paper, attributed there to Rosa; no minimality over
  trees is verified anywhere, and §2 notes that the appeal to Rosa is not needed. §6's readings of the
  literature — that no weak-tensor-product theorem in Gallian's survey addresses α-free factors, and
  the accounts of Ahmed–Snevily, El-Zanati–Fu–Shiue and López–Muntaner-Batle — are outside the
  program's reach.
* **One residual is definitional rather than arithmetical.** §7 says so: that "near α-valuation" and
  "weak tensor product" are read as Anglin–Huczynska–McCourt intend. §7 records that every quotation
  is verbatim from that paper's LaTeX e-print with only formatting altered (its displayed formulae and
  its two-item list set inline, its citation and cross-reference macros rendered as bracketed
  references), and that the reading was pinned by reproducing that paper's own printed valuation,
  equation (1), and its stated properties — which the program does check. A definitional misreading is
  not something arithmetic can exclude, and auditing the quotations means consulting the cited
  e-print, which is not part of this folder.
* **Nothing is proved in the opposite direction.** §7 states that it is *not* claimed that every such
  product has an α-valuation, that `K` is not claimed to be the smallest counterexample, and that
  nothing is claimed about α-free near-α trees beyond the two factors used here, nor about non-tree
  α-free bipartite factors. §6 repeats that the two products exhibited are single instances.

## 4. How to check it

```sh
shasum -a 256 verify.py
python3 verify.py
```

The first prints

    88a1e82db678e47be5781a2ccdce1e93c370e9e0a2ef92ec7589a4baf7122e34

and the same value stands in the header of `verify.output.txt`, on the `sha256:` line beside
`program: verify.py`, so transcript and program can be paired. The second reproduces the 44 `PASS`
lines, the SCOPE block and the verdict; it needs Python 3.9 or later and nothing else (the recorded
run used Python 3.9.25). It exits 0 if and only if every check passes, printing `FAIL` lines and a
failing verdict otherwise; the recorded status is 0. In `verify.output.txt` only the header above
`=== program output follows ===` and the closing `=== program exited with status 0 ===` line are not
program output, as that file itself says.

No program is needed for the headline witness. §3 is a pencil check: 25 distinct labels in {0,…,36},
maximum 27 on the 16-vertex side against minimum 28 on the 9-vertex side, and the four listed blocks
of 36 differences sorting to 1,…,36.
