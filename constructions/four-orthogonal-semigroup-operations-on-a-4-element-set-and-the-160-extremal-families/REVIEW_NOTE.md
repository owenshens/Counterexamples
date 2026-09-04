# Referee note

*Four Pairwise Orthogonal Semigroup Operations on a Four-Element Set, and the 160 Extremal Families*

Files here: `paper.tex`, `paper.pdf`, `verify.py` (a checking program), `verify.output.txt` (one
recorded run of it), and this note. There is no data file; the program's docstring says its input is
the objects the paper prints and "Nothing else: no file, no artifact, no data."

## 1. What the paper claims

Araújo, Bentz, Cameron, Hendrey and Kinyon (arXiv:2608.25092v1) close their paper on complete
mappings of semigroups with thirteen problem environments. The seventh asks, for the set S_n of
semigroup operations on a fixed *n*-element set, for max(*n*), the maximum size of a pairwise
orthogonal family, and for the extremal families, as functions of *n*; Section 1 quotes it character
for character (lines 4606–4611 of the e-print). The paper settles the cell *n* = 4.

* **Theorem 1.** max(4) = 4: a pairwise orthogonal family of four operations on {1,2,3,4} is
  exhibited in Section 2, no five semigroup operations on a 4-element set are pairwise orthogonal,
  and there are exactly 160 labelled extremal families, in 10 orbits under the diagonal action of
  Sym(4) of sizes 8, 8, 8, 8, 8, 24, 24, 24, 24, 24, listed in Section 4.
* **Proposition 2.** max(1) = 1, max(2) = max(3) = 3, with 1, 2 and 3 labelled extremal families, in
  one orbit each.

The division of labour matters. **max(4) ≥ 4 is a hand proof.** Section 2 identifies X with (F_2)²,
and **Lemma 3** proves, for an arbitrary elementary abelian 2-group A, that f_M(a,b) = b + M(a) is
associative iff M² = M, that f_M ⊥ f_N iff M + N is invertible, and that f_M ⊥ L for L(a,b) = a.
Applied to the three rank-1 idempotents M_2, M_3, M_4 printed there this gives all six
orthogonalities "with no computation involved" (Section 2), which adds that Section 3 "supplies only
the matching upper bound max(4) ≤ 4". **Lemma 4** (Section 3) collects the elementary facts —
balance 4(i), commutative operations never orthogonal 4(ii), max(*n*) ≥ 3 by 4(iii), 4(iv), and the
classical ceiling 4(v) max(*n*) ≤ N(*n*) + 2, giving only max(4) ≤ 5. The one statement left to a
computation is therefore the negative *no five associative operations on a fixed 4-set are pairwise
orthogonal*, obtained by exhaustion over all 3492 labelled associative tables on [4].

Section 2 also states that the fourth exhibited table f_4 is the Cayley table published at lines
706–720 of the same e-print (the Mace4 example, with E(S) = {1,3}) and is cited, not claimed.
Lemma 3(ii) and Lemma 4(v) are attributed to the classical literature in Section 5, (C) and (A).

## 2. What the program checks

`verify.output.txt` ends `VERDICT: ALL 155 CHECKS PASS` and `=== program exited with status 0 ===`.
The header above the program output records the program name `verify.py`, its SHA-256 and
`Python 3.9.25`; the banner states standard library only and exact integer arithmetic. The 155 checks
fall in five labelled blocks (34 + 26 + 39 + 19 + 37).

* **A, "the family exhibited in Section 2" (34).** For f_1…f_4: array form agrees with the printed
  flat row-concatenated string; associativity over all 64 triples; value-balance (each value in 4 of
  the 16 cells); all six pairwise orthogonalities; the four self-pairs *not* orthogonal;
  Rees–Sushkevich type (4,1,1) for f_1 and (1,2,2) for f_2, f_3, f_4; f_4 identical to the source's
  published table, with idempotents {1,3}. This is the lower half of Theorem 1 as one instance, plus
  the attribution paragraph.
* **B, the F_2 picture (26).** The xor forms f_1(a,b) = a and f_i(a,b) = b + M_i(a), with M_2, M_3,
  M_4 sending 0,1,2,3 to (0,0,2,2), (0,1,1,0), (0,3,0,3) as printed; each idempotent of rank 1; the
  three sums of determinant 1 over GF(2). Then Lemma 3 exhaustively at this size: associative for
  exactly the 8 idempotent M of the 16 linear maps, orthogonal iff M + N invertible over all 256
  ordered pairs, the left-zero band orthogonal to every linear f_M, exactly six rank-1 idempotents
  whose orthogonality graph is 3-regular with 9 edges, two triangles and no 4-clique (M_2, M_3, M_4
  being one triangle), and M = 0, N = I to show higher-rank idempotents are not excluded. Three MOLS
  of order 4 are exhibited so the ceiling N(4) + 2 = 5 is not vacuous, only one of them associative.
  (The run heads this block "the F_2 picture of Section 3"; the F_2 material is Lemma 3, in Section 2.)
* **C, the unpruned census (39).** Every associative table on a fixed *n*-set is generated for
  *n* = 1, 2, 3, 4 with no lemma and no symmetry reduction: counts 1, 8, 113, 3492 against A023814,
  diagonal orbit counts 1, 5, 24, 188 against A027851, associativity re-tested, no repeats, and at
  *n* = 2, 3 the generator cross-checked against a sweep of all 2⁴ and all 3⁹ tables. Balanced counts
  1, 4, 5, 48, with an explicit over-full-value certificate for every unbalanced table (3444 at
  *n* = 4), so Lemma 4(i) is instantiated rather than assumed; value-balanced ⇔ completely simple as a
  set identity in both directions over each whole census, which is the unnumbered Remark of Section 3;
  edge counts 0, 5, 7 and 354 of C(48,2) = 1128; and a control of 100000 seeded random pairs with an
  unbalanced member, none orthogonal.
* **C2, the cliques (19).** Clique counts by size at *n* = 4, {1: 48, 2: 354, 3: 542, 4: 160, 5: 0},
  giving max(4) = 4 with 160 labelled extremal families, and at *n* = 2 ({1: 4, 2: 5, 3: 2, 4: 0})
  and *n* = 3 ({1: 5, 2: 7, 3: 3, 4: 0}), each of maximum 3 with 2 and 3 extremal families — that is
  Proposition 2. The last two *n* = 4 entries are obtained a second time by an algorithm-free scan of
  all 194580 four-subsets and all 1712304 five-subsets pair by pair, agreeing on the same 160 families
  and on the empty 5-level, so the load-bearing negative rests on no clique algorithm. Also the
  Rees-type breakdown of the 48, and that the 16 group tables on [4] (12 labelled Z_4, 4 labelled
  Klein) span no edge, as Lemma 4(ii) says.
* **D, the classification (37).** The second half of Theorem 1 and the table of Section 4: 10 orbits
  of sizes [8, 8, 8, 8, 8, 24, 24, 24, 24, 24] summing to and covering the 160; each printed
  representative an extremal family of its printed size, the ten meeting 10 distinct orbits; the
  Section 2 family is the sixth representative (the row marked (∗)); L ⊥ R, and the 16 / 16 / 128 / 0
  split with none containing both, as Lemma 4(iv) forces; exactly 8 of the 160 contain a group
  operation, those being the four Klein tables in one orbit of 8, including the printed
  `1234214334124321`; and {L, R, Klein} is a 3-clique in no maximum family.

## 3. What the program does not check

The run closes with a SCOPE block, items (a)–(d), carried over here.

1. **The lower bound is a hand proof and the program is a control on it.** Lemma 3 is proved for an
   arbitrary elementary abelian 2-group A; the program tests it only at |A| = 4 — all 16 linear maps,
   all 256 ordered pairs on (F_2)². The quantifier over A is proved, not checked.
2. **Completeness of the vertex set is transcribed from a cited source, not proved.** That the 3492
   generated tables are all of S_4 is pinned by agreement with the published A023814(4), and the
   orbit count 188 with A027851(4). Section 3 says so and adds that this "cannot be eliminated by any
   argument inside this note".
3. **The ceiling is classical and is not reproved.** SCOPE (b): three MOLS of order 4 are exhibited
   only to show N(4) + 2 = 5 is not vacuous; N(4) = 3 is not proved, nor N(6) = 1 (Tarry). Section 6
   repeats that Lemma 4(v) is quoted from the classical literature.
4. **Nothing at *n* ≥ 5.** SCOPE (a): the cells *n* = 1, 2, 3, 4 and nothing beyond; no table of
   order above 4 is generated anywhere, so no value of max(*n*) for *n* ≥ 5 is computed or supported.
   Section 6, "What is not settled", says the result should not be quoted as closing the Problem,
   which asks for max(*n*) as a function of *n*.
5. **The attributions of Section 5 are textual and bibliographic.** SCOPE (c): facts about
   arXiv:2608.25092v1 and about the Belousov school, covered by no check. The identification of f_4
   with the table at lines 706–720 is checked only at the level of the table itself and of
   E(S) = {1,3} — not the quotation, not the line numbers; likewise for the Problem quoted in
   Section 1 and its line numbers.
6. **Novelty and minimality are not checked.** SCOPE (d) records them as judgements.

No computation outside this folder is relied on: `verify.py` regenerates from the printed objects
alone every count the paper reports for *n* ≤ 4, and the only external quantities are the two cited
sequence values of point 2, flagged in the paper at the point of use.

## 4. How to check it

```sh
shasum -a 256 verify.py     # 097444922d37d99ea191c1f7a489e071fcafae6863412ad36633fa4b6f4e6560
python3 verify.py           # Python 3.9+, standard library only; exit status 0 iff every check passes
```

The digest was computed from `verify.py` as shipped here, and it is the value carried in the header of
`verify.output.txt` beside the program name, so transcript and program can be paired before the
transcript is trusted. The program takes no arguments and reads no input. It is short: its NOTE lines
report the *n* = 4 census as 3492 tables from 136152 backtracking nodes, and 1628 ms elapsed in total.
