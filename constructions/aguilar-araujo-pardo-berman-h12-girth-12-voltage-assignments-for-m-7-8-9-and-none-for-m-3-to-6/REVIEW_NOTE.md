# Review note

Paper: *Voltage assignments on a 44-vertex reconstruction of $H_{12}$: girth 12 for $m=7,8,9$, and
none for $3\le m\le 6$.*

Folder: `paper.tex`, `paper.pdf`, the program `verify.py`, and its recorded run `verify.output.txt`.
The base graph, all six voltage 20-tuples and the three exhibited 12-cycles are printed in full in
the paper, so nothing outside this folder is needed to re-check them. The one thing that does need the
cited source is the *identification* discussed under "What the program does not check" below, and the
paper says so (its Sections 7 and 9).

## 1. What the paper claims

Aguilar, Araujo-Pardo and Berman (*Semicubic cages and small graphs of even girth derived from
voltage graphs*, Ars Math. Contemp. **26** (2026) #P3.03) ask, in the first sentence of the third
sub-bullet under "Related to semicubic graphs" of their Section 6 "Open questions", p. 29 — quoted
verbatim in Section 1 — for voltage assignments on their voltage graph $H_{12}$ giving girth 12 for
$m\in\{3,4,\dots,9\}$, or a proof that none exists.

For **one explicit 44-vertex voltage graph**, printed completely in Section 2 (43 tree edges at
voltage 0; 20 labelled arcs in the fixed coordinate order
$[a,b,c,d,e,j,q,u,w,s,l,g,f,k,r,v,\alpha,t,p,h]$), the paper proves that girth-12 lifts exist for
$m=7,8,9$, exhibited at orders **290, 331, 372** (**Theorem 3**), and that none exists for
$m=3,4,5,6$ (**Theorem 4**) — with **Proposition 1** (the base combinatorics: 63 edges, cycle rank
$63-44+1=20$, bipartite with classes 23 and 21, hence lifts of order $41m+3$ with $63m$ edges, three
vertices of degree $m$ and $41m$ of degree 3), **Lemma 2** (girth $\ge12$ forces a nonzero voltage sum
on every base cycle of length $<12$, over the space $\mathbb{Z}_m^{20}$ that Section 3 argues is
complete up to switching), **Corollary 5** (the $(3,12)$-cage is therefore not a $\mathbb{Z}_3$ lift of
that graph) and **Corollary 6** ($n(\{3,7\};12)\le290$ and $n(\{3,8\};12)\le331$). A subsection of
Section 4 prints three further witnesses that additionally keep the source's Table-4 values
$a,b,c,d=1,-1,-1,1$; Theorem 3 itself is unrestricted.

Corollary 6 is what settles something in the literature: per Section 6, the $m=8$ value 374 is
Corollary 2.4 of the source and 331 improves it, while the $m=7$ value 334 belongs to Araujo-Pardo,
Balbuena, López-Chávez and Montejano, Theorem 2.2, so nothing here improves the source's authors at
$m=7$.

**The paper does not claim to settle the cited question, and says so repeatedly.** Section 2 fixes the
convention that "$H_{12}$" means the graph printed there "and nothing else; it is not asserted to be
the graph of [Figure 16]". Section 7 records the two places where the source's Tables 3–4 must be read
against p. 27's own deletion list, states that no arc-by-arc transcription of Figure 16 is supplied,
and concludes that Proposition 1, Theorems 3 and 4 and Corollaries 5 and 6 are statements about the
printed graph — Corollary 6 being the one consequence unaffected, its two witnesses being checked as
graphs in their own right. Section 8's first bullet repeats that the question "is not claimed to be
answered, in either direction". Three further narrowings, all the paper's own: the $m=9$ witness
carries **no** record claim (372 realises the order the source names on its p. 28, but Goedgebeur,
Jooken and Van den Eede publish $n(\{3,9\};12)\le360$); the sub-bullet's second question, a *different*
graph $H'_{12}$, is untouched; and Table 1 of the source prints $\mathrm{rec}(\{3,8\};12)=304<331$,
which Section 6 reads as a misprint for 374 while noting that no erratum is verified.

## 2. What the program checks

`verify.output.txt` closes with

    VERDICT: ALL 117 CHECKS PASS
    === program exited with status 0 ===

Every object consumed is printed in the paper; no input file, no third-party package; exact integers
and `fractions.Fraction`. The 117 checks fall into five blocks.

* **Base graph, 15 checks** — Proposition 1 in full (44 vertices, 43 tree edges, 20 arcs, 63 edges, no
  loop or repeated edge, histogram $\{1{:}3,\,3{:}41\}$ with the degree-1 vertices exactly
  $x^{*},y^{*},z^{*}$, the tree spanning, cycle rank 20, bipartite with class sizes 21 and 23, every
  arc crossing); three of them cover Section 2's leaf restoration — the eight deleted leaves absent,
  $T_{12}$ on 52 vertices and 51 edges, 49 non-pinned, reproducing the source's $49m+3$ and
  $441=49\cdot9$.
* **Short-cycle census, 11 checks** — Section 3: shortest base cycle 6, no odd cycle below 12,
  $18+27+82=127$ undirected cycles of lengths 6, 8, 10, i.e. **254 directed**, reproducing the source's
  "There are 254 such cycles" (p. 27); every coefficient $\pm1$; the arity histogram; 16 arcs forced
  nonzero and the four unforced ones $w,s,p,h$.
* **The six lifts, 51 checks** — Theorem 3 and the restricted witnesses: seven per witness (order
  $41m+3$, $63m$ edges, no loop or multi-edge, degree histogram, connectivity, girth exactly 12, all
  127 necessary constraints satisfied) for the three at 290, 331, 372 and the three keeping $a,b,c,d$;
  three that each displayed 12-cycle is 12 distinct lift vertices with all 12 consecutive pairs edges
  of the lift; six that the restricted witnesses keep $a,b,c,d=1,-1,-1,1$ and differ from the headline
  ones. Girth is measured by breadth-first search from every vertex, depth capped at 6.
* **The four emptiness claims, 21 checks** — Theorem 4: per $m$, that the core enumeration is the full
  flat cell ($3^5=243$, $4^7=16384$, $5^8=390625$, $6^8=1679616$), that its constraints are exactly
  those of the 127 whose whole support lies in the core (10, 22, 37, 37), that there are **zero
  survivors**, and the order that cell would have had (126, 167, 208, 249); plus that the ten relations
  printed for $m=3$ are genuine constraints with no solution over $\mathbb{Z}_3$, that the four-case
  argument closes all four branches, and that the $m=3$ lift would be cubic on 126 vertices.
* **Published arithmetic, 19 checks** — Section 6 and Corollary 6: $\lceil109m/3\rceil+17$ at
  $m=3,7,8,9,10$; the three orders and the source's 413 for $(H_{12},10)$; $272<290<334$ (gap 18) and
  $308<331<374$ (gap 23); that 372 does not improve 360; Corollary 2.4 of the source in exact rationals,
  reproducing five Table 1 entries and giving 374 at $m=8$ against the printed 304, itself below the
  lower bound of its own row; that Table 1's lower bounds at $m=5,6$ likewise disagree with the source's
  own formula; and six touching Theorem 4.2 of Goedgebeur, Jooken and Van den Eede (their 360 at $m=9$
  and 243 at $m=6$, two single case values, and 290 and 331 against 362).

## 3. What the program does not check

* **The identification with the source's Figure 16 is not checked at all.** The program is handed the
  graph of Section 2 and checks statements about that graph; the readings of Section 7 are assertions
  about the source's Tables 3–4 and the prose on its p. 27, which are not reprinted, so a referee who
  wants them checked must fetch them and compare arc by arc (Section 9 says exactly this). Every
  headline statement is therefore conditional on that identification, except Corollary 6.
* **Three load-bearing arguments are hand proofs; the program controls their inputs only.** (i) The
  gauge-completeness argument of Section 3 — that normalising the tree to 0 makes $\mathbb{Z}_m^{20}$
  the whole space up to isomorphism of lifts — is a switching argument, of which the program checks
  only the two combinatorial inputs, that the 43 tree edges span and that the cycle rank is 20.
  (ii) Lemma 2 is a hand proof: the program computes the constraint system but does not establish the
  lemma. (iii) The depth-6 cap in the girth measurement is justified by the hand argument in the proof
  of Theorem 3; the program merely implements the capped search.
* **Given Lemma 2, Theorem 4 is decided in full, not sampled**: the four flat enumerations run over
  the printed cores with no search and no propagation, and no enumeration over $\mathbb{Z}_m^{20}$
  itself is performed or needed. But Theorem 4 does **not** say $n(\{3,m\};12)>41m+3$ for those $m$ —
  see the remark closing Section 5 — and, per Section 8, says nothing about the published
  $n(\{3,4\};12)\le220$, $n(\{3,5\};12)\le230$, $n(\{3,6\};12)\le243$.
* **Theorem 3 rests on the exhibited tuples, not on a search.** The paper states that the measurement
  of the lifts "is the only computation the existence half depends on, and it uses no solver". How the
  tuples were arrived at is not part of this folder, and the paper does not rest on it.
* **Transcribed, not recomputed.** The closing scope line of `verify.output.txt`, verbatim:

  > NOTE SCOPE: this program checks the objects and the arithmetic printed in the paper. NOT
  > RE-RUN here: the uniqueness of the (3,12)-cage, which the paper quotes from the literature and
  > flags as an inherited external fact; the existence of the s/2 = 4 edges at pairwise distance
  > >= 6 that case 3 of GJV Theorem 4.2 needs, nor the applicability conditions of its cases 1 and
  > 2, nor any minimisation over its parameters -- the two case-value checks above are single
  > evaluations, not floors; the numerical values 334, 374 and 360 of the published record table,
  > which are transcribed from the literature and not recomputed; and the total number of
  > girth-12 assignments at any m, which the paper does not claim. The four emptiness claims are
  > re-proved here in full, by flat enumeration over the printed cores with no search.

  So the uniqueness of the $(3,12)$-cage, on which **Corollary 5** depends, is inherited (flagged in
  the paper's remark after that corollary and in Section 8), and the record values 334, 374, 360 that
  make Corollary 6 an improvement are transcribed. The two Theorem 4.2 case values are single
  parameter choices, so the two "outside the reach" checks are relative to 362 and are not floors;
  Section 9 says the same. **No check verifies an erratum**, in particular none verifies that Table 1's
  printed 304 is a misprint for 374.
* **Neither claimed nor checked**: the number of girth-12 assignments at any $m$; optimality of 290
  and 331 (the lower bounds 272 and 308 are untouched, and no smaller $(\{3,7\};12)$ or
  $(\{3,8\};12)$ graph outside this family was sought); whether the witnesses of Section 4 are
  pairwise non-isomorphic or isomorphic to known graphs; the sub-bullet's second question
  ($H'_{12}$); and the *second* sub-bullet of the same page, on the companion family $G_{12}$ of
  order $49m+3$.

## 4. How to check it

```sh
python3 verify.py            # 117 PASS lines, then the verdict; exits 0
shasum -a 256 verify.py      # d3ccb747bc57fb11c0d83bba92c2b650270c4b50f04aaef10a016d8b0ca2ca9a
```

Python 3.9 or later, standard library only. The header of `verify.output.txt` carries the SHA-256 of
the program beside its name, so transcript and program can be paired; the digest above was computed
from the shipped `verify.py` and is the value that header prints. A rerun reproduces the recorded
output line for line, the only expected difference being the `python <version>` line the program prints
about itself; the recorded run used Python 3.9.25.

Two claims need no program at all: the existence half reduces to one breadth-first search on the
290-vertex lift, and the $m=3$ nonexistence is the ten-relation, four-case argument printed in
Section 5.
