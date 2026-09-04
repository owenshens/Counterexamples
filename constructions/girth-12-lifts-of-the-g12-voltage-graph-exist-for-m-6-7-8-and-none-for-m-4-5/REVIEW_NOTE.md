# Referee note

**Girth-12 Lifts of the Voltage Graph $G_{12}$ for $m = 6, 7, 8$, and Nonexistence on Its Printed
Arc Set for $m = 4, 5$**

The folder holds four files: `paper.tex` and `paper.pdf` (the paper), `verify.py` (the program of the
paper's Section 7), and `verify.output.txt` (a recorded run of it). Nothing else is required: the
program is Python 3.9, standard library only, no solver and no external data file, and the recorded
run states that the program does not parse the paper.

## 1. What the paper claims

Aguilar, Araujo-Pardo and Berman (reference [1], *Ars Math. Contemp.* **26** (2026) #P3.03; preprint
arXiv:2305.03290v1 as reference [2]) construct in their Theorem 5.5 a voltage graph $G_{12}$ whose
$\mathbb{Z}_m$-lift is a $(\{3,m\};12)$-graph of order $49m+3$ for every $m \ge 9$, and ask, in the
**second sub-bullet** under "Related to semicubic graphs" of their Section 6 "Open questions", for
voltage assignments producing $(G_{12};m)$-graphs for $4 \le m \le 8$, or a proof that none exist.
That sub-bullet is the target; the paper's Section 1 quotes it and Theorem 5.5 verbatim.

* **Theorem 1 (existence).** For each $m \in \{6,7,8\}$ there is a $\mathbb{Z}_m$ voltage assignment
  on the $24$ arcs of the paper's Table 2 whose lift is a $(\{3,m\};12)$-graph of order $49m+3$,
  namely $297$, $346$, $395$. One assignment per $m$ is printed in Section 4 as $24$ residues in the
  symbol order $a,b,c,\dots,\varepsilon$ of Table 2 (the $m=8$ vector is the $m=7$ vector read in
  $\mathbb{Z}_8$ with $g = 7$ in place of $6$). This answers the sub-bullet at those three values.
* **Theorem 2 (nonexistence, narrow).** With the $24$ arcs of Table 2 held fixed as printed, so that
  only the $24$ voltages vary, no $\mathbb{Z}_m$ voltage assignment has a lift of girth $12$ at
  $m = 4$ (**Proposition 6**, a pigeonhole on four arc voltages, no computation) nor at $m = 5$
  (Section 6, **Proposition 7**, by exhaustion of all $5^{24}$ assignments).

Supporting: **Lemma 3**, the criterion that the lift has girth $12$ iff no simple base cycle of length
at most $10$ has voltage sum $\equiv 0 \pmod m$; **Proposition 4**, that this cycle set has $126$
members, $12$ of length $6$, $24$ of length $8$, $90$ of length $10$, matching the source's published
"252 short cycles" counted in both orientations; **Proposition 5**, that the $24$ arc voltages
parametrise the $\mathbb{Z}_m$ lifts up to isomorphism, so the cells are exactly $\mathbb{Z}_4^{24}$
and $\mathbb{Z}_5^{24}$. The Remark in Section 2 argues that the eighth "Starting leaf" entry of the
source's Table 3 is a misprint for $x_{1011}$; the paper flags that correction as its own (Section 8,
item 3), and all its numbers use the corrected reading.

**No record and no bound is claimed.** Section 1 records the published upper bounds
$n(\{3,6\};12) \le 243$, $n(\{3,7\};12) \le 334$, $n(\{3,8\};12) \le 374$, quoted at second hand from
the girth-$12$ rows of the bound table of arXiv:2411.17351v1 (reference [3]), and states that the
witnesses, of orders $297$, $346$, $395$, are dominated in every case.

## 2. What the program checks

`verify.output.txt` records **80 checks, all passing**, closing `VERDICT: ALL 80 CHECKS PASS`, plus
three `NOTE` lines that are not checks. By block, with the claim each block serves:

* **19 checks — the base object of Section 2.** The $51$ tree edges built twice and shown equal
  (transcribed Table 1 against the source's bit-string pruning rule); $|V(T_{12})| = 52$,
  $|E(T_{12})| = 51$, $24$ arcs, $|E(G_{12})| = 75$, cycle rank $24$; degree multiset
  $\{1^3,3^{49}\}$ with the three degree-$1$ vertices exactly the pinned $x^*,y^*,z^*$; the arc
  incidences at the $x$-, $y$- and $z$-leaves; connectedness and bipartiteness; both printed voltage
  rows reproduced by the symbol order, largest printed absolute voltage $3$; the three
  pinned-to-pinned tree distances all $6$; the four tree distances $9$ that Proposition 6 uses; the
  count identity $48 + 3 = 51$ of Proposition 5; and the impossibility of the *literal* printed
  Table 3, under which $x_{1001}$ would have degree $5$ — one of the three grounds given in the
  Remark of Section 2.
* **9 checks — Proposition 4 and the bookkeeping of Lemma 3's criterion.** The census $126$ with
  length histogram $\{6\!:\!12,\ 8\!:\!24,\ 10\!:\!90\}$; a second census by XOR of fundamental
  cycles over the co-tree arcs, agreeing on edge sets; $2 \times 126 = 252$; no $4$-cycle; every
  short cycle carrying at least one arc; the support-size histogram of the $126$ forms; and the split
  of Section 3 into $14$ arcs forced nonzero and ten whose fundamental cycles have length $12$.
* **2 controls on objects outside the problem.** The Heawood graph as a $\mathbb{Z}_7$ lift (cubic,
  order $14$, girth $6$ by both girth routines), and the arc-free lift of $T_{12}$ at $m = 4,5,9$, on
  which the check "must stay SILENT here, and does".
* **8 checks — the source's Theorem 5.5 and the profile quoted in Section 1.** Order $49m+3$ and
  size $75m$ for $m = 3,\dots,14$; three vertices of degree $m$ and $49m$ of degree $3$; girth
  $6,6,8,8,8,10$ at $m = 3,\dots,8$ and $12$ at $m = 9,\dots,14$ under the printed voltages, so the
  source's threshold $m \ge 9$ is exactly attained; girth $8$ and not $12$ at $m = 6$ there, so
  Theorem 1's assignments are genuinely different ones on the same arc set; and three checks on that
  assignment's cycle sums, including the source's set $\{1,\dots,8\}$ attained exactly.
* **24 checks — Theorem 1.** For six vectors, one per $m \in \{6,7,8\}$ in each of two families:
  residues already lie in $\mathbb{Z}_m$; the full lift has $49m+3$ vertices ($297$, $346$, $395$),
  $75m$ edges, three vertices of degree $m$ and $49m$ of degree $3$; girth exactly $12$ **by two
  algorithms** (BFS from every vertex, and deleting each edge then BFS between its endpoints); and
  none of the $126$ forms vanishes mod $m$, checked independently of the lift. The transcript's
  "first family" is the three vectors printed in Section 4; the second family is supplied by the
  program and appears nowhere in the paper.
* **10 checks — Proposition 6.** Each of $k,p,v,\alpha$ closes a fundamental cycle of length $10$ and
  so is forced nonzero; the six difference forms arise from cycles of lengths $6,6,8,8,10,10$ as
  tabulated in the proof; the system is exactly ten forms, each one of the $126$ enumerated cycle
  forms; it is UNSAT over all $256$ four-tuples in $\mathbb{Z}_4$ and also at $m = 3$; it is
  **vacuous** at $m = 5,6,7,8$, so it cannot have killed those cells, and all six witnesses satisfy
  it; and a second independent core $b,d,f,h$ is likewise ten forms, UNSAT at $m = 4$, vacuous at
  $m = 5$.
* **8 checks — Proposition 7 and the exhaustions.** UNSAT over the complete cells $\mathbb{Z}_4^{24}$
  and $\mathbb{Z}_5^{24}$, all $24$ variables branched, all $126$ forms enforced, no gauge quotient
  and no unit-orbit cut; the same search returning solutions at $m = 6,7,8$ whose full lifts have
  girth $12$, a control that the search can find a solution where one exists; and a subset search
  finding that at $m = 4$ the only UNSAT cores on at most four variables are the three $4$-cliques
  $\{b,d,f,h\}$, $\{j,k,l,p\}$, $\{k,p,v,\alpha\}$, none on three or fewer, while at $m = 5$ there is
  none on four or fewer. Two `NOTE` lines report the search sizes — $184$ nodes at $m = 4$, $66505$
  at $m = 5$ — expressly as "reported, not claimed", node counts depending on the variable ordering.

Only the *necessity* half of Lemma 3 is used for the two negatives, so an incomplete constraint set
could not manufacture them; Theorem 1 uses Lemma 3 not at all, each witness being verified on its full
lift.

## 3. What the program does not check

* **Proposition 6 is a hand proof and the program is a control.** The $m = 4$ result is a
  four-variable pigeonhole to be read. The program re-derives its ten constraints, confirms each is
  one of the enumerated short cycles, and confirms the $256$-case UNSAT of that subsystem; it does not
  verify the argument. The separate full-cell $\mathbb{Z}_4^{24}$ exhaustion confirms the conclusion,
  not the argument.
* **Proposition 7 has no hand proof.** Section 6 says so: "This case is machine-only here: we know of
  no short hand argument for it." A referee must accept a branch-and-prune exhaustion of
  $5^{24} = 59{,}604{,}644{,}775{,}390{,}625$ points against the $126$ necessary conditions; that
  search is the only evidence offered.
* **Lemma 3 is not verified as a statement.** Its sufficiency half rests on a verbatim quotation of
  the source's proof of Theorem 5.5 (lollipop walks; non-reversing non-cycle walks) plus the paper's
  own lollipop count $2d + c + 2 \ge 12$; only the ingredients — no $4$-cycles, the pinned-to-pinned
  distances $6$ — are computed. Neither theorem depends on sufficiency.
* **Proposition 5's normalisation is not re-verified, only its counts** $49$, $48$ and
  $48 + 3 = 51$. The run says exactly that in the text of the check itself, naming the Gross–Tucker
  voltage-graph normalisation as the cited argument. Theorem 2 does not depend on the "up to
  isomorphism" clause: both negatives quantify over the whole cell with no quotient taken.
* **Quoted material and the published bounds are transcriptions, not computations.** The source's
  Theorem 5.5, the sub-bullet, the "252 short cycles" sentence, the commented-out preprint material
  that is one of the three grounds for the Table-3 correction, and the bounds $243$, $334$, $374$ are
  read out of references [1], [2], [3]; those sources are not in this folder and cannot be checked
  from it. The *derived* quantities — $252 = 126 \times 2$, the set $\{1,\dots,8\}$, the order
  $49m+3$, the girth profile — are recomputed.
* **The wide reading is untouched at $m = 4$ and $m = 5$.** The paper does not settle whether the
  sub-bullet fixes the printed arc set. Its "Scope of Theorem 2" paragraph and Section 8, item 1,
  require Theorem 2 always to be read as "no $\mathbb{Z}_m$ voltage assignment on the arc set of
  $G_{12}$ gives girth 12" and never as the unqualified "no such graphs exist", and state that under
  the wider reading, in which each $x$-leaf may be re-paired to some $y$-leaf and some $z$-leaf, the
  cells $m = 4,5$ remain open. The issue does not arise at $m = 6,7,8$. The group is $\mathbb{Z}_m$
  and nothing else, so the $m = 4$ cell is $\mathbb{Z}_4^{24}$ and not also
  $(\mathbb{Z}_2 \times \mathbb{Z}_2)^{24}$.
* **The run states its own scope.** Its closing `NOTE SCOPE` line lists what it re-derives, says it
  does not parse the paper and "makes no claim to have enumerated everything the paper states", and
  lists as `NOT RE-RUN HERE`: the wide variant in which the $24$ arcs themselves are re-paired, which
  is untouched at $m = 4,5$; the gauge/normalisation argument, of which only the counts are checked;
  and any statement about published order records, of which the paper makes none.
* **The transcript contains more than the paper claims.** Three of the six witness vectors are the
  second family, which the paper neither prints nor uses, and the small-core searches support no
  printed statement; Section 7 flags this ("The recorded run also prints checks on matters this note
  makes no claim about"). The $m = 5$ core search covers subsets of at most **four** of the $24$
  variables and, in its own words, "says nothing about larger subsets".
* Nothing is claimed about $m \ge 9$ beyond reproducing Theorem 5.5, nor about the sibling $H_{10}$
  and $H_{12}$ sub-bullets of the same open-questions section, nor about the $(3,14)$-cage bullet
  (Section 8, item 2).

## 4. How to check it

```sh
python3 verify.py          # one PASS line per check, then the verdict; exits 0 iff all pass
shasum -a 256 verify.py
```

    8417a0874d45ee32e3a1bbd4a5dccce1fcb7145b1427d63fcb583a9e3fb47d36

The transcript opens with a provenance block carrying that same SHA-256 beside the program name, and
the Python version of the recorded run (3.9.25), so program and transcript can be paired. That block
and the closing status line are marked in the file itself as not being the program's own output;
everything between them is.

## 5. One editorial defect

Section 2 and Section 8, item 3, both refer to "Remark 2", but the paper's Remark environment is
unnumbered and prints simply as "Remark"; the "2" is the number of the preceding Table 2. Both should
read "the Remark in Section 2". No mathematical claim is affected.
