# Referee note — *A 24-Permutation Hamiltonian Cycle of P(4) That Does Not Extend to a u-Cycle*

Files a referee has here, and nothing else is needed: `paper.tex`, `paper.pdf`, `verify.py`,
`verify.output.txt`.

## 1. What the paper claims

$P(n)$ is the digraph of overlapping permutations: its vertices are the $n!$ permutations of
$\{1,\dots,n\}$, with an arc $p\to q$ exactly when $\operatorname{st}(p_2\cdots p_n)=\operatorname{st}(q_1\cdots
q_{n-1})$. Given a directed Hamiltonian cycle $C=(p_0,\dots,p_{N-1})$ of $P(n)$, $N=n!$, the
Chung–Diaconis–Graham procedure places undetermined letters $u_0,\dots,u_{N-1}$ on a circle and requires every
cyclic window $(u_i,\dots,u_{i+n-1})$ to be order-isomorphic to $p_i$; if the digraph $D(C)$ of forced strict
inequalities is acyclic, a linear extension of it yields a u-cycle for $n$-permutations.

**Theorem 2** prints one explicit 24-entry cyclic sequence (two tables in Section 2, closed by the arc
$3124\to1234$), proves it is a directed Hamiltonian cycle of $P(4)$ using neither self-loop of $P(4)$, and
exhibits inside $D(C)$ the directed cycle of strict inequalities

$$u_2<u_5<u_7<u_9<u_{12}<u_{14}<u_{17}<u_{20}<u_{23}<u_2,$$

displayed as equation (1). So that Hamiltonian cycle has no linear extension and extends to no u-cycle for
$4$-permutations.

What this settles is the universally quantified statement **as the two citing sources state it**: the belief
reported in `arXiv:2408.05984v1` for "*any* Hamiltonian cycle in $P(n)$", and the parenthetical "possibly all,
which is conjectured" of `arXiv:2603.01005v1`. **Lemma 1** ($P(n)$ is the line digraph of $H_n$, so
$\operatorname{Ham}(P(n))\leftrightarrow\operatorname{Eul}(H_n)$ is a bijection, not an inclusion) is what makes
the two quantifiers range over the same set, so the Eulerian route of the 2026 source is a reparametrisation of
the whole search space and not an escape. Section 1 states expressly that the 1992 paper was not read and that
no claim is made about what it states; Section 4 records that u-cycles for $n$-permutations exist for every $n$
and that Hurlbert's construction is untouched.

The proof is by hand: a 24-entry table, 24 arc checks, a first-letter/suffix tabulation, and nine readings of a
window. The abstract says it "needs no computation".

## 2. What the program checks

`verify.output.txt` ends `VERDICT: ALL 73 CHECKS PASS`. The 73 `PASS` lines fall into the transcript's four
steps; four `NOTE` lines carry no verdict.

- **Step 1 — 9 checks.** For each of two 24-entry objects (prefixes `w1_`, `w2_`): 24 entries; entries exactly
  the 24 permutations of $S_4$; all 24 arcs of $P(4)$ hold including the wrap arc $23\to0$; neither self-loop
  used. Plus one check that $P(4)$ has exactly two self-loops. The four `w1_` lines are the "Every arc holds"
  and "It is Hamiltonian" parts of the proof of Theorem 2.
- **Step 2 — 42 checks.** Twenty `w1_` lines: each of the nine links of (1) is forced by the window the paper's
  table names (each line also printing every window that forces the link); the nine links close into one
  directed cycle; the step sizes $[3,2,2,3,2,3,3,3,3]$ sum to 24; the forced digraph is cyclic, a topological
  sort placing only 7 of 24 positions; the consecutive-rank reading has 46 arcs and the literal
  all-ordered-pairs reading 72, the latter with zero antisymmetry violations, containing all nine links, and
  cyclic too; no forced arc spans cyclic distance more than 3; the shortest directed cycle has length exactly 9
  under both readings. These are the "contradiction" half of Theorem 2 together with the "Independence of the
  reading" paragraph of Section 3, where the counts 72 and 46 appear. Twenty further `w2_` lines repeat all of
  it for a second 24-permutation cycle with chain $u_0<u_3<u_5<u_8<u_{11}<u_{14}<u_{16}<u_{19}<u_{22}<u_0$
  (topological sort places 8 of 24). Two remaining lines record $\lceil 24/3\rceil=8$ and that the two objects
  are not rotations of one another.
- **Step 3 — 12 checks**, the controls of Section 3, pinned to the $n=3$ output published in
  `arXiv:2408.05984v1`: that source's cycle of $P(3)$ has all six arcs and is Hamiltonian on $S_3$; its
  published word $U_3=142342$, read as a circular word, gives exactly that cycle, *including* the two wrap
  windows $(4,2,1)$ and $(2,1,4)$; the forced order there is acyclic under both readings and $142342$ satisfies
  every forced strict inequality; the transitive reduction reproduces the published Hasse diagram; the word uses
  only 4 letters on 6 positions; $P(3)$ has exactly 8 directed Hamiltonian cycles and all 8 have acyclic forced
  order, $P(2)$ has 1 and it is good; an anti-control of three identical $123$ windows is correctly reported
  CYCLIC; and dropping the three wrap-around windows would make both exhibited objects look acyclic — the single
  failure mode Section 3 says the wrap-around convention exists to exclude.
- **Step 4 — 10 checks.** $|\operatorname{Ham}(P(4))| = 280\cdot 6^6 = 13\,063\,680$ from the BEST theorem with
  an exact integer Matrix–Tree cofactor: $H_4$ has 6 vertices, 24 edges, all in- and out-degrees 4; the
  $\operatorname{Ham}(P(4))$ and $\operatorname{Eul}(H_4)$ predicates agree on all $24\times24=576$ ordered
  pairs; the arborescence count 280 comes out the same from each of the six roots and again by brute-force
  enumeration with no linear algebra; the same code path reproduces the exhaustive count 8 at $n=3$ and the de
  Bruijn Eulerian-circuit counts 2, 16, 2048, 67108864.

The recorded run declares "exact integer arithmetic only"; `verify.py` imports only the Python standard library.

## 3. What the program does **not** check

- **Theorem 2 is a hand proof and the program is a control.** Section 5 opens by saying everything in
  Theorem 2 is checkable from the printed tables with no computer; the program re-derives that certificate, and
  the paper does not rest on it.
- **Lemma 1 is stated for all $n$ and proved by hand in Section 1.** The program tests the
  $\operatorname{Ham}(P(4))$/$\operatorname{Eul}(H_4)$ agreement only on the 576 ordered pairs at $n=4$.
- **The verbatim quotations are not checked.** They are byte-level facts about two external e-print sources;
  item (c) of the run's closing `SCOPE, NOT RE-RUN HERE` note says exactly this, and Section 5 repeats it.
- **Nothing at $n\ge5$** — item (b) of the same note; Section 4 says "we report no computation at $n\ge5$".
- **Whether every $n\ge4$ has such a cycle** — item (d). Only $n=4$ is settled (Section 4).
- **No census of the bad cycles at $n=4$** — item (a): counting the Hamiltonian cycles of $P(4)$ with cyclic
  forced order, or any refinement such as a shortest-chain histogram or a count of symmetry orbits, would
  require enumerating all $13\,063\,680$ cycles, and no such count is a claim of the program. Section 4 likewise
  claims no such count, and says the arc counts 72 and 46 are about the exhibited cycle alone. Relatedly,
  $\lceil 24/3\rceil = 8$ is checked only as arithmetic: the `NOTE` line after it states that whether 8 is
  attained by any Hamiltonian cycle of $P(4)$ with cyclic forced order is not re-derived.
- **Published data are reproduced by the program, not checked against it.** The header of `verify.py` names them
  as inputs: the $n=3$ Hamiltonian cycle, the word $U_3=142342$, the Hasse diagram, and the four de Bruijn
  counts. They are transcribed from the cited sources; being reproduced by an independent recomputation is what
  makes them controls rather than checks.
- **Part of the run is surplus to the paper**, and Section 5 says so, naming the second 24-permutation cycle
  (which the paper does not print), the exhaustive verdicts at $n=3$ and $n=2$, and the total. The 24 `w2_`
  lines and the total $13\,063\,680$ therefore support no claim of the paper.
- **Three cited sources were not read** (Section 4): Chung–Diaconis–Graham 1992 (no open-access copy),
  Hurlbert's 1990 thesis, and the 2009 survey. Neither journal version of the two 2024/2026 e-prints was
  obtainable, so preprint and published text were never diffed; Section 1 flags this, and Section 4 judges the
  risk low because the 2026 e-print independently restates the conjecture as open.

One bookkeeping discrepancy, not mathematical: the header of `verify.py` describes its two input cycles `W1`,
`W2` as "printed in Section 3 of the paper", whereas the paper prints one cycle, in Section 2, and does not
print the second at all.

## 4. How to check it

```sh
python3 verify.py
```

Python 3.9 or later, standard library only, no arguments and no external data file. The program prints one line
per check and the verdict at the end, and exits non-zero if any check fails. The recorded run used Python
3.9.25 and ends with status 0.

The first lines of `verify.output.txt` are a header rather than program output, and they carry the SHA-256 of
the program, so transcript and program can be paired. Recompute it with

```sh
shasum -a 256 verify.py
```

    0c2d660cd7fd54e095d834631d690b2943b15e064b1926d1babfd3b4baa9f117

which is the digest that header records.
