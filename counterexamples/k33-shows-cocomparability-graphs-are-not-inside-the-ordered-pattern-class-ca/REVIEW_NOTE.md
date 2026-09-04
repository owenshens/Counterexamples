# Review note

Paper: *$K_{3,3}$ is a cocomparability graph with no $P_a$-avoiding vertex order*

Files in this folder: `paper.tex`, `paper.pdf`, `verify.py`, `verify.output.txt`, and this note.
The cited sources are not here; the bibliography of `paper.tex` gives their arXiv numbers and
journal data.

## 1. What the paper claims

Feuilloley and Habib (arXiv:2112.00629, quoted from v2) classify grounded intersection graph
classes by forbidden ordered patterns. Their Open problem 2, quoted verbatim in Section 1, asks
whether the class $C_a$ — the graphs admitting a vertex order with no occurrence of the four-vertex
ordered pattern $P_a$ — and the cocomparability graphs are incomparable, offering in a parenthetical
the alternative that the cocomparability graphs are *contained* in $C_a$. Section 1 quotes the
source's pattern definition, the names $a=(1,2)$, $b=(2,3)$, $c=(3,4)$, $d=(1,4)$ of the four
optional pairs, and the drawing convention under which a pair outside $S$ is unconstrained, and
reads $P_a$ as $v_1<v_2<v_3<v_4$ with $v_1v_3, v_2v_4\in E$ and $v_1v_2\notin E$, the pairs
$(2,3)$, $(3,4)$, $(1,4)$ unconstrained.

**Theorem 1** settles the direction the question leaves open: $K_{3,3}$ (parts $A=\{0,1,2\}$,
$B=\{3,4,5\}$, $n=6$, $m=9$) is a permutation graph, hence a cocomparability graph, and none of the
$720$ linear orders of its vertices avoids $P_a$. So the cocomparability graphs are not contained
in $C_a$, and the parenthetical alternative is false. Both halves are proved by hand in Section 2.
The first half exhibits $\pi=4\,5\,6\,1\,2\,3$, whose nine inversions are exactly the cross pairs,
and cites Pnueli–Lempel–Even for the identification of the permutation graphs with the graphs that
are simultaneously comparability and cocomparability. The second half derives the displayed
equivalence (1) — $P_a$ occurs iff two vertices of one part sit at positions $i<j$ and two vertices
of the other part both sit at positions $>j$ — and then contradicts $q_2<p_2$ with $p_2<q_2$ by
pigeonhole. A concrete occurrence is printed for the alternating order $0\,3\,1\,4\,2\,5$ at
positions $1<3<4<6$. The abstract states that the argument "uses no computer".

**Corollary 2** concludes that $C_a$ and the cocomparability graphs are incomparable, the answer
Open problem 2 asks for, and states that only the Theorem 1 direction is new here: the other
direction — $C_6\in C_a$ via the cycle order, $C_6$ not a cocomparability graph by Gallai — is the
source's own assertion, quoted in Section 1. **Corollary 3** (Section 3) concludes
$C_a\subsetneq{}$ interval filament graphs, the inclusion the source's lead-in singles out as the
lowest edge of its Figure 1.

Section 1 further contains a Remark on two conventions of the source (the two spellings of the class
name, and a "connected" clause the paper treats as a normalisation), and an argument that the
published incomparability result for the sibling class $C_b$ — Observation 6.9 of the max
point-tolerance paper, restated by Jelínek and Töpfer — does not settle Open problem 2.

## 2. What the program checks

`verify.output.txt` is a recorded run of `verify.py`. It closes with `VERDICT: ALL 43 CHECKS PASS`
and status 0. The 43 `PASS` lines fall into seven steps; the counts below are read off that
transcript, and its Step 2 and Step 3 headings name the two halves of the paper's Theorem 1
explicitly.

* **Step 1, 5 checks** — the object of Section 2: the printed edge list is exactly the cross pairs
  of the printed parts; $n=6$, $m=9$; 3-regular; triangle-free; complement $=01\,02\,12\,34\,35\,45$,
  two components of order 3.
* **Step 2, 5 checks** — the second half of Theorem 1. `0 of 720 orders avoid P_a` twice over, by an
  incremental decider and by a prune-free scan of all $C(6,4)=15$ position quadruples in each order,
  with 0 disagreements between the two on all 720 orders; the occurrence printed for the alternating
  order is confirmed to be one (non-edge 0–1, edges 0–4 and 1–5); and the pigeonhole step of the
  hand proof holds in every one of the 720 orders.
* **Step 3, 4 checks** — the first half of Theorem 1. The route the paper takes: $\pi=4\,5\,6\,1\,2\,3$
  has exactly the nine inversions $E(K_{3,3})$. Three routes the paper does not take: the order
  $0\,1\,2\,3\,4\,5$ is umbrella-free among its $C(6,3)=20$ triples; one explicit orientation of the
  complement is transitive (6 arcs $=|E(\text{complement})|$, 0 reversed pairs, 0 violations); and a
  transitive orientation of the complement is found again by exhaustive search over all $2^6=64$
  orientations.
* **Step 7, 13 checks** — controls of both polarities. Three bear on the paper: $C_6\in C_a$ via the
  cycle order and $C_6$ not a cocomparability graph, the two facts the source's parenthetical
  asserts and Corollary 2 quotes; and the check that a certain order of an 8-vertex graph *built by
  the program* (the transcript prints it as $C_8$ with the chords 04 and 26, $m=10$) avoids
  $P_{ab}$ and contains $P_a$, which the transcript says "fails if a and b are swapped" — this is
  the check Section 4 describes as the one that would fail if the optional pairs $a$ and $b$ were
  read the wrong way round.
* **Steps 4, 5, 6 ($9+2+5=16$ checks) and the remaining ten controls of Step 7** concern statements
  the paper does **not** make: a membership criterion for complete multipartite graphs, that every
  complete multipartite graph is a permutation graph, that no graph on at most five vertices lies
  outside $C_a$ (with three pairwise non-isomorphic order-6 witnesses — $K_{3,3}$, the prism, the
  octahedron — exhibited), and controls on $K_4$, $C_5$, $C_7$, $K_{2,3}$, $K_{2,m}$, $K_{3,4}$,
  $P_\emptyset$, and edgeless and complete graphs. Section 4 of the paper warns of this in advance,
  including that the program's headings refer to numbered statements that do not appear in the
  paper: the transcript's Step 4 "Theorem 2" and Step 5 "Corollary 3" are **not** the paper's
  Theorem 1, Corollary 2 or Corollary 3. Section 4 states that nothing above depends on any of those
  checks; by the counts here they are 26 of the 43.

## 3. What the program does not check

* **Theorem 1 is a hand proof; the program is a control.** Step 2's 720-order exhaustion does
  independently confirm the theorem's second half for this one graph, but it is not the argument the
  paper gives, and the abstract states the argument uses no computer.
* **Quoted, not re-derived.** That the permutation graphs are exactly the graphs that are both
  comparability and cocomparability is quoted from Pnueli–Lempel–Even. The inclusion
  cocomparability $\subseteq$ interval filament, on which Corollary 3 rests, is quoted from
  Feuilloley and Habib: the proof of Corollary 3, Section 4 and Section 5 all say so, and the
  transcript's closing SCOPE / NOT RE-RUN block says the same. That $C_6$ is not a cocomparability
  graph is attributed to Gallai; the program reproduces it for $C_5$, $C_6$ and $C_7$ as a control,
  but Corollary 2's second direction is presented as the source's assertion, not as proved here.
* **A quoted equivalence is only sampled.** Per the SCOPE / NOT RE-RUN block, the agreement of the
  umbrella characterisation of cocomparability with "the complement has a transitive orientation" is
  checked exhaustively at $n=4$ and on every named graph in the run, and is otherwise quoted from
  the literature, not proved.
* **The quantifiers of the non-paper statements are sampled.** From the same block: the program's
  own "Theorem 2" is cross-checked against exhaustive search for $2\le n\le 7$ only and its explicit
  sufficiency orders for $2\le n\le 10$; the minimum-order statement is exhaustive over all labeled
  graphs on at most 5 vertices, and at order 6 the run exhibits three pairwise non-isomorphic
  witnesses and does NOT claim they are the only ones.
* **Nothing about orders 7 and above**, and nothing about a finite forbidden-subgraph
  characterisation of $C_a$ — disclaimed both by that block and by Section 5, "What is not settled",
  which also records that no characterization of $C_a\cap{}$cocomparability and no structural
  (non-pattern) description of $C_a$ is offered (the latter is Open problem 3 of the source,
  untouched), and that $K_{3,3}$ is not claimed to be the only witness of its order.
* **The quotations are transcriptions.** Section 1's verbatim quotations from Feuilloley and Habib
  (including a word missing from the source's lead-in, left unrepaired), and the statements Section 1
  takes from the max point-tolerance paper and from Jelínek and Töpfer — among them that $K_{3,3}$
  lies inside $C_b$, which the paper explicitly does not re-derive — rest on sources not in this
  folder. The program neither can nor does check that any quotation is faithful.

## 4. How to check it

```sh
python3 verify.py
shasum -a 256 verify.py
```

`verify.py` takes no arguments, reads no input file, imports only `sys` and `itertools`, and uses
integer arithmetic only, so no decision depends on a rounding mode. It prints one
`PASS <name> [detail]` line per check and exits 0 if and only if every check passes. The header of
`verify.output.txt` carries the SHA-256 of the program that produced that output, so the transcript
and the program can be paired. Computed with `shasum -a 256` from the `verify.py` shipped here it is

    1fb55abd49ea121dc986708c812a8668727567b3b9c40aa7f0ec5b51a70dfdff

which is the digest that header records. The header also records the interpreter as Python 3.9.25;
Section 4 of the paper states Python 3.9+, standard library only.
