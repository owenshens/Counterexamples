# Review note

Files in this folder: `paper.tex`, `paper.pdf`, `verify.py`, `verify.output.txt`, and this note.
Nothing outside the folder is needed to check the paper.

## 1. What the paper claims

The paper answers the closing suggestion of the final remarks of Bridi, Martins, Marquezino and de
Figueiredo, *The only Class 0 Flower snark is the smallest* (arXiv:2505.22941v2), quoted verbatim in
§1: whether the pebbling number of the House of Graphs graphs `#1395` (the truncated tetrahedron),
`#44170` and `#44172` is indeed 13, and whether `#6698` is Class 0 — equivalently, whether the Flower
snark $J_3$ is the only Class 0 cubic 12-vertex graph of girth 3 and diameter 3. Both halves are
settled, for the four adjacency lists printed in §2 ("The graphs"):

* **Theorem 1**: $\pi(\#1395)\ge 14$, $\pi(\#44170)\ge 14$, $\pi(\#44172)\ge 14$; in particular the
  suggested value 13 does not hold for any of the three. Each bound is carried by one 13-pebble
  unsolvable configuration ($W_1,W_2,W_3$, Table 1), certified by **Lemma 4** through the move-closed
  reachable set $\mathrm{Reach}(C)$, of sizes 6358, 1912 and 1711.
* **Theorem 2**: $\pi(\#1395)=\pi(\#44170)=\pi(\#44172)=14$ and $\pi(\#6698)=12$, so `#6698` **is**
  Class 0; it is not isomorphic to $J_3$ (**Proposition 3**), so $J_3$ is not the only Class 0 graph
  in that parameter cell. The upper bounds come from the level-indexed exhaustive sweep of §4 ("The
  exact values"), organised by **Lemma 5** (downward closure) and **Lemma 6** (one-move recurrence);
  the level sizes are Table 2.

§1 and §2 are explicit that `#1395`, `#6698`, `#44170`, `#44172` are *names for the printed lists*,
that the lists were served by the House of Graphs REST API for those ids but that this identification
is not verified (§1 records that the database's search endpoint returned HTTP 401), and that nothing
in the paper is a claim about the database's holdings. Three of the four lists are pinned instead to
the source's own figure of 12-pebble unsolvable configurations by an exact edge-set identity, $18=18$
edges with empty symmetric difference under API index $=$ label $-1$; for `#6698`, which that figure
does not show, there is no such pin. Nothing in the source is contradicted: §3 ("Three certificates")
re-certifies the source's own three 12-pebble configurations unsolvable (closure sizes 1262, 1251,
483), confirming its $\pi\ge 13$, and Table 2 records $\pi(J_3)=12$ for the control $J_3$.

## 2. What the program checks

`verify.output.txt` records a run of `verify.py` in `FULL` mode under Python 3.9.25, standard library
and exact integer arithmetic only, ending `VERDICT: ALL 110 CHECKS PASS` and `program exited with
status 0`. By block:

| block in the transcript | checks | paper claim it corresponds to |
|---|---|---|
| Step 1, the exhibited graphs as objects | 10 | the ambient clauses of Proposition 3 for all five graphs: simple, cubic, 12 vertices, 18 edges, connected, girth 3, radius 3, diameter 3 |
| Step 2, invariants by two routes | 22 | triangle counts, 4-cycle counts by direct enumeration and by the trace formula, automorphism-group orders and vertex orbits — support only; the paper asserts none of these numbers |
| Step 3, identification and non-isomorphism | 5 | Proposition 3: `#1395` $\cong$ the truncation of $K_4$, the graph read off the source's $J_3$ figure $\cong$ the standard $J_3$ recipe, and the five graphs pairwise non-isomorphic (10 pairs). The additional check that the figure's $J_3$ is Petersen with one vertex expanded to a triangle is not a claim of the paper |
| Step 4, certificates | 19 | Theorem 1 and §3: each $W_i$ has total 13 and $C(r)=0$, each closure is move-closed and never puts a pebble on the target, sizes 6358 / 1912 / 1711, hence $\pi\ge 14$ for the three graphs; likewise the source's three 12-pebble configurations, sizes 1262 / 1251 / 483. Two further 12-pebble certificates $X_1$ (size 198) and $X_2$ (size 313) appear here and are not in the paper |
| Step 5, the two lemmas tested on data | 3 | Lemma 5 on 30 one-pebble removals from the certificates; the brute-force odometer against $\binom{s+10}{10}$; Lemma 6 against forward reachability on 2799 configurations of `#44172`, target 5, 0 disagreements |
| Step 6a, the exhaustive sweep | 38 | Theorem 2 and every printed row of Table 2: the per-target vectors at totals 12 and 13, the sums $(1980,24)$, $(734,52)$, $(542,12)$, the all-zero rows for `#6698` and $J_3$, $\pi=14$ for the three graphs, $\pi(\#6698)=\pi(J_3)=12$, emptiness of level 14, and constancy of the per-target counts on automorphism orbits |
| Steps 6b and 6c, sixteen graph6 strings | 4 + 3 | no claim of the paper; see §3 below |
| Step 7, the arithmetic the paper prints | 6 | $\binom{22}{10}=646646$, $\binom{23}{10}=1144066$, $12\cdot 646646=7759752$, the printed row sums, and 13 pebbles per witness. The sixth, $41186376$, is not printed in the paper |

A transcript note after Step 4 records that Theorem 1 is complete at that point and needs no
exhaustive search.

## 3. What the program does **not** check

* **The upper bounds are the program's, not a hand proof; the program here is not a control.** §5
  ("Scope") says so itself: that $U_{12}=\emptyset$ for `#6698` and $U_{14}=\emptyset$ for the other
  three is proved by exhaustion over a search space whose completeness is itself proved (Lemma 5),
  but it is exhaustion, and no weight-function or LP-dual certificate for $\pi(\#6698)\le 12$,
  checkable with integer arithmetic and no program, is offered. Theorem 1 asks for less trust, resting
  on the closure of three single configurations, but §5 notes those closures are not printed either,
  so it too is a machine check, only a small and local one.
* **The paper is not parsed.** The closing SCOPE note of the run states this: the adjacency lists,
  configurations and graph6 strings are transcribed *into* `verify.py`, so the transcription of the
  paper's printed objects and figures into the program's inputs, and the coverage of every quantity
  the paper asserts, are not checked. Table 1, Table 2 and Proposition 3 must be compared with the
  transcript by eye.
* **Quantifiers proved but only sampled.** Lemmas 5 and 6 are proved in §4 for all configurations;
  Step 5 only samples them — 30 one-pebble removals for Lemma 5, and Lemma 6 up to level 5 on one
  graph and one target. Their general validity is the hand proof, not the run's.
* **The edge-set pin of §1 is not among the 110 checks.** The $18=18$ identity tying three of the
  lists to the source's figure is asserted in the paper; it is not recomputed in the transcript.
* **The database identification is not checked.** A closing NOT RE-RUN note records that what the
  House of Graphs database holds, and the identification of the four adjacency lists with the ids
  `#1395`, `#6698`, `#44170` and `#44172`, are not re-derived — the lists are inputs. The paper's
  verification remark at the end of §5 says the same, and §1 says the identification is not verified
  there either.
* **No count of the parameter cell.** The second NOT RE-RUN note records that no completeness
  statement is verified: from the sixteen transcribed graph6 strings the run re-derives only that they
  are sixteen pairwise non-isomorphic cubic 12-vertex girth-3 diameter-3 graphs, i.e. at least
  sixteen. §5 likewise says how many such graphs there are is not settled here, and that the source's
  report of five is neither used nor contradicted.
* **Steps 6b and 6c lie outside the paper.** They concern those sixteen graph6 strings, their
  pairwise non-isomorphism (120 pairs), the fact that the five named graphs are among them, and their
  $\pi$-distribution $\{12:2,\ 13:5,\ 14:9\}$, the eleven unnamed graphs being swept one target per
  vertex orbit. The paper's verification remark in §5 says the run reports the pebbling numbers of
  eleven further such graphs "which nothing here claims".
* Also outside the paper: the Step 2 invariants, the certificates $X_1$ and $X_2$, the Tietze
  identification in Step 3, and $41186376$ in Step 7.

## 4. How to check it

```
shasum -a 256 verify.py
python3 verify.py            # the recorded run: FULL mode
python3 verify.py --quick    # skips Steps 6a and 6c; the upper bounds are then NOT verified
```

The digest of the shipped `verify.py` is

```
423f7598635fd85d9472fd65c347f2ec0319195ee1d2a8374f0d9925e86965ae
```

and the header of `verify.output.txt` carries the same `sha256:` line beside `program: verify.py`, so
the transcript and the program in this folder pair up. `-jN` sets the number of parallel workers;
`--quick` prints an extra SCOPE note saying $\pi\le 14$ and $\pi(\#6698)\le 12$ are not verified by
that transcript. Only the closing `VERDICT` line and the exit status decide a run.
