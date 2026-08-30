# A 15-Vertex Refutation of the Kimura--Matsumoto--Sato Eternal Domination Conjecture

`a-15-vertex-refutation-of-the-kimura-matsumoto-sato-eternal-domination-conjecture`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package, no external data file, no
network and no randomness. It runs in about a minute on one core. The program prints one line per
check and a closing verdict, and exits 0 only if every check passes. The recorded run reports
**43 checks, all passing**:

    VERDICT: ALL 43 CHECKS PASS

Its inputs are the objects printed in the paper, together with published graph6 strings and
further ear placements carried in the program itself: the graph6 string and the 27-edge list of
the witness `G0` from Section 2, the graph6 string of the seed `H*`, and the fifteen published
graph6 strings of MacGillivray, Mynhardt and Virgile. Every quantity the paper states is
re-derived from those and compared: the degree sequence, triangle-freeness, maximum clique size 2,
2-connectedness, the non-planarity edge bound, `mu = 7`, `theta = 8` with an explicit clique
cover, factor-criticality, `theta(G0-v) = 7` at all fifteen vertices, `gamma^inf(G0) = 7`,
`gamma^inf(G0-v) = 7` at all fifteen vertices, and the fifteen dominating-6-set counts displayed
in Section 3. The run also covers material the paper does not claim: `alpha = 6`, `gamma = 4`, two
vertex-disjoint 5-cycles of which one is induced, `|Aut(H*)| = 48` with its eight orbit sizes on
the 78 endpoint pairs, all 78 ear attachments including tightness at every one of their deletions,
and two sparser witnesses at orders 17 and 19.

Both directions of every eternal-domination verdict carry a certificate that a **separate pass**
re-verifies from the final data only, so no answer depends on the fixed-point loop having been
iterated correctly. In the positive direction the surviving family is re-checked member by member
to consist of dominating `k`-sets closed under every attack. In the negative direction the round
at which each dominating `k`-set died is kept as a rank, and a fresh pass checks that every
dominating `k`-set admits an attack whose every legal response is non-dominating or of strictly
smaller rank; induction on that well-founded rank, not the loop, is what proves
`gamma^inf > k`. All arithmetic is on Python integers; there is no floating-point decision
anywhere.

The program is calibrated in both polarities before any conclusion is drawn, on published objects.
As **forced positives** it re-derives `gamma^inf <= 6 < 7 = theta` for all thirteen order-13
graphs of MMV's Table 10, and — the calibration that matters, because it is the exact shape of the
fifteen negative verdicts the refutation needs — for the two order-14 graphs of their Table 11 it
returns YES at `k = theta - 1 = 6` on 14 vertices, with surviving families of sizes 1741 and 1749.
So the decider demonstrably *can* answer YES in precisely the regime where Section 3 needs it to
answer NO. As **proved-silent** controls it reproduces `gamma^inf(C_m) = ceil(m/2) = theta(C_m)`
for `m = 5,7,9,11,13,15` (C_15 being a genuinely sparse member of `G0`'s own order) and
`gamma^inf = 6 = theta` for the Grötzsch graph. As **anti-controls** it confirms that two 5-cycles
sharing a vertex are factor-critical despite the cut vertex and are correctly *not* a
counterexample, and that `P_5` and `C_6` are correctly rejected as not factor-critical, so the
parity and degree gates are not vacuous.

The recorded run was produced on a cloud compute instance, CPython 3.9.25 on Linux, exit status 0.
A rerun on a different machine and operating system (macOS, CPython 3.9.6) produced
**byte-identical** output, 8350 bytes in both cases.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an exit
status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    f56835dc5e64199f330d4f4213a19468bb8367471cc1961cb097311e05ea98ff

## Scope

Four of the paper's five clauses are hand computations on the printed edge list and need no
program at all: triangle-freeness, the edge count `27 = 2n-3`, `theta(G0) = 8` with
vertex-criticality (a matching argument, since a triangle-free graph's cliques are its vertices
and edges), and the failure of the conjecture's left-hand side, `gamma^inf(G0) <= 7 < 8`, which
follows from the subadditivity of `gamma^inf` over a vertex partition together with the *published*
bound `gamma^inf(H*) <= 6`. The machine ingredient is sixteen lower-bound tests: the fifteen
bounds `gamma^inf(G0 - v) > 6`, together with `gamma^inf(G0) > 6`, which is what makes
`gamma^inf(G0) = 7` exact rather than only `gamma^inf(G0) <= 7`. The program re-derives the hand
clauses as well, but the paper does not rest on it for them.

The program's closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the MINIMALITY of the order 15. This program does not census graphs. That the least
> order of a counterexample is 15 rests on a separate exhaustion of all connected triangle-free
> factor-critical graphs of order at most 14 with n <= |E| <= 2n-3, which is not repeated here and
> is not a claim of the paper's theorem.

> NOT RE-RUN: whether some order-15 counterexample has FEWER than 27 edges. The order-15 cell was
> never censused, only the structured 78-graph family above, so the sharp edge threshold at n = 15
> is OPEN and nothing here bounds it.

> NOT RE-RUN: the asymptotic reading. Clause 5 (gamma^inf < theta) holds for every k in the ear
> tower by Lemma S, but the tightness clause is verified here only for k = 1, 2, 3, so 'no bound
> |E| <= cn+d with c > 3/2 rescues the conjecture' is a CONJECTURE, not proved; the general case
> needs a lower-bound technique for gamma^inf that does not exist in print.

> NOT RE-RUN: the two censuses that produced the sparser witnesses -- an exhaustive sweep of 3081
> two-ear configurations at order 17, and 200 of the 25740 disjoint ear triples at order 19. Only
> the two named witnesses are re-checked above; the remaining 25540 triples are UNTESTED, not
> failures, and no completeness follows from this run.

> NOT RE-RUN: planarity beyond the edge-count argument. G0's non-planarity is certified here only
> by |E| = 27 > 2n-4, which is decisive for triangle-free graphs; no Kuratowski subdivision is
> exhibited and no embedding is searched for.

> NOT RE-RUN: bibliography. No page number, conjecture number, DOI or publication date is checked
> by this program; it checks mathematics only.

Three further points a referee is entitled to.

**Only the "if" half is refuted.** Conjecture 2.1 is an equivalence; the paper makes its
right-hand side true and its left-hand side false. Nothing here bears on the forward implication
`gamma^inf(G) = theta(G) => gamma^inf(G-v) = theta(G-v)`.

**The value `gamma^inf(H*) = 6` is ours, not published.** What is in print is only
`gamma^inf(H*) < theta(H*) = 7`, i.e. `gamma^inf(H*) <= 6`, from the caption of KMS Figure 2.1;
that inequality is all the proof uses. The exact value 6 is re-derived by the program but is not
needed.

**Conjecture and statement numbering.** The label "Conjecture 2.1" and the locator *journal page
47* come from the published article, which is open access; the statement was read there and is
reproduced in full in the paper, so what is refuted does not depend on the label. The program
checks no bibliographic fact.
