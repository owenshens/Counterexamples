# A Four-Vertex Counterexample to Conjecture 7.6 of Devriendt

`diamond-counterexample-devriendt-resistance-curvature`

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

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **35 checks, all passing**:

    VERDICT: ALL 35 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    39067b5ccbe0c28cddb1c2e5871826b2c61773a83aa995d56dab1dcba354875c

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the vertex-minimality Proposition quantifies over ALL positive conductances on the order-<=3 graphs, a continuum. The graph enumeration here is exhaustive (K2, P3, K3 are the only connected simple graphs on 2 or 3 vertices) but the conductances are a finite exact rational grid, so this census CORROBORATES the Proposition and does not establish it; the two structural identities checked above, which hold at every grid point, are checked here rather than proved. The paper does not rest on this census: its proof of the Proposition is a complete argument (c_e omega_e is the spanning-tree inclusion probability of e, so it is 1 on a bridge and less than 1 on a nonbridge; K2 and P3 are trees; every vertex of K3 meets two nonbridges), and the grid is corroboration of that argument, not a substitute for it. Also not re-run: that the definitions of normalization, resistance capacity and resistance curvature used here, and the numbering of the statement refuted and of the companion statement left standing, agree with the cited source; those are transcribed. And a caveat on the count: 8 of the 35 checks are SOLVER SELF-TESTS rather than tests of the paper -- resistance_matrix_is_a_metric, resistance_matrix_via_kirchhoff_cofactors, resistance_matrix_via_spanning_forest_ratios, edge_probabilities_sum_to_n_minus_1, edge_probabilities_from_spanning_tree_enumeration, capacity_defined_on_all_16_subsets, capacity_on_pairs_matches_closed_form, elementary_slack_census_is_24 -- each asserting a theorem true of every positively weighted connected graph, or comparing two internal recomputations, or counting something fixed by n=4. They would catch a corrupted solver but cannot fail however the exhibited object is perturbed, so the refutation is carried by the remaining 27 checks. Everything else in the paper is reproduced exactly, in exact rational arithmetic.
