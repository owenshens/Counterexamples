# A Counterexample to a Cartesian-Product Conjecture\\ for Node Resistance Curvature

`rooted-trees-refute-node-curvature-product-conjecture`

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
passes. The recorded run reports **26 checks, all passing**:

    VERDICT: ALL 26 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    b98fb7cc508e4e44fa40f66e6d043343b8c58b1a42130ae7d88459688ef577fc

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: pairs of rooted trees with more than 13 vertices in total, or with a factor on more than 8 vertices, were not enumerated; nor were factors that are not trees; nor factors of root curvature strictly negative (root degree at least 3), which also satisfy the hypothesis, so the census above is a statement about curvature-zero pairs only; the paper makes no minimality claim.
> NOT RE-RUN: the conjecture refuted above is TRANSCRIBED BY HAND from the source it is attributed to (Conjecture 4 of Dawkins et al., arXiv:2403.01037v1 (hand transcription)), namely "for graphs G1, G2 and vertices i in V(G1), j in V(G2): if p_i(G1) <= 0 and p_j(G2) <= 0, then p_(i,j)(G1 box G2) < 0". This program never reads that preprint, so none of the following is machine-checked here: that the statement is numbered 4; that it appears in version v1; that its authors are the ones cited; and, above all, that its logical form is the one encoded above -- a NON-STRICT inequality in each hypothesis (so that the exhibited p_i = p_j = 0 qualifies) and a STRICT inequality in the conclusion (so that p_x = 0 would already refute it). The check transcribed_conjecture_predicates_are_the_printed_inequalities pins those two inequalities in the CODE to the two inequalities in the TEXT printed at the head, and nothing more: text to source is unchecked. Every check above would still pass if that transcription were wrong, and the paper's arithmetic would still be right while its attribution was not; a referee must compare the statement printed at the head of this transcript against the source by eye. Everything else established above is derived from the two edge lists alone.
> NOT RE-RUN: the curvature p_v itself is likewise a transcribed definition, p_v(G) = 1 - (1/2) sum_{u ~ v} omega_G(u,v) with unit conductances and omega the effective resistance; it is computed here by three independent routes, but that it is the quantity the conjecture speaks of is not checked against the source either.
