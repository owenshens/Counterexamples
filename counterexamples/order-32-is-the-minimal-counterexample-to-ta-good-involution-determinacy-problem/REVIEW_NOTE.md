# Ta's Determinacy Problem for Good Involutions of Conjugation Quandles Fails at Order 32, and, Granting the Published Group Counts, at No Smaller Order

`order-32-is-the-minimal-counterexample-to-ta-good-involution-determinacy-problem`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |
| `REVIEW_NOTE.md` | this file |

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
About five seconds, exact integer arithmetic throughout, no randomness. The program prints
one line per check and a closing verdict, and exits 0 only if every check passes. The
recorded run reports **237 checks, all passing**:

    VERDICT: ALL 237 CHECKS PASS

It reads the objects exhibited in the paper --- the two order-32 presentations, the
conjugacy-class lists of §3, the presentations of the census groups, and the group-count row
of OEIS A000001 quoted in §4 --- and re-derives every quantity the paper claims. In
particular it re-proves the three set-level equivalences of Lemma 2.1 by exhaustion over all
`x, y` in `G` and all pairs of central values, so no published theorem is assumed; and it
counts `Good(Conj G)` for the two order-32 witnesses a fourth, independent way, by building
`rho` from *every* function `Cl(G) -> Z(G)` and testing the definition of a good involution
literally. Its Part 4 does the same for a pair of groups of order 48; the paper does not
discuss that pair and nothing in the paper depends on it.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    7abfe2ae5a4c3d3884f7c3badedc74cd6d1fb3a0cc31329bfd4e45a6790034e3

## Scope

The program's own statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the isomorphism-type counts A000001(n) are taken from OEIS as external input
> and are not proved here; the program only exhibits enough pairwise non-isomorphic groups to
> meet them.
>
> NOT RE-RUN: the third bucket reported in the source row, (32,4,14) with 21600 against
> 552960, is not re-derived here and nothing above depends on it.
>
> NOT RE-RUN: the search that FOUND the witnesses -- a sweep of 2450 nonabelian groups of
> orders 6..48 -- is not reproduced; this program verifies exhibited objects only.
>
> NOT RE-RUN: orders 33..47 and every order above 48 are not examined, so nothing here claims
> 32 is the ONLY order at which the triple fails.
>
> NOT RE-RUN: the definition-only scan over ALL self-inverse permutations of G is done at
> orders 6, 8, 10 only; above that the count relies on Lemma 2.1, which is itself verified by
> exhaustion on every group used.
>
> NOT RE-RUN: the table of Appendix B of arXiv:2505.08090v5 is used as an external
> cross-check for orders 6..22; it is never an input to any claim above.

Two consequences are worth stating plainly. First, the completeness of each order's census
in Proposition 1.2 --- and therefore the *minimality* half of the result --- rests on the
published isomorphism-type counts of OEIS A000001; everything else about those groups,
including their pairwise non-isomorphism, is established inside the program. Second, the
refutation itself (Theorem 1.1) rests on nothing external at all: the two order-32 groups are
rebuilt from the presentations printed in §3, their class lists are compared element by
element with the lists printed there, and the counts 128 and 1024 are each obtained four
independent ways, one of which uses no structure theory beyond the definition.
