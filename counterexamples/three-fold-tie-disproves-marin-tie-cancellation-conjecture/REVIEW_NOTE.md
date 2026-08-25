# A Counterexample to Mar\'in's Conjecture 5.14

`three-fold-tie-disproves-marin-tie-cancellation-conjecture`

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
passes. The recorded run reports **89 checks, all passing**:

    VERDICT: ALL 89 CHECKS PASS (INCOMPLETE -- see G6/G7 above)

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    ac8f7b5498e2f669953415c994a47e8aa37a8d771d5c76faff40df95ad217579

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> GAP REGISTER -- what a green verdict above does NOT establish
> G1 IMPORTED, NOT DERIVED: c_top(Lambda) = c(Lambda,T(Lambda)) lies in
> {+1,-1}  (Marin Prop. 5.10(iv)).  This program does not build
> tau^B_t, the alternant expansion or the Delta_t division, so it
> cannot confirm the membership.  MITIGATION: all 2^3 sign patterns
> are enumerated, so the refutation is independent of WHICH signs
> occur -- but if any c_top were 0, or +-2, the argument fails.
> G2 PARTLY IMPORTED: the rule that reads T(Lambda) off V(Lambda) (Prop.
> 5.10(v) plus the 5.19(ii) indexing).  DERIVED HERE (check 9b), once
> supp nu is indexed as the source's proof fixes it -- V restricted to
> the complement of a transversal, sorted decreasingly, which at t = 3
> means deleting one entry of V from the nonzero folded class: that
> X(Lambda) is the coordinatewise dominant index (the content of Prop.
> 5.10(v)), and that every index COORDINATEWISE above X(Lambda) -- so
> every term of the inversion sum but the shift (1,1) -- misses supp nu
> TOGETHER WITH ALL FOUR OF ITS W(D_2) IMAGES.  That last clause is the
> reflection route back into the support which the paper's four-
> sentence justification of eq. (4) leaves implicit: the summand of the
> 5.19(ii) inversion is the W(D_2)-ANTISYMMETRIC extension tilde_nu, so
> coordinatewise domination alone would not settle it.
> STILL IMPORTED: the 5.19(ii) inversion formula itself, the global
> sign s_r, and that supp nu contains no index outside the deletion set
> just described.  The only external test of the composite rule remains
> the calibration: Marin's own Remark 5.13 reports a two-fold tie at
> T = (2,2) for beta = (8,7,6,3,2,1,0), and the same code reproduces
> it.  One data point.
> G3 NO LONGER IMPORTED (closed by check 5b): 'Lambda is a dominant B_3
> weight, hence has at most 3 parts' is not needed -- every partition
> with 4 or more parts and |Lambda| <= 18 is swept and has
> a^B_Lambda = 0, so no such Lambda can join the tie.  What remains
> imported is only that Lambda is a PARTITION (dominant, integral),
> i.e. that no weight with a negative or fractional coordinate enters
> the branching sum.
> G4 IMPORTED: Littlewood's stable orthogonal restriction formula and
> its validity for 2*ell(lambda) <= N.  Only its parts-even vs
> columns-even CONVENTION is discriminated here (check 6).
> G5 IMPORTED: M(beta) = max{T(Lambda) : a^B_Lambda != 0}, the maximum
> restricted to the SUPPORT of a^B.  That restriction is load-
> bearing, and the next check shows it:
> Lambda = (8,8,8) is a perfectly good 3-row even partition with
> V = (21, 19, 17), X = (21, 19), T = (8, 8), which exceeds M = (6, 6) coordinatewise.
> It is excluded only because a^B_(8,8,8) = 0 (it is not inside
> lambda = (6,6,6)).  So if Marin's max ran over all Lambda
> rather than the support, M(beta) would not be (6,6).
