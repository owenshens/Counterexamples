# The Face-Angle Map of a Tetrahedron Is Not Injective on the Danger Cylinder

`the-face-angle-map-of-a-tetrahedron-is-not-injective-on-the-danger-cylinder`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

The paper answers Problem 3 of Nikitenko and Nikonorov, *On face angles of tetrahedra with
a given base* (Results Math. **81** (2026), no. 3, Paper No. 60; arXiv:2505.22374v4), in the
negative, by exhibiting two distinct points of `Cyl ∩ Π⁺` with the same image under the
face-angle map `F`.

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

Python 3.9 or later, standard library only (`fractions`, `decimal`, `itertools`): no
third-party package and no external data file. The program prints one line per check and a
closing verdict, and exits 0 only if every check passes. The recorded run reports **84
checks, all passing**:

    VERDICT: ALL 84 CHECKS PASS

It reads the object exhibited in the paper -- the four rational points `A = (1,0,0)`,
`B = (15/17, 8/17, 0)`, `C = (0,-1,0)`, `P = (4/5, 3/5, 0)` -- and nothing else, and
re-derives every quantity and every inequality the paper asserts.

Every decision is made in exact arithmetic. The base data are rationals. The two apex
heights are quadratic irrationals, so the decisive comparisons are carried out in the real
quadratic field `Q(sqrt 35)`, implemented in the program with an exact sign test; no
floating-point number is ever compared. The `decimal` module appears only to confirm that
the truncated decimal expansions *printed* in the paper are correct digit strings for
numbers already identified exactly; those checks are labelled `digits-*` and nothing in the
refutation rests on them.

Beyond the witness itself the program runs further checks that bear on statements the paper
does not make: that the pair `(sigma, pi)` solving all three coordinate conditions at once is
the same for every concyclic configuration, on 4004 concyclic rational configurations; the
identity `sigma = 2(P·H - R²)` and `sigma < 0` on every acute base in that scan; and, as a
control on the coordinate order and the arc convention, the value
`F(D) = (-cos BAC, cos ABC, cos ACB)` that the source paper records for `D` on the
circumcircle. None of these is needed for the refutation.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    e454c9a28c4ea042410249dd8ed74a155e469c79548c5893e29dc94e76b0aefd

## Scope

The refutation in the paper needs no computer at all: the object is six rational numbers and
the decisive step is three rational identities plus three sign checks, all displayed in the
paper. The program is a check on those, not a substitute for them.

**Which of the two it is.** The program **confirms the exhibited object**. Theorem 1 -- the
whole refutation -- is re-derived end to end from the four points and nothing else, so on the
paper's only claim the program is a full independent check. Lemma 2 is proved by hand in the
paper and is not re-derived by the program.

The program's output was recorded when the paper carried further material, since removed;
its scope notes therefore refer to statements the paper no longer makes, and to section and
statement numbers the paper no longer uses. Quoted from that output:

> NOT RE-RUN: pairs of points of Cyl cap Pi^+ lying over DIFFERENT base points. The paper
> settles only the same-vertical-generator cell, and nothing here touches the rest of
> Problem 3.
>
> NOT RE-RUN: the dimension, topology, triple points and global structure of the
> self-intersection locus of the surface FC.
>
> NOT RE-RUN: the behaviour on Cyl cap Pi^0 (the plane r = 0), which the source paper
> already settles and which is a different subdomain.
>
> NOT RE-RUN: Proposition 3 and Lemma 2 as THEOREMS. This program CONFIRMS THE EXHIBITED
> OBJECT: it re-derives every number and inequality of Theorem 1 from the four points A, B,
> C, P, and it checks the universal (sigma, pi) of Remark 3 on a finite rational scan. The
> general statements -- the reduction modulo the quadratic (Lemma 2), the uniqueness of
> (sigma, pi) from linear independence, and hence the "only if" half of Proposition 3 -- are
> proved in the paper, by hand, and are not re-derived here.
>
> NOT RE-RUN: the necessity proof of Section 4 as a THEOREM. Part 4 above verifies the
> identity sigma = 2(P.H - R^2) and confirms sigma < 0 on every acute base in a finite
> rational scan; the step from that scan to "every acute base" is the Cauchy-Schwarz and
> Euler argument in the paper, which is not re-derived here.
>
> NOT RE-RUN: the claim that strictly obtuse is NOT sufficient. The three auxiliary bases
> the paper names in Section 4(iii) -- the empty cells of (5,3,3) and of the isosceles
> (120,30,30) family, and the nonempty cell of the scalene (7,5,3) -- are not re-derived
> here. Nothing in Theorem 1 depends on them.
>
> NOT RE-RUN: every statement in Section 5 about the P3P literature. Those are quotations
> from other papers, not computations, and this program does not fetch or check any external
> source.
>
> NOT RE-RUN: the conversions of the exact cosines into degrees quoted in the paper.
>
> NOT RE-RUN: any minimality or uniqueness claim. The paper makes none, and this program
> exhibits one witness and one universal identity, not a classification.

Two further limits belong to the paper rather than to the program, and its scope section
states them there. First, both apexes of the witness lie on one vertical generator of
`Cyl`; pairs over different base points are not considered, and nothing here describes the
non-injectivity locus. Second, the cylinder over the circumcircle is the P3P *danger
cylinder* under another name, and that literature has not been surveyed exhaustively: items
not obtained include Rieck's papers on repeated solutions (2012, 2014, 2018), Rieck's
J. Math. Imaging Vis. **66** (2024), no. 1, 75--91 (Zbl 1537.51020), the 2022
Wang--Zhang--Hu article on the deltoidal surface (Zbl 1554.68181), Zhang--Hu 2006, and the
zbMATH review bodies, which were licence-blocked to us. No priority is claimed; the theorem
is checkable on its own terms.
