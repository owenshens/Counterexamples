# A Ten-Element Graded Poset Refutes Contractibility Monotonicity for Filtered Order Complexes

`a-ten-element-graded-poset-refutes-contractibility-monotonicity-for-filtered-order-complexes`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex` | the paper, LaTeX source (amsart; compiles with `tectonic`, no external `.bib`) |
| `paper.pdf` | the paper, compiled from `paper.tex` |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run of `verify.py`, with the program's SHA-256 in the header |
| `REVIEW_NOTE.md` | this file |

Those five files are the whole folder; nothing here reads any other input.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **92 checks, all passing**:

    VERDICT: ALL 92 CHECKS PASS

It reads the object exhibited in the paper -- the fifteen cover relations displayed in
Section 2, transcribed into the program as the text block it parses -- and re-derives from
those fifteen pairs alone every quantity the paper asserts: the transitive closure and its
rank-gap profile, the recomputed cover relation, purity of rank 3, the f-vectors and Euler
characteristics of all three filtration stages, their integral homology by Smith normal form
(so torsion is visible rather than assumed absent), every intermediate object of the two
Mayer-Vietoris arguments
(the star decomposition of the order complex, the complex `Y`, the pieces `A`, `B`, `A cap B`,
the tree along which `B`'s two cones meet, and the assertion that the comparison map `phi` is
an isomorphism, checked as an equality of integer lattices), two chain identities in `A`, two
coning computations, and five controls in both polarities. The paper itself states only part
of this; the program checks more than the paper claims, never less.

The 27-step collapse of the paper's Section 4 is handled in the opposite direction, and this
is the one place where the program checks the paper rather than merely agreeing with it. The
27 rows *as printed in the paper* are transcribed into the program as a second text block
(`COLLAPSE_TABLE`) and parsed, and that table is what gets replayed: at row `i` the printed
`tau` must have exactly one proper coface in the complex remaining at that step, and it must
be the printed `sigma`; after row 27 exactly one face must remain, and it must be the vertex
`d1`. An independent free-face search of `Delta^(2)(P)`, run without reference to the table,
is then required to return the same 27 rows in the same order, so the paper and the program
cannot diverge silently. This was checked adversarially in a scratch copy: exchanging rows 23
and 24, changing the `sigma` of row 20, and deleting row 27 each make the program report
`FAIL` and exit 1.

Three of the controls are worth naming, because they are what makes the positive answers
above measurements rather than defaults. The same collapsing routine is run on
`Delta^(3)(P)`, where it stalls with 26 faces left, so it is not answering yes
indiscriminately. The homology engine is run on the six-vertex minimal triangulation of the
real projective plane and reports the torsion `Z/2`, so the "no torsion" results above are
measured. And Kitajima's own published wedge formula for his example with four levels of two
is reproduced independently by the same engine.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    aff3d7777f6082a182a453ea227ed42bc25e63894bf691beb26c8c2ec580a9b3

The run was on the control plane (darwin 25.6.0), Python 3.9.6, a single process, well under
a second of wall time, exit status 0. There is no randomness anywhere in the program: the
free-face search is a fixed sort order, and all arithmetic is exact integer arithmetic.

## Scope

**What the program does, in one sentence:** it *confirms the exhibited object* -- it takes the
ten-element poset printed in Section 2 of the paper and re-derives, from those fifteen cover
relations alone, every property of it that the paper asserts, including replaying the paper's
own printed collapse table. It does **not** re-derive *how* that object was found, and it does
not decide the underlying open problem in any direction other than the one the exhibited
counterexample settles. Concretely: the refutation of the implication of Remark 2.13 at
`k = 3` is fully re-derived; the census that located the object, minimality, priority, and the
case `k < n` are not touched at all, and the paper claims none of them.

The program's own closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the exhaustive census behind the minimality remark of Section 5 (all rank-3
> graded posets with every level of size at most three: 22,930,191 labelled posets, 49,843
> isomorphism classes, 201 witness classes). That search ran elsewhere, its stdout was not
> preserved, and nothing here reproduces it. The paper therefore asserts no minimality, and
> the refutation does not depend on it.
>
> NOT RE-RUN: the correction of Section 6 to Kitajima's Remark 4.8, p. 15 (three printed
> groups Z^12 that must read Z^11). That is a long-exact-sequence
> argument on the source's own printed integers, is logically independent of the refutation,
> and is not a computation this program performs.
>
> NOT RE-RUN: any search of the literature. Whether the k = 3 case was settled elsewhere is
> not a question a program can answer, and no claim of priority is checked here.
>
> NOT RE-RUN: the open case k < n. Nothing here bears on whether the implication of Remark
> 2.13 holds at a proper filtration stage; the witness has k = n = 3.

Two things deserve to be read as limits on the paper and not only on the program.

1. The witness has `k = n = 3`, so the complex that fails to be contractible is the full order
   complex; the paper says so in Section 5, and in the abstract, and states plainly that the
   case `k < n` is untouched.
2. The minimality of ten elements is *not* claimed anywhere in the paper, precisely because
   the search that would support it is not reproducible from this folder. Neither is the
   minimality of `k = 3` or of `n = 3`, and neither is any claim of priority.

**Position in the literature.** The paper makes no novelty claim and should not be read as
making one. It records only that the filtration is the poset case of the filtered nerve of
Di--Ivanov--Mukoseev--Zhang (reference [1]) and that the sentence it refutes is Kitajima's
Remark 2.13, p. 8, quoted verbatim.

Everything the paper actually asserts -- Lemma 1, Proposition 3, Proposition 4 and
Theorem 5 -- has a proof in the paper that a reader can follow with no machine at all: two
applications of one suspension lemma for the non-contractibility, and a printed table of 27
elementary collapses for the contractibility. The program is a second, independent pass over
the same claims, not the only route to them.
