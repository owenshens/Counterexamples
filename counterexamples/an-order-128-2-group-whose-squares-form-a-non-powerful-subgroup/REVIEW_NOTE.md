# An Order-128 2-Group whose Squares Form a Non-Powerful Subgroup

`an-order-128-2-group-whose-squares-form-a-non-powerful-subgroup`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |
| `REVIEW_NOTE.md` | this file |

Those five files are the whole folder; nothing here depends on anything outside it.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package, no group-theory library
and no external data file. The program prints one line per check and a closing verdict, and
exits 0 only if every check passes. The recorded run reports **55 checks, all passing**:

    VERDICT: ALL 55 CHECKS PASS

Its whole input is the object exhibited in the paper — the seven permutations of the eight
points printed in Section 2, the 28 relators of the presentation, and the sixteen square
roots of Table 1 — transcribed verbatim into one block at the top of the file. Every
quantity it compares against the paper's statements it derives from that input, in exact
integer and permutation arithmetic with no floating point anywhere. In particular it
rebuilds the group by closure rather than looking it up, and it verifies the presentation by
collecting all 16384 products of normal words using the printed relators and nothing else,
which is the step that makes the 28 relators a presentation of this group and not of a
larger one. It runs in about a second on a laptop.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by
an exit status, both written by the run harness. The header records the SHA-256 of the
program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    dbab58e2472790f91c5985900e326f0cd40fa922cfe7fac3408e41d4ad56cb22

## Scope

**What the program establishes.** It re-derives the paper's claim in full, and not merely the
consistency of the exhibited object. The paper's claim is an existence statement — *there is a
finite 2-group of exponent 8 whose set of squares is a subgroup that is not powerful* — and a
witness plus a re-derivation of its properties is a complete proof of it. Every property the
theorem asserts is recomputed from the printed permutations by closure (order 128, exponent 8,
the order distribution, closure of the square set under multiplication and inverses, its
non-abelianness, |S'| = 2, S^4 = 1, and hence the failure of powerfulness), together with the
hypotheses of Kourovka 21.137 themselves and the two answers of `NO`. What the program does
**not** do is search: it does not look for the witness, and it does not check that this witness
is minimal or unique — but the paper makes no minimality or uniqueness claim, so no claim in
the paper is left unverified by that omission. Corrupting any part of the transcribed object
(a permutation, a relator right-hand side, a Table 1 root, or one of the printed invariants)
makes the program exit non-zero, so it is not passing vacuously.

The program's own closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the CENSUS claims of the accompanying record are outside this program. It does
> not verify that 128 is the least order at which the phenomenon occurs, that exactly ten
> groups of order 128 are witnesses, or the counts 34 at order 256 and 2094 at order 512;
> those need the SmallGroups library. Nothing here asserts that G is the unique witness of
> its order.

> NOT RE-RUN: the library label SmallGroup(128,928). This program uses no group-theory
> library, so it can neither confirm nor deny that catalogue number; the object is named here
> only by its permutations, by Syl_2(S_8) = C2 wr C2 wr C2, and by its 28 relators.

> NOT RE-RUN: the general lemma that S' <= S^2 holds in EVERY finite group. That is proved in
> one line in the paper; the program checks it only on the 410 distinct 2-generated subgroups
> of this G.

> NOT RE-RUN: the second sentence of Kourovka 21.137, the odd-prime branch. It is a different
> question, no object here bears on it, and it remains open.

> NOT RE-RUN: every bibliographic and priority claim. Whether this observation is new is not
> a computation and is not tested here.

The same limits are stated in the paper's Section 3, and two of them bound the claim itself
rather than only the program: **only the headline sentence and the third sentence of Kourovka
21.137 are answered here — the second sentence, the odd-prime branch, is a different question
and remains open**, and **no minimality, uniqueness or priority claim is made**.

The prior-art search behind the paper's last point returned **NEAR MISS rather than a clean
negative**, and the paper's Section 3 names the near miss: Alharbi and Alghamdi, *The wreath
product of powerful p-groups*, Symmetry 15 (2023) 1987, whose Remark 3 establishes 21.137's
hypothesis — closure of the set of p-th powers — for wreath products of two cyclic p-groups,
but which never asks whether the resulting power subgroup is powerful, does not treat iterated
wreath products, and excludes p = 2. Two standard monographs on groups of prime power order
(Berkovich and Janko, Volumes 1–3; Leedham-Green and McKay) were not reachable in full text,
and since the Sylow 2-subgroup of the symmetric group on 8 letters is the most obvious 2-group
of exponent 8 to test, an exercise-level prior appearance is possible.
