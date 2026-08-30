# A Non-Standard (5,4) Broadcast of the Infinite Grid of Density 1/20

`a-non-standard-5-4-broadcast-of-density-1-20-refutes-the-blessing-et-al-grid-conjecture`

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

Python 3.9 or later, standard library only: no third-party package and no external data file, and
no floating-point comparison anywhere — every decision is between integers or `Fraction` objects.
It runs in under a second. The program prints one line per check and a closing verdict, and exits 0
only if every check passes. The recorded run reports **43 checks, all passing**:

    VERDICT: ALL 43 CHECKS PASS

It reads the object exhibited in the paper — the integer basis `[[10,0],[4,2]]` and the membership
predicate `y % 2 == 0 and (x - 2*y) % 10 == 0` — and re-derives from those two lines every quantity
the paper displays, comparing each against the printed value. Five controls of both polarities are
included, so that the checker is shown able to return each answer: the standard broadcast `T(18,4)`
is confirmed valid at (5,4), the published maximum `d = 18` and its lowest-`e` tie-break are
reproduced, no `T(d,e)` with `19 <= d <= 22` validates, the reception vector of `T(18,5)` at (4,2) is
matched entry by entry against a third party's published vector, and the known non-standard optimum
at (2,3) is accepted.

The checks are also sensitive to the object, which a reader can confirm in a scratch copy: change the
basis from `[[10,0],[4,2]]` to `[[10,0],[2,2]]` — a *different* sublattice of the same index 20, and one
the index-20 census rejects — together with the matching predicate `(x - y) % 10 == 0`, and nine checks
fail and the program exits 1, reporting minimum reception 0. Flipping a single entry of the printed
coset table fails one check and exits 1.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    9ac0d861f012352d3f807673e5fec65130d6b1a0f7fde7d8ab336034978d284c

## Scope

**What the program does and does not settle.** It does more than confirm the exhibited object. It
re-derives the object itself from the two printed lines by three independent routes and checks they
agree; it re-derives every number the paper prints about the object; and beyond the object it also runs two
exhaustions the paper does not claim — the index censuses at $19,20,21,22$ (all $\sigma(D)$
sublattices of each index, finding the witness and its mirror and nothing else) and five members of
the doubling family $2L_{t'}$, $2\le t'\le6$, each over all $4D$ classes of its own quotient. What
it does **not** do is reprove the two published lower bounds the
paper's inequality is measured against; those are quoted, and only their arithmetic at the relevant
parameters is checked. So the paper's headline inequality
$\delta_{5,4}\le1/20<1/18=\delta_{4,2}$ is machine-checked on its left half and rests on a cited
theorem on its right. The program's own closing statements, quoted from its output:

> NOT RE-RUN: the lower bound delta_{t,2} >= 1/(2(t-1)^2) of Drews-Harris-Randolph is QUOTED, not
> reproved. Only its value at t = 4 is arithmetic-checked. The strict inequality 1/20 < delta_{4,2}
> rests entirely on that published theorem.
>
> NOT RE-RUN: Shlomi's lower bound delta_{t,r} >= r / C_{t,r} is QUOTED, not reproved; only the
> arithmetic C_{5,4} = 84 and 4/84 = 1/21 is checked here.
>
> NOT RE-RUN: no search above index 22. Nothing in this program excludes a (5,4) broadcast of
> density below 1/20, so delta_{5,4} is bounded, not determined: 1/21 <= d <= 1/20.
>
> NOT RE-RUN: the censuses at index 19..22 enumerate SUBLATTICES ONLY (one tower per period).
> Periodic sets with two or more towers per period are not enumerated at any index.
>
> NOT RE-RUN: the doubling family is enumerated only for 2 <= t' <= 6. The general lemma
> "S_{t'} = 2 L_{t'} is a (2t'+1,4) broadcast for every t' >= 2" is NOT proved, so a
> Herrman-van Hintum threshold t_0 >= 13 is untouched.
>
> NOT RE-RUN: no literature search. In particular Harris-Insko-Johnson, "Projects in (t,r) Broadcast
> Domination" (Springer 2020, doi 10.1007/978-3-030-37853-0_8), is paywalled and unread.

One further limit belongs to the paper rather than to the program, and the paper's scope section
states it: the conjecture was already refuted in print at other parameters, so what the paper adds
is the instance `(t,r) = (4,2)`.

The recorded output is left exactly as run, so that it still matches the SHA-256 above. Several of
its checks concern statements the paper does not claim — the index censuses, the further members of
the doubling family, and the remarks on published tables and on the Drews–Harris–Randolph Further
Question. Those are outside the scope of the paper and nothing in the paper depends on them.

The decisive claim of the paper — Theorem 1 — needs no program at all: it is twenty sums of at most
three terms each, with the checksum lemma as a joint check against the single number 85.
