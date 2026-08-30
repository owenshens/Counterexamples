# Supermodular Log-Submodular Valuations Obey the Jaccard Triangle Inequality on Every Lattice

`supermodular-log-submodular-valuations-obey-the-jaccard-triangle-inequality-on-every-lattice`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

The paper answers the second open problem of Section 7 of Badica--Badica, arXiv:2608.18194v1,
affirmatively and in full generality: for every lattice -- no distributivity, no relative
complementation, no least or greatest element, no finiteness -- and every strictly positive
monotone supermodular log-submodular valuation `f`, the generalized Jaccard distance
`1 - f(A^B)/f(AvB)` satisfies the triangle inequality. The decisive argument is a five-step
hand proof that collapses to `(P-p)(Q-q) >= (b-p)(b-q)`; the program below is a check on the
paper's numbers, not the source of the result.

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
Exact integer, `Fraction` and `Q(sqrt 5)` arithmetic throughout; no floating-point value ever
decides a check. The program prints one line per check and a closing verdict, and exits 0 only
if every check passes. The recorded run reports **27 checks, all passing**:

    VERDICT: ALL 27 CHECKS PASS

It reads the objects exhibited in the paper -- the cover relations of `B_2` and the two
valuations tabulated in the paper's sharpness section -- together with further finite objects
of its own (`M_3`, `N_5`, and a near-miss valuation `g`), and re-derives every number the paper
states. In particular it

* checks that those covers really define lattices, and that `M_3` is non-distributive and
  relatively complemented while `N_5` is neither;
* enumerates the small lattices: every lattice of order at most 4 is a chain or is
  distributive and relatively complemented, and exactly two of the five order-5 classes are
  non-distributive, namely `M_3` and `N_5`. The enumerator is controlled against a published
  sequence rather than against itself: its isomorphism-class counts for `n = 2..6` are
  `1, 1, 2, 5, 15`, which is OEIS A006966;
* verifies the two identities carrying Steps 4 and 5 of the proof as **exact polynomial
  identities in `Z[p,P,q,Q,b]`**, by sparse expansion rather than by sampling, including the
  nonnegativity certificate `(P-p)(Q-q)-(b-p)(b-q) = (P-b)(Q-q)+(b-p)(Q-b)`;
* confirms the four order facts used in Steps 1 and 2 at all 4029 ordered triples of all 25
  isomorphism classes of lattices of order at most 6;
* reproduces the two violation margins `7/30` and `1/4` of the forced-positive controls, and
  the step each one breaks (Lemma 2 and Lemma 1 respectively), so that neither hypothesis on
  `f` can be dropped;
* re-runs the finite sweep and reproduces its three totals exactly (1814 admissible
  valuations, largest margin 0, zero violations);
* re-checks a near-miss valuation `g` in exact `Q(sqrt 5)`: monotone, supermodular, reciprocal
  supermodular, **not** log-submodular, and satisfying the triangle inequality with equality.
  This is why the paper explicitly declines to call the result a converse of the source's
  necessity theorem, or half of a characterisation.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    36bfe417e3afc7cf11fb4dbc9777ee2e0bc65ea656c529103b03cc6f103aac9f

## Scope

The program's own closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the theorem of this paper quantifies over ALL lattices and all strictly positive
> monotone supermodular log-submodular valuations, and no program can range over that class.
> What is machine-checked above is (i) the two algebraic identities that carry Steps 4 and 5,
> as exact polynomial identities, and (ii) the order facts of Lemmas 1 and 2, on lattices of
> order <= 6 only. The finite sweep of section F is a CONTROL over orders 4,5,6 on a six-value
> integer grid; it is NOT an exhaustion of the admissible cone, and "largest margin 0" means
> "nothing beat the equality case", not "the cone was searched".

> NOT RE-RUN: the source paper's own theorems (th:mod, th:super-log-sub, th:super-nec-cond)
> are quoted, not re-proved; and Step 3 of the proof (the sign split on p+q-b) plus the
> statement that d_{f,J} is in general only a PSEUDOmetric are hand arguments with no
> computational content. Nothing above touches the source's open problems 1, 3 or 4, or any
> literature claim.

One further limit belongs to the paper rather than to the program, and is stated in its
sharpness section: the result is the **sufficiency counterpart** of the necessity theorem in
Section 6 of the source, for a *different* hypothesis class, and therefore is **neither its
converse nor half of a characterisation**; the source itself calls requiring `f` and `1/f`
supermodular "a softer mathematical constraint" than requiring `f` supermodular and
log-submodular, and the valuation `g` checked above keeps the two classes apart. The source's
other open problems stay open and are untouched here. The paper makes no priority claim; the
adjacent items it names, with what each leaves unsettled, are
Badica--Badica--Logofatu--Neremzoiu (Mathematics 13 (2025) 384), Kosub (Pattern Recognition
Letters 120 (2019) 36--38), Simovici and McDiarmid.
