# A Fractional Weak Tiling Measure in Z_24 That Is No Convex Combination of Tiling Complements

`a-fractional-weak-tiling-measure-in-z-24-that-is-no-convex-combination-of-tiling-complements`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

The third item of the `Open problems.` list of Kadir and Fan, *Tiles and weak tiles in Z_{pq}*
(arXiv:2607.02149), is **false** for cyclic groups. That item is quoted verbatim in Section 1 of the
paper, together with its locator in the source file — `Weak_Tiling_7_05.tex`, v2, 39,297 bytes; the
list header at line 710 / byte offset 34177, the item itself at lines 716–735, its `\item` token at
byte offset 34823 and its text beginning at byte offset 34829. (The source's Definition, which
supplies the normalisation, is the `definition` environment at lines 86–88, the statement being
line 87.)
It carries no printed number, no rendered designator and no LaTeX label, which is why it is quoted
in full rather than cited by number.

The refuting object is `M = 24`, `E = {0,4,6,10}` and the measure `g` that is 1 at 0, 8, 16 and 1/2
at 7, 9, 11, 19, 21, 23. Two independent certificates are given: a separating linear functional
(`L(g) = 3` against `L ≤ 2` on all sixteen tiling complements) and an extreme-point certificate
(the nine support columns of the circulant are independent). Section 3 extends the same
construction to `E = d·{0,2,3,5}` in `Z_{12d}`, refuting the item for **every** modulus divisible by
12 and larger than 12.

Two points a referee should note about scope, both stated in the paper rather than buried:

* The item's inline "that is" clause drops the normalisation `f(0) = 1` that the source's own
 Definition carries. The paper reads that normalisation back in and says why (without it every set
 weakly tiles, and the source's main theorem would be vacuous). Under the literal unnormalised
 reading the item is false already at `M = 7`, for an unrelated reason — an *empty* hull rather than
 a *non-integral* polytope — and Section 4 states that case explicitly instead of hiding it.
* The result concerns **cyclic** groups. A published near miss, Kiss–Londner–Matolcsi–Somlai,
 *A lonely weak tile*, refutes the general finite-abelian analogue by the opposite mechanism (an
 empty hull) in a non-cyclic group. The paper cites it and says how the two differ.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only (`fractions`, `itertools`, `sys`): no third-party
package and no external data file. All arithmetic is exact integer or `Fraction` arithmetic; no
floating point value is formed anywhere, so no decision depends on a rounding mode. Runtime is
about two tenths of a second. The program prints one line per check and a closing verdict, and
exits 0 only if every check passes. The recorded run reports **44 checks, all passing**:

 VERDICT: ALL 44 CHECKS PASS

It reads the objects exhibited in the paper — `M`, `E`, and the nine printed values of `g` — and
derives from them every quantity the paper asserts. Nothing is read from a data file and nothing is
taken on trust from the paper's prose. What it covers, part by part:

* **A** — `g ≥ 0`, `g(0) = 1`, `g` is not 0/1, the twenty-four convolution values `(1_E * g)(y) = 1`
 (reported as `violations: []`), the mass `6 = M/|E|`, `E − E ⊆ 2Z_24`, the pin lemma of Section 1
 applied to `supp(g)`, the even-translate partition and the odd-translate double cover that the
 paper's proof uses, and that `g` is not positive definite (it is real but not symmetric).
* **B** — all sixteen tiling complements of `E`, found by **exhaustive enumeration of all
 C(24,6) = 134,596 six-element subsets** of `Z_24`, checked against the structural description
 `(2a + <8>) ∪ (2b+1 + <8>)`, and the four through 0 checked against the list printed in the paper.
 It also checks that those four are points of the pinned polytope, i.e. that the hull the paper
 separates `g` from is **nonempty** — which is what makes the refutation non-trivial.
* **C** — `L(g) = 3`, the multiset of `L(1_T)` over all sixteen complements, and the separation.
* **D** — both rank statements (rank 9 on 9 columns over all 24 rows; rank 8 on the 8 free columns
 of the reduced pinned system over the 20 rows of `Z_24 \ E`) and the determinant 4 of the odd
 block. These are two figures about two different matrices and the program checks both and prints
 which is which.
* **E** — four **anti-controls**: the pinned polytopes of `(3,{0,1})`, `(7,{0,1,3})`, `(4,{0,1,2})`
 and `(6,{0,1,3})` are *empty*, so those instances are not counterexamples. This is the check that
 the machinery reads the statement the source actually made rather than a weaker one.
* **F** — the unnormalised `M = 7` instance: `|det| = 24` for the Fano circulant, the unique
 solution `1/3` everywhere, and an exhaustive search over all 127 nonempty subsets of `Z_7`
 confirming that `{0,1,3}` has no tiling complement at all.
* **G** — the `Z_12` seed: the pinned polytope of `{0,2,3,5}` in `Z_12` is the single **integral**
 point `1_{{0,4,8}}`, which is why the phenomenon is invisible at `M = 12`.
* **H** — the infinite family, twice over: by the two block ranks that constitute its whole proof
 (valid for all `d`), and end to end for `d = 2,…,8`, i.e. `M = 24, 36, 48, 60, 72, 84, 96`. It
 also checks that the `d = 2` member is literally the object of Part A.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an exit
status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 f7da82095f83bf7ffb10d2872450a2c9ba9de0ec01bfd0c570c6c923e7c3515a

That run is the one recorded here: `verify.py` was executed in this folder under CPython 3.9.6 on
macOS (arm64), exit status 0, and the transcript is its complete captured stdout — nothing was
elided or truncated.

**Where the rest of the record lives, and where it has holes.** The decisive mathematics is a hand
proof and needs no machine; the paper's Sections 2 and 3 are checkable with pencil and paper. Three
further programs were run during the investigation and are filed, byte-for-byte and with SHA-256
digests, in this project's artifacts directory `runs/wave23/artifacts/t8101/` alongside its
`MANIFEST.json`. They are *not* part of this folder, they share no code with `verify.py`, and this
note does not restate their figures as claims of the paper. Quoting only what that manifest itself
records:

* `t8101_verify.py` — an independent verifier of the `M = 24` object. Run on a fleet slot, status Success, RC 0. Its
 captured stdout **is** filed (`t8101_verify.out`, 3,373 bytes, recorded as complete).
* `t8101_census.py` — an exhaustive census of all 16,777,088 pairs `(M,E)` with `8 ≤ M ≤ 24` and
 `0 ∈ E`. Run on a fleet slot, status Success, RC 0. **Its raw stdout is NOT in the
 record**: the manifest states that the dispatch filed no artifact and that the planned S3 objects
 were reported absent. Only transcribed figures survive.
* `t8101_a2_seeds.py` — a sweep for further seed orbits. Run on a fleet slot, dispatched
 2026-08-31T07:08:41Z, which the manifest records as **after** the decision on this problem; it did
 not feed the verdict. **Its raw stdout is also NOT in the record.**

Two consequences are carried into the paper rather than papered over. First, the census's
conclusion that `M = 24` is the least modulus carrying such an object rests on a run whose
transcript is missing, so **the paper claims no minimality at all** (Section 4). Second, the seed
sweep's further orbits at moduli 18, 20 and 24 are likewise not reproduced here, so the family the
paper states and `verify.py` checks is exactly `{12d : d ≥ 2}` and no more.

## Scope

What is refuted is the quoted item as a statement about cyclic groups `Z_M`, under the source's own
Definition (a nonnegative `f` with `f(0) = 1` and `1_E * f = 1_G`). The program's own closing
statement of what it does not cover, quoted from its output:

> NOTE SCOPE: what is verified above is item 3 of the source's open problems under the source's OWN Definition -- a nonnegative f with f(0) = 1 and 1_E * f = 1 -- in a CYCLIC group. NOT RE-RUN and NOT CLAIMED here: (i) that M = 24 is the least modulus carrying such an object -- an exhaustive census over 8 <= M <= 24 was run elsewhere in this project and its raw transcript is not part of this folder, so minimality is left open here; (ii) the positive-definite variants asked in the same list, which this object cannot reach at all (check A10: g is not even symmetric, hence not positive definite); (iii) the two continuous Kolountzakis-Lev-Matolcsi versions, for convex bodies in R^d and for finite unions of intervals in R, which a finite cyclic object says nothing about; (iv) further seed orbits at moduli 18, 20 and 24 recorded in this project's artifacts but not reproduced here, so the family verified above is exactly {12d : d >= 2} and no more; (v) whether a cyclic group contains a LONELY weak tile, i.e. a set with a weak tiling measure and no tiling complement at all, which is untouched by everything above.

One further point of novelty, named rather than dismissed: that an exact-cover linear program can
fail to be integral is textbook. What is new here is the *identification* of a non-integral instance
whose `E` is nevertheless a genuine tile with sixteen tiling complements — so the hull in the quoted
item is nonempty and the failure is not the trivial one — together with the infinite family that
follows from it.
