# A Complete Set of Mutually Unbiased Bases in Dimension 27 Whose Projective Toric Design Is Neither a Group nor a Coset of One

`a-complete-mub-set-in-dimension-27-whose-projective-toric-design-is-not-a-group`

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
It runs in a few seconds. The program prints one line per check and a closing verdict, and
exits 0 only if every check passes. The recorded run reports **35 checks, all passing**:

    VERDICT: ALL 35 CHECKS PASS

Its whole input is what the paper prints: the field presentation
`GF(27) = F_3[t]/(t^3 - t - 1)` with the coordinate index of equation (3), the 27 digit words of
Table 1, and the three generators of `Lambda` displayed after Table 1. Those words are
transcribed verbatim into the program's `PRINTED_FORMS` and `PRINTED_LAMBDA_GENS`. Everything
else -- the field arithmetic, the trace, the spread maps of (6), the design `X`, the character
sums and the closure sweeps -- is rebuilt from that data. The counts the program reports include
`265,356`, `9,477`, `255,879`, `142,884`, `142,506`, `351`, `598` of `729`, `435,942` of
`531,441`, `0` of `729` and `19,683`; the paper states those it needs.

All arithmetic is exact in `Z`, `F_3` or `Z[zeta_8]`; no floating-point value decides any check.
Sums of cube roots of unity are settled by the integer identities
`n_0 + n_1 w + n_2 w^2 = 0 <=> n_0 = n_1 = n_2` and
`|n_0 + n_1 w + n_2 w^2|^2 = sum n_i^2 - sum_{i<j} n_i n_j`, and the `d = 2` remark by exact
multiplication in `Z[zeta_8] = Z[z]/(z^4+1)`.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    9868086c4b1568cd97bc7055af363a998efa240efbe52a4d909bddadb9efc6f1

## Scope

The nonclosure half of the refutation needs no program: Table 1, the three generators of
`Lambda`, and the three facts about `GF(27)` in Section 3 of the paper decide it by hand, and the
machine's own first failing sum is `Qm2 + Qm2`, i.e. `m = 2 = -1`, the very nonsquare the hand
argument names. The design and angle conditions, by contrast, are exhaustive finite checks over
142,506 characters and 265,356 pairs, and it is those the program carries out. The
program's closing statements of what it does **not** cover, quoted from its output:

> NOT RE-RUN: nothing here concerns d = 6 or d = 8, or any characteristic-2 dimension. The
> polar-form step 'polar-form-is-trace-pairing' needs char F odd, and no even-characteristic
> (Galois-ring / Kerdock) construction was built or tested.

> NOT RE-RUN: no MUB matrices are constructed. The program verifies the spread and
> quadratic-form data from which Theorem 4.4 of arXiv:2311.13479 and Abdukhalikov's construction
> produce a complete set of 28 MUBs in C^27; it does not multiply out the 28 bases.

> NOT RE-RUN: MUB-equivalence is untouched. The program does not decide whether some complete
> set of MUBs MUB-equivalent to this one has a subgroup (or coset) projective toric design; the
> MUB-equivalence group is a continuum and was not searched. Nor does it recoordinatise the
> spread by sending each of the 28 components to infinity in turn.

> NOT RE-RUN: no minimality is claimed or tested. Nothing here says d = 27 is the least
> dimension, or the least prime-power dimension, admitting a non-group projective toric 2-design
> of size d^2 that satisfies the angle condition; no census of other dimensions was run.

> NOT RE-RUN: the bridging dictionary between group toric designs and semifield spreads is used
> only in the direction needed here, and only for this X. The program verifies steps (i)-(iii)
> of Lemma 3.1 computationally on the printed object, not as a general theorem.

Two further limits belong to the paper rather than to the program. First, the **literal**
printed statement of Conjecture 4.9 also fails trivially to any translate of a group design
that misses the identity -- the paper exhibits the `d = 2` instance in Section 1 and the
program verifies it. What makes the `d = 27` object
load-bearing is that it contains the identity, so it also refutes the free repair "coset of a
subgroup". Second, Lemma 3.1 is **not new**: it is the known relative-difference-set /
planar-function / commutative-semifield dictionary of Godsil-Roy and Klappenecker-Roetteler,
with Kantor's published correction restricting it to *semifield* spreads, and the witness spread
together with its non-semifield property is printed by Abdukhalikov. The paper says so in
Remark 3.2. What is new is only the conclusion drawn for Conjecture 4.9.
