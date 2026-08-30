# A Stabilizer Code Admitting No Logical Z That Meets the Hypothesis of Zheng--Liu Proposition 3

`the-smallest-magic-state-code-refutes-the-zheng-liu-fixed-plane-universality-claim`

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

Python 3.9 or later, standard library only: no third-party package, no external data file, no
network. Every stabilizer group, logical frame, spoiling product and coset representative it
consumes is a literal in the program's own source; those relevant to this paper are copied from
`paper.tex`, and the rest belong to the longer draft described below.
All arithmetic is exact integer or `Fraction` arithmetic; no decision is taken on a float. It
prints one line per check and a closing verdict, and exits 0 only if every check passes. The
recorded run reports **211 checks, all passing**:

    VERDICT: ALL 211 CHECKS PASS

It re-derives four named witness codes, of which only the first, the `[[3,1,1]]` code, appears in
this paper (abelian, independent, printed element list,
`|N(S)| = 2^(n+k)`, `4^k` logical classes, each printed spoiling product closing as Pauli
arithmetic with a multiplier genuinely in `S` and a Z-free right-hand side, the representatives
pairwise in distinct cosets and exhausting the classes) and then confirms the clean-class counts a
second time by an exhaustive sweep of `N(S)`, independently of the printed certificates. It also
carries checks on material this paper does not state: the other three witness codes, a `k=2`
certificate, the block family `S(n,k)` for every member with `2n <= 16`, an exhaustive census of the
six cells with `n <= 4` including a clean-class histogram, a counting bound, and positive controls.
Its labels and commentary name sections, lemmas, propositions and tables of a longer draft, not of
this paper. The paper claims only the `[[3,1,1]]` certificate and the coset lemma; the remaining
checks are not claims of the paper. Runtime is about
four seconds on one core.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    666f544ee013410f2ee38c6d290f0161e8a2d3f7d984a824006a26c75ef7d2d8

## Scope

The program's closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: the census here is exhaustive for n <= 4 only (all 17478 stabilizer groups with
> 1 <= dim S <= n-1). The n = 5 cells and every n >= 6 row quoted in the discovery record are NOT
> recomputed here; a miss there would be inconclusive rather than negative.

> NOT RE-RUN: the target paper [[15,1,3]] and [[14,2,2]] protocols, and its [[6,1,2]] code, are not
> decided here -- only its [[3,1,1]], [[4,1,1]], [[5,1,3]] and [[7,1,3]] are.

> NOT RE-RUN: nothing here computes the dynamical map of the target paper. The claim checked is
> exactly that Proposition 3 hypothesis is unsatisfiable; whether z = 0 is nonetheless an invariant
> plane for these codes is a different question, and the y-direction checks above show these codes
> DO have an invariant coordinate plane.

> NOT RE-RUN: the lambda = (3/2)^n / 2^k rule is checked only for sign agreement on the six cells
> above. It is an empirical heuristic fitted to census data and no proof of it is claimed or
> attempted.

Two further limits belong to the paper rather than to the program, and the paper states both. The
witness code is not new: it is imported by the target paper from Howard and Dawkins (2016). What is
new here is the statement about it, never the object. And
the target's Proposition 3 is a sufficient
condition, stated with ``if''; it is correct and is not touched --- only the separate, unnumbered
sentence asserting that its hypothesis can always be arranged is refuted.
