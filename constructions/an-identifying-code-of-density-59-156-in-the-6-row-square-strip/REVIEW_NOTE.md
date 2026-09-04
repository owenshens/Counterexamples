# An Identifying Code of Density 59/156 in the 6-Row Square Strip

`an-identifying-code-of-density-59-156-in-the-6-row-square-strip`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

`d*(S_6) <= 59/156`, where `d*` is the minimum density of an identifying code (radius 1) of the 6-row
infinite square strip `S_k = Z x {0,...,k-1}`. **The paper claims the upper bound only**; the matching
lower bound is neither claimed nor reported in the paper, and its witness was never filed (see below). Theorem 30 of
the Lobstein--Hudry--Charon survey gives exact values for `k <= 5` and, for `k >= 6`, only the
bracket `7/20 + 1/(20k) <= d* <= 7/20 + 3/(10k)`, which at `k = 6` is `[43/120, 2/5]` of width `1/24`.
The paper quotes part (f) of Theorem 30 with its locator (hal-02916929v1, Section 6.4.1.1, PDF
page 11 = manuscript page 10). That quote was re-checked against the HAL PDF at bundle-review time:
the quoted text matches the source symbol for symbol, including the chapter's density symbol `∂`.

**The two halves are not of equal standing, and the paper is written to keep them apart.**

* **Upper half, `d* <= 59/156`.** Complete and self-contained. The witness is an explicit period-52
 code, printed in the paper twice (52 six-bit column masks, and the same object as six 52-character
 row strings). Verifying it from the definition is a finite check: 312 domination tests and 1612 pair
 tests, because two vertices at distance `>= 3` have disjoint closed neighbourhoods and therefore
 distinct identifiers as soon as both are nonempty. No solver is involved, and the resulting strict
 inequality `59/156 < 2/5` already shows the published upper endpoint is not tight at `k = 6`.
* **Lower half, `d* >= 59/156`: not claimed, and no longer even reported in the paper.** The
 computation behind it produced an integer feasible potential `q` on Jiang's transfer
 digraph `H_6` satisfying `52*w(u,v) - 118 + q(v) - q(u) >= 0` on every edge, which a telescoping
 argument would convert into `lambda(H_6) >= 59/26`, hence `d* >= 59/156`,
 and hence the equality; the revised paper asserts none of this. **The 8,072,634-entry vector `q` is
 printed nowhere and is not in this folder or in the row's artifact record, and the standard output of
 the run that produced it was never captured.** A referee cannot check the lower half from this
 folder; they would have to recompute it, and the paper no longer relies on it.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only (`fractions`): no third-party package and no external data
file. It reads the 52 masks and six row strings printed in Section 2 of the paper and re-derives the
quantities involved. It prints one
line per check and exits 0 only if every check passes. The recorded run reports **28 checks, all
passing**:

 VERDICT: ALL 28 CHECKS PASS

Runtime is about 0.07 s. The three parts of the output are: (1) the two encodings of the period-52
tile agree row by row; (2) the tile is an identifying code --- 0 domination failures over 312 vertices,
0 separation failures over 8372 pairs at L1 distance `<= 6`, and all 6760 swept pairs at distance
3 to 6 do have disjoint closed neighbourhoods, so that part of the locality argument is checked and
the general statement at every distance `>= 3` is taken from the note; (3) the bracket endpoints
`43/120` and `2/5` and the gaps `17/780` and `31/1560`.

Nothing floating-point decides anything: all densities and comparisons are `Fraction`, and the two
`%.12f` prints are display only. There is no randomness.

I also confirmed, outside the program, that the 52 masks and the six row strings in
`paper.tex` are character-identical to the literals in `verify.py`, so the program really reads the
paper's object.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an exit
status, both written by the run harness. The header records the SHA-256 of the program that produced
the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 d7a888e01fa2f7d6c51dddb50fa10daf9d7c324929cad879a3f920ea6b937af2

That run was executed for this folder on the local control plane (macOS, arm64) under Python 3.9.6,
`python3 verify.py`, exit status 0.

**`verify.py` is new code written for this folder, not a copy of the scripts that produced the
result.** The provenance of the underlying computation, as recorded in the row's artifacts manifest
(the run record's `MANIFEST.json`, 13 files, sha256 per file):

* the upper-half check (`verify_tile2.py`, sha256 `75d7590900c29567d3918b7da7319b24068d9f052f071ae2bebf99dcccf5ab88`),
 the period-13 witness and minimum-period core (`final.py`), the symmetry scan (`sym.py`), the tile
 equivalence (`cmp.py`) and the encoding cross-check (`xcheck.py`) were run locally as
 `python3 <file>`, each in under 0.08 s, with their outputs filed alongside;
* the lower-half computation was the run record's `verify6.py` (sha256
 `11b4a5e1abb2529d7d4d8c1b11bf149cc9f528dd0373bd95a7bcaedd2f8baff5`), dispatched by the fleet
 dispatch script in detached mode, with the slot chosen automatically and a 3600 s timeout, onto
 one fleet slot, the recorded dispatch id `(a dispatch id, redacted)`,
 231 s, Success with RC=0 verified, 11,448 B of stdout; and the optimal-core follow-up was `core.py`
 on a second, different slot, recorded dispatch id `(a dispatch id, redacted)`, 172 s,
 4,370 B of stdout.

Three gaps in that record, stated rather than papered over. (i) **Neither slot stdout was written to a
file**; the manifest records the byte counts but not the bytes, the compute instances self-terminate,
and the session's cloud credential has expired, so the lower-half violation counts survive only as
transcriptions. (ii) The instance *types* were not recorded and are not guessed. (iii) Of three
independent lower-bound computations, only one is re-runnable from what was filed: a second attack
angle's `k = 6` script (recorded dispatch id `(a dispatch id, redacted)`, exit 0) and a third
implementation were never filed, because that dispatch omitted the artifact-filing flags. Nothing in
this folder depends on any of the three.

One correction I made while building this folder, worth flagging to a reader of the row's own
document: that document states "6,136 pairs at distance 3 and 4". Recomputed here, the number of such
pairs is **3,068** unordered (6,136 is the same set counted as ordered pairs), while the same
document's companion figure 8,372 is an unordered count. The paper and the program print 3,068 and
8,372, which are consistent with each other. No claim of the result depends on either figure.

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE: the checks above establish that the hard-coded period-52 masks have density 59/156,
> dominate all 312 vertices in one fundamental domain, and separate all 8,372 swept pairs at L1
> distance at most 6. Every swept pair at distance 3 through 6 has disjoint closed neighbourhoods.
> NOT ESTABLISHED HERE: disjointness or separation at L1 distance greater than 6, the exact value of
> d*(S_6), any lower bound, the optimality of the printed tile, and the cases k >= 7.

Beyond that, and beyond the three provenance gaps above:

* **`k = 7` is open**, and so is the exact value at `k = 6`; the paper says so, and claims only the
 upper bound `d*(S_6) <= 59/156`, with no claim that it is close to optimal.
* **Theorem 30(f) posed no question.** It is a proved theorem with a gap, so the paper does not present
 itself as answering anything; the revised title and abstract claim only the upper bound.
* **The result contradicts nothing.** `59/156` lies strictly inside the theorem-backed bracket
 `[43/120, 2/5]`. It narrows that interval from above; it is not a refutation.
* **The method is Jiang's, not ours.** The periodicity theorem, the transfer digraph `H_k` and the
 reduction to a minimum cycle mean are all from arXiv:1607.03848v2; both bracket endpoints belong to
 Bouznif--Havet--Preissmann and to Daniel--Gravier--Moncel. What the paper offers is the period-52
 code and the bound it gives.
* Every density in the paper is an identifying-code density, not a locating-dominating one.
* The revised paper drops material the reviewers judged excessive for a short note: the HAL byte count,
 the full census arithmetic and the three unfiled rebuilds, the tight-edge pass, the `k = 7` timeout,
 the stabiliser experiment, the period-13 tile, the interpolation discussion, and the two
 bibliographic remarks. Some checks in the recorded run therefore concern quantities the paper no
 longer prints; the program and its transcript were left untouched so that their hashes still pair.
