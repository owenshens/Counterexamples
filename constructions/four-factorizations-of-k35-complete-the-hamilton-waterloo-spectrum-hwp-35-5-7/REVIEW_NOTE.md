# Four Factorizations of K_35 Realizing beta in {1,2,4,6} in the Hamilton–Waterloo Problem HWP(35;5,7)

`four-factorizations-of-k35-complete-the-hamilton-waterloo-spectrum-hwp-35-5-7`

Supporting material for this paper: the program that re-derives its computational claims from
the objects printed in it, and a record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

Those four files are the whole deliverable and they are self-contained: `paper.tex` compiles with a
stock TeX distribution and `verify.py` needs only Python 3.9+ and the standard library. Nothing
here reads an external data file. The *Provenance* section below names further paths under
`runs/wave23/artifacts/t7941/` — those are internal search-and-archive records of this project,
**not files shipped in this folder**, and nothing in the paper or in `verify.py` depends on them.

## What is claimed

1. `(16,1)`, `(15,2)`, `(13,4)` and `(11,6)` lie in `HWP(35;5,7)` — Theorem 1 of the paper.
 Each is witnessed by one explicit factorization of `K_35` into 17 two-factors, printed in full in
 Section 4 as cycle lists on the vertices `0..34`.
2. Corollary 2, stated conditionally: *if* the theorem of Wang, Lu and Cao is as paraphrased in
 Section 1 of the paper from its e-print (reference [2] — quoted, not displayed, not re-proved
 and not verified by `verify.py`), then `HWP(35;5,7)` is completely determined, all 18 pairs
 with `alpha,beta >= 0` and `alpha+beta = 17` being realizable.

Nothing else. The lemma on the `p=7` lane that earlier drafts carried (a `rho`-invariant
`C_7`-factor is a union of five intra-class edge orbits, so the `p=7` lane cannot reach
`beta = 4, 5, 6`) has been **removed from the paper**: it motivated a design choice and was not
needed for Theorem 1. `verify.py` and its recorded run still refer to it as "Lemma 2" (check
`p7-lane-omits-4-and-6` and item (2) of the program's scope note); those two files are left
byte-for-byte as run, so that reference now has no counterpart in the paper. The check itself is
arithmetic about `beta = 7*A7 + B7` and no claim of the paper depends on it.

No priority is claimed for the four cells, and no literature check bearing on priority is
reported. Burgess–Danziger–Traetta's exception clause for `v = mn` names `beta in {1,2,4,6}` plus
the listed quadruple `(5,7,9,8)` at `v = 35`, and Wang–Lu–Cao proved the other 14 pairs, `(9,8)`
among them; falling outside the positive range of those two theorems is a fact about the two
theorems and not evidence that no other source realizes these cells. Whether the four cells had
been realized elsewhere is not settled here.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package, no external data file, no
network. All arithmetic is exact integer arithmetic on the vertex set `{0,...,34}`; there is no
floating point and no tolerance anywhere in the program. It prints one line per check and a
closing verdict, and exits 0 only if every check passes. The recorded run reports **66 checks,
all passing**:

 VERDICT: ALL 66 CHECKS PASS

The 28 generating objects embedded in the program are the cycle lists of Section 4 of the
paper, transcribed character for character; the program re-derives every quantity the paper
asserts about them, and nothing is taken from the search that produced them. Its five groups
of checks are:

- **A (8 checks)** — the two prescribed maps `rho` (order 7) and `sigma` (order 5) are
 permutations of the 35 vertices preserving the stated blocks, and their induced actions on
 `E(K_35)` have exactly the orbit census the paper's bookkeeping uses: 85 orbits of size 7
 (15 intra-class + 70 inter-class) and 119 of size 5 (14 + 105), with `85*7 = 119*5 = 595 =
 |E(K_35)|`, every intra-class orbit a `p`-cycle on its class and every inter-class orbit a
 perfect matching between two classes. Plus the necessity `alpha+beta = (35-1)/2 = 17` and the
 arithmetic consequence of the withdrawn `p=7`-lane lemma (`beta = 7*A7 + B7` with `B7 <= 3`
 misses 4, 5, 6), which the paper no longer states and no claim of the paper uses.
- **B (28 checks)** — one per printed object: the cycles are vertex-disjoint, span `0..34`, all
 have the declared length, and give 35 distinct edges.
- **C (24 checks)** — six per cell: every object declared invariant is fixed setwise by the lane
 map; every base object's orbit closes after exactly `p` steps into `p` distinct factors; the
 expansion has exactly 17 factors; every factor is 2-regular on all 35 vertices with cycle type
 `[5]*7` or `[7]*5` as declared (walked component by component, not assumed); the 17 factors are
 pairwise edge-disjoint with union **equal to `E(K_35)` as a set** (595 edges, set equality, not
 a count); and the counted `(alpha,beta)` is the claimed cell.
- **D (4 checks)** — a negative control on the checker itself: the same routine applied to a
 one-transposition tamper of each cell must refuse, and does (4 of 4).
- **E (2 checks)** — the arithmetic of assembling the spectrum: our four `beta` are exactly
 `{1,2,4,6}` and together with the other 14 values cover every `beta` in `0..17`.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 02be092132f7920b3d1585828cfd23a065fc96832fffda6265e7551da3493ab6

Run locally, single process, Python 3.9.6, under 0.2 s, no randomness. The recorded exit status
is 0.

Where the objects came from, as recorded in this row's artifacts manifest
(`runs/wave23/artifacts/t7941/MANIFEST.json`, 14 entries):

- The four certificates were found by five CP-SAT search scripts (`job2.py`, `job3a.py`,
 `job4a.py`, `job4b.py`, and `job3b.py` for the inherited `(9,8)` control) dispatched to
 32-vCPU AWS slots, CP-SAT with 30 workers, `python 3.9.25` and `ortools 9.15.6755`. The
 manifest records CP-SAT status `OPTIMAL` and, for four of the five, a wall time; it records
 **no** SSM `CommandId`, slot label or instance id for those five dispatches — those were not
 captured, and the manifest carries `null` there rather than a guess. **The stdout of the five
 search jobs was never captured to a file and no S3 URI for it exists**, so the only surviving
 record of those runs is the manifest's transcription plus the objects themselves.
 Consequently: the search is **not** reproducible from this folder, and CP-SAT with 30 workers
 is not deterministic in any case — a re-run may return a *different* factorization. What the
 claim rests on is the objects printed in the paper, and those are re-verified from scratch by
 `verify.py`.
- The archived artifacts also contain a second, longer certificate for `(16,1)` (11 objects, one
 base orbit) and a reconstruction of the published cell `(9,8)`. Neither is printed here: the
 first is redundant, and `(9,8)` belongs to Wang–Lu–Cao.
- `verify.py` in this folder is a fresh, self-contained program written for the paper's cycle
 presentation. It is not the artifact-side checker (`artifacts/t7941/verify_t7941.py`, which
 reads the same factorizations as edge lists); the two agree on all four cells.

## Scope

Quoted from the program's own closing statement:

> SCOPE NOTE -- what this program does and does not decide.
> DECIDED HERE, in full: that the four printed certificates are factorizations of K_35
> into 17 2-factors realizing (alpha,beta) = (16,1), (15,2), (13,4) and (11,6), i.e.
> that those four cells lie in HWP(35;5,7). Each is a finite, self-contained check on
> the objects printed in the paper.
> NOT DECIDED HERE (1): the fourteen remaining pairs with alpha+beta=17. Those are
> Wang-Lu-Cao's theorem, quoted, not re-proved; check E only does the arithmetic of
> putting the two halves together.
> NOT DECIDED HERE (2): Lemma 2 (a rho-invariant C_7-factor is a union of five
> intra-class orbits). Its hand proof is in the paper; this program checks only the
> finite orbit census the proof consumes (check A) and the arithmetic consequence.
> NOT DECIDED HERE (3): anything about which cells the PRINTED journal version of
> Burgess-Danziger-Traetta excepts. That is a bibliographic question; the completion
> above rests on Wang-Lu-Cao plus these four objects and does not depend on it.
> NOT CLAIMED ANYWHERE: uniqueness, enumeration, or minimality of these factorizations.

(The program's "Lemma 2" is the `p=7`-lane lemma, which has since been removed from the paper;
see "What is claimed" above. The program is left exactly as run.)

Further limits a referee should know, beyond what the program can see:

- **Bibliographic risk.** The exception clause of Burgess–Danziger–Traetta was read from the
 e-print `arXiv:1510.07079v2` (single source member `HW-submit2.tex`, 114,292 bytes, member
 timestamp 2015-11-20, theorem at line 1022) and Wang–Lu–Cao from `arXiv:1605.00818`
 (source member `ARCS.tex`, 66,252 bytes). **The typeset journal pages of neither paper were
 read**; a copy-editing difference in a printed exception list is therefore not excluded. It
 would change which cells one *describes* as open, not the validity of the four objects.
 Corollary 2 is stated with the Wang–Lu–Cao paraphrase as an explicit hypothesis, so a reader who
 does not grant that paraphrase loses the corollary and keeps Theorem 1.
- The same authors' companion paper (`doi 10.1002/jcd.21586`) is paywalled, has no e-print, and
 **was not read**. Its non-overlap with these cells is inferred from the 2023 survey's account
 that it excepts `v = lcm(m,n)` — an inference from a survey, not from the paper.
- *MathSciNet* was **not consulted** (no access); for a 2017 journal paper a review there is a
 plausible home for a prior remark on these cells.
- The prior-art search that did run covered nine channels (arXiv, zbMATH, OpenCitations,
 Semantic Scholar, Crossref, plus full-source reads of the two key e-prints) and found no
 source stating the completion of `HWP(35;5,7)` and none stating the withdrawn `p=7`-lane lemma.
 One channel **failed rather than answered**: the OpenAlex citer list for the
 Burgess–Danziger–Traetta work returned HTTP 429 and then a fetch deadline on two attempts. Not
 finding a source is not a demonstration that none exists, and the paper claims no priority.
- `(9,8)` is **not** claimed here. It is a lemma of Wang–Lu–Cao (`arXiv:1605.00818`,
 J. Combin. Des. 26 (2018) 27–47); our independent reconstruction of it was a calibration
 control and is credited entirely to them.
- No minimality claim, and none about any layer other than `v = 35`. The paper's earlier
 assertion that eight cells at `v = 33` are open has been removed as unsupported.
