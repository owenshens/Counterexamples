# The Binomial Edge Jacobian of (K_3,K_3) Loses Rank in Characteristic Two

`binomial-edge-jacobian-rank-drops-by-one-in-characteristic-two-answering-landsittel-nevo`

Supporting material for this paper: the program that re-derives its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

Landsittel and Nevo, *Analytic spread via linear matroids* (arXiv:2607.07458v1), prove that
`rk H(G,H) = l(J_{G,H})` over any field of characteristic `0` or larger than `2^{|V(G)||V(H)|}`
(the statement labelled `propB` in their e-print source, quoted verbatim in Section 1 of the paper
with a byte-level locator), and ask in their Question `q:char` whether the characteristic
hypothesis can be removed.

The claim of this paper is that it cannot, and that characteristic `2` is the entire obstruction:

- at `G = H = K_3`, over **every** field, `l(J_{K_3,K_3}) = 9`, while `rk H(K_3,K_3) = 9` in every
 characteristic other than `2` and `= 8` in characteristic `2`.

That single pair is the whole claim of the paper. The value `l(J_{K_3,K_3}) = 9` is not the paper's:
the source itself proves `l(J_{K_n,K_m}) = nm` for `m,n >= 3` and remarks that its proof needs no
assumption on the field. A general triangle-covered family, statements about other pairs, and any
novelty claim were removed from the paper and are not asserted here either.

⛔ **What must not be taken away.** Nothing printed in the source is false. `propB` as printed
assumes `char F = 0` or `char F > 2^9 = 512` at this pair, and `F_2` lies outside that hypothesis.
What is refuted is the affirmative answer to an **open question** the authors themselves pose, not
a published proposition. The degree-two representability half of `q:char` is untouched: at
`(K_3,K_3)` over `F_2` the nine generators are algebraically independent, so their algebraic
matroid is free and therefore representable over every field. What fails is the Jacobian's ability
to compute it.

## What was checked and how

The decisive verification is carried out in the paper itself, on the objects printed there, and
needs no program at all. Table 1 of the paper is the `9x9` matrix `H(K_3,K_3)`; Table 2 is its
specialisation at `X = I`. Six of the nine rows of Table 2 carry a single `±1`, in six distinct
columns; deleting those six rows and columns leaves

 A = [[1,1,0],[1,0,1],[0,1,1]], det A = -2.

So `rk H(I) = 6 + rk A` over any field: `9` when `2` is invertible, and `8` in characteristic `2`,
where the three rows of `A` sum to zero. A specialisation cannot raise rank, so that is a lower
bound for the rank over the function field; the matching upper bound in characteristic `2` is the
Euler vector `v = (x_{i,k})`, which satisfies `H v = 2f = 0` because every generator is a quadric.
The lower bound `l = 9` is characteristic-free and comes from two classical identities over `Z`,
`adj(adj X) = det(X) X` and `det(adj X) = det(X)^2`, which put `F(x)^2` inside `F(g)`. That is the
whole argument, and it is a `3x3` determinant.

`verify.py` re-derives all of it from the printed edge lists rather than reading it, and adds:

- the polynomial identity `det H(K_3,K_3) = 2 det(X)^3` in `Z[x_11,...,x_33]` (54 monomials each
 side), its integer content `2`, its vanishing mod `2` and its non-vanishing mod
 `3, 5, 7, 11, 13, 509` — the last chosen below the source's own threshold `512 = 2^9`;
- the determinants `-1, 2, -3, 4, -5` of the map `E -> tr(E) I - E` on `n x n` matrices for
 `n = 2..6`. The line the program prints from these, that characteristic `2` degenerates only at
 `n = 3`, is **wrong as stated** — `2` divides `n - 1` for every odd `n` — and nothing in the paper
 depends on it; the paper draws only the case `n = 3`;
- seven printed pairs built from `K_3`, `K_4`, the diamond and the bowtie — corroboration only, corresponding to no claim of the paper — each by a
 two-sided squeeze: an evaluation supplies the lower bound and the Euler vector the upper, so
 every rank reported is exact rather than a bound;
- the non-chainable pair `(2K_3, 2K_3)`, where the matrix is verified to be a direct sum of four
 copies of `H(K_3,K_3)` and the characteristic-`2` deficiency is `4`, not `1` — chainability buys
 the exact value, not the refutation;
- the must-stay-silent control `(K_n, K_2)` for `n = 3..6`, the Theorem `thm-hyp` side of the
 source. The three trace-zero `2x2` matrices are written out as kernel vectors, verified to
 annihilate every generator as identities over `Z`, and shown independent mod `2`; with the
 evaluations this pins `rk H = 2n-3` in characteristic `0` and in characteristic `2` alike. The
 mechanism is silent exactly where it should be;
- one negative control, so the checker can return "no": the deliberately wrong identity
 `det H = 2 det(X)^2` is rejected.

**One correction to our own earlier draft, recorded rather than smoothed over.** An earlier
internal write-up of this result stated `det H(K_3,K_3) = -2 det(X)^3`. With the row and column
orders of Table 1 — the source's own ordering — the determinant is `+2 det(X)^3`; the sign is an
artefact of the chosen orders and is not an invariant, and `verify.py` demonstrates this by
swapping two rows and recovering `-2 det(X)^3`. Only the integer content `2` is used anywhere, and
it is ordering-free. The paper states the identity with the sign pinned to the printed ordering.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file. It
runs in under a second. The program prints one line per check and a closing verdict, and exits 0
only if every check passes. The recorded run reports **75 checks, all passing**:

 VERDICT: ALL 75 CHECKS PASS

All arithmetic is exact — `int`, `fractions.Fraction`, and a `GF(2^16)` whose modulus
`x^16 + x^12 + x^3 + x + 1` is proved irreducible inside the program by trial division against
every `GF(2)[x]` polynomial of degree 1 to 8. There is no floating point. There is no randomness:
the one search in the file (for an evaluation point of a corroborating pair) walks the fixed linear
congruential sequence `s_{t+1} = (1103515245 s_t + 12345) mod 2^31` from the written seed
`20260901`, so a re-run is byte-identical.

Its inputs are exactly the objects printed in the paper: the edge lists of `K_3`, `K_4`, the
diamond, the bowtie, `2K_3` and `K_2`; Table 1 cell by cell; the point `X = I`; and the identities
of Lemmas 2 and 3.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an exit
status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 e2ce2c8a389b23c012601dc42b0e38dc3a462eac19ce3c4b98a2ada28e6be1aa

`verify.py` was written for this folder and run locally (Python 3.9.6, exit status 0). It is **not**
one of the programs that produced the result; it is an independent re-derivation of the paper's
printed claims from the printed objects.

The result itself came out of a set of recorded jobs, transcribed below from the run's artifacts
manifest and not reconstructed. **None of them is shipped in this folder**, and nothing in the
paper depends on re-running any of them.

- The manifest's own summary of the headline is that it needs no computation: the certificate is a
 prose file, `angle2_THEOREM.md`, sha256 `671d51a904327aceafcde867b36cc94f1b2a80bd20b875f52bfa7d4b3f887851`,
 7,798 B, holding the separability dictionary, the `K = L(d)` construction and both polarities of
 the rank. The manifest states in its own words that "the unconditional refutation needs NO
 computation at all ... the computations corroborate it".
- The command the manifest records as the way to re-run the corroborating census is
 `python3 job7_t11440.py`. That script is `job7_t11440.py`, sha256
 `fa5088abbccbf41c261b0ddca6a816951a948580be29b5b6ba529ebf474cb2e3`, 11,155 B, dispatched as
 `aws/slot_run.sh --detach --force-decided AUTO job7_t11440.py 900 t11440-angle2-job7` on slot
 S24, CommandId `(a dispatch id, redacted)`,
 `STATUS Success`, `RC=0`, `vcpu=1`. Its captured stdout is `angle2_job7_OUTPUT.txt`, sha256
 `93321e60d00e833130003d59aa7fc9e9d8275b836b16f3a4a4ef3865e790a2c3`, 1,687 B, recorded complete.
- Six sibling jobs are indexed the same way, each with its own CommandId in the manifest's
 `invocation` field: `angle2_job1_t11440.py` (S06, `(a dispatch id, redacted)`),
 `angle2_job2_census_t11440.py` (S13, `(a dispatch id, redacted)`),
 `job3_t11440.py` (S29, `(a dispatch id, redacted)`),
 `job4_t11440.py` (S08, `(a dispatch id, redacted)`),
 `job5_t11440.py` (S07, `(a dispatch id, redacted)`) and
 `job6_t11440.py` (S42, `(a dispatch id, redacted)`), all `STATUS Success`, `RC=0`.
 Interpreters and solvers: CPython 3 with sympy 1.14.0 for exact-rank cross-checks; everything
 else hand-rolled exact integer, `GF(2^16)` and `GF(2^61-1)` linear algebra, no floating point.
- The run's artifact gate returned `ARTIFACT-GATE t11440 OK files=21 indexed=15 sha_ok=15 notes=0`,
 exit 0: 7 CODE, 7 OUTPUT and 1 CERTIFICATE indexed, every digest re-hashed.

Five limits of that record are stated by the manifest or by the run's own review, and are repeated
here rather than smoothed over.

1. **Six further scripts in the artifacts directory are deliberately not indexed**
 (`census_t11440.py`, `census2_t11440.py`, `census3_t11440.py`, `census4_t11440.py`,
 `charcensus_t11440.py`, `controls_t11440.py`). The manifest's `note` field says why: they were
 written by a different arm of the same investigation, and only the arm that ran a script can
 attest that its output is complete. Their outputs are therefore **not** part of the attested
 record, and nothing in this folder rests on them.
2. **One indexed output is truncated and says so.** `angle2_job2_OUTPUT_TRUNCATED_PREFIX_LOST.txt`
 (21,110 B, `complete: false`) lost its *front* to the SSM stdout cap; its object-store copy is
 absent (`head_object` 404 — the job never wrote it), so the lost prefix is not recoverable, and
 the lost counts were re-derived and re-printed by job 3. That channel is unread, not empty.
3. **One printed line of `job7`'s own output is a vacuous check, and it is disclosed rather than
 quoted.** The line `ARITHMETIC refuted 91 + holds 29 = 120 of 435 : False` comes from a source
 line that evaluates `len(ref) + len(hold2) == len(res)`, i.e. `120 == 435`; the label is wrong
 and the comparison is not the one it names. It signals nothing about the census either way, and
 no number in this folder comes from it.
4. **The censuses are not part of this paper's claims.** The run's censuses over 435 and 561 pairs
 are corroboration; a substantial fraction of their cells are compute-limited rather than
 negative (the run records 163 of 561 as INCOMPLETE and 82 as genuinely open), and their
 characteristic-`2` ranks are point-evaluation lower bounds except where they meet a proved
 upper bound. None of those numbers is quoted in the paper, and `verify.py` re-runs no census: it
 verifies seven printed pairs by a two-sided squeeze instead, which tests the conclusion rather
 than the search.
5. **Wall-clock seconds and instance types beyond `vcpu=1` were not captured** by the dispatch
 logs for most jobs and are not guessed here.

## Scope

- **Only the `propB` half of Question `q:char` is settled.** The Theorem `thm-hyp` half — the same
 question for `J_G = J_{G,K_2}` in `2n` variables — is **untouched**, and the authors' own remark
 at line 233 of their source, that they have no counterexample there, remains true after this
 note. Section 10 of `verify.py` is the control that shows why: on `(K_n, K_2)` the corank is
 already `3 = dim sl_2` in every characteristic, so the Euler direction is not binding and no drop
 occurs.
- **No odd-characteristic failure is exhibited anywhere,** and no pair other than `(K_3,K_3)` is
 claimed, in either direction.
- **`l(J_{G,H})` is not computed from its definition anywhere in this folder.** It enters only
 through the source's own equigenerated lemma `l = trdeg_F F(g)` (their `lem0`, no characteristic
 hypothesis, taken from Heinzer–Kim Proposition 4.8), and the value `9` at `(K_3,K_3)` is
 certified by the adjugate identities, not by a Rees-algebra or special-fibre computation.
- **The published integer `l(J_{K_n,K_2}) = 2n-3` is not recomputed.** The control verifies only
 that the rank is the same in characteristic `0` and in characteristic `2` there.
- **Three ingredients are other people's** and are cited in the paper as such: Ingleton's
 non-representable algebraic matroid over `Z/2Z` (p. 166, Example 15 — cited by the source
 itself), the inequality `rk H <= l` over an arbitrary field (Beecken–Mittmann–Saxena, Lemma 9,
 which is the source's own `lem-ineq-any_characteristic`), and the identity
 `det H(K_3,K_3) = ±2 det(X)^3`, which is the classical Hessian of `det_3`. No novelty claim is
 made for anything in the paper.
- **The spread half is the source's own, and the paper now says so in its abstract.** The source
 proves `l(J_{K_n,K_m}) = nm` for `m,n >= 3` and remarks that the proof needs no assumption on the
 field; at `n = m = 3` that is Lemma 2 of the paper. The paper reproves that case from the adjugate
 identities only because the source's proof routes through Dey–Ofir–Grussler (arXiv:2605.27682),
 whose field generality was not checked here. The value is the source's.
- **Prior-art search answered, with named gaps.** Nothing found states this result. The citation
 axis of the source was closed three times with controls: Semantic Scholar returned
 `citationCount = 0` on the paper record and an empty (not missing) citation list; OpenAlex work
 `W7167821243` gave `cited_by_count = 0`; OpenCitations returned `[]` against a positive control
 that returned 24,956 records. All 54 distinct cite keys of the source's `biblio.bib` (55,940 B)
 were enumerated and resolved. Three full LaTeX sources were read end to end — the direct
 predecessor arXiv:2510.05915 (zero occurrences of "characteristic", "Jacobian", "matroid",
 "inseparab"), Tyomkyn arXiv:2602.11892, and Dey–Ofir–Grussler arXiv:2605.27682. arXiv metadata
 sweeps and ten zbMATH searches answered with their result counts read each time. The source's own
 suppressed Macaulay2 experiment (inside an `\iffalse` block) tests a characteristic-`2` criterion
 at `(K_5,K_2)` and reports that the two ranks *coincide* — a pair structurally outside the
 triangle-covered family, so there is no self-scoop. The gaps that remain are gaps, not
 clearances:
 - **MathSciNet was not consulted at all.** For a 2026 preprint a review there is a plausible
 place for a prior observation.
 - **arXiv full-text search does not exist.** The arXiv API searches title, abstract, authors and
 comments only, so a result buried in the body of a paper whose abstract avoids these phrases
 would be invisible to every sweep run. This is the honest residual limit.
 - **The OpenAlex `cites:` filter failed** (HTTP 429, then a fetch deadline). The same question
 is closed from the other side three times over, but that channel is unread.
 - **One zbMATH query returned HTTP 404 and 0 B** (`vanishing Hessian positive characteristic`)
 and is recorded as unread in either direction.
 - **Cross-result duplication was not checked.** No clustering was run against the other results
 of the same batch, so nothing here asserts that this result is unique within it.
