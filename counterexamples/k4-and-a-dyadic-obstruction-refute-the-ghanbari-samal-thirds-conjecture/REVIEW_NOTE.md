# Review note

Folder: `paper.tex`, `paper.pdf`, `verify.py`, `verify.output.txt`, and this note. Nothing outside the
folder is needed to read the paper or to re-run its computation.

## 1. What the paper claims

The paper refutes Conjecture 1 of Ghanbari and Šámal, *Facial diagrams and cycle double cover*
(arXiv:2605.01410v1), quoted verbatim in §1 from `main.tex` lines 306–313: in a random embedding
`(π,λ)` of a bridgeless cubic graph with `m` edges, the expected numbers of bad singular, good
singular and regular edges are each exactly `m/3`. §1 records that the source defines no probability
measure, so the paper fixes one: the uniform measure on the `2^{n+m}` signed rotation systems
`Ω = {(π,λ)}`. Under that measure all three parts fail, by two independent arguments.

* **Theorem 1** (§2, "A dyadic obstruction"): for a finite connected edge-transitive cubic graph,
  `Aut(G)`-equivariance forces the per-edge count `A_{e,i}` to be independent of `e`, so
  `E[#class i] = m·A_i/N` with `N = 2^{n+m}`, and `E[#class i] = m/3` would require `3A_i = 2^{n+m}`.
  Each of the three parts therefore fails for every such graph. §2 also bounds the deviation from
  below by `m/(3N)`, which on `K_4` is `1/512`, improved to `1/32` using Lemma 2.
* **Theorem 3** and **Corollary 4** (§3, "The `K_4` computation"): an exact census of the
  `2^{10} = 1024` signed rotation systems of `K_4` gives, per edge, bad singular in 288 systems, good
  singular in 368 and regular in 368, hence `P[bad] = 9/32`, `P[good] = P[reg] = 23/64`, and
  `E[#bad] = 27/16`, `E[#good] = E[#reg] = 69/32` against `m/3 = 2`. Since `2` is not a member of that
  multiset, no relabelling of the three classes rescues any part, and `K_4` is a simple, 3-connected,
  bridgeless cubic graph.

Two auxiliary statements support this. **Lemma 2** (§3) shows switching (invert `π_v`, flip `λ` on the
three edges at `v`) is a free `Z_2^V` action on `Ω` preserving all three edge classes, so averaging
over `Ω` equals averaging the `2^m` signatures at one fixed rotation, exactly; the census is done at
the single rotation `π_0` displayed in §3, and the class column sums over the 64 signatures are 108,
138, 138. **Theorem 5** (§4, "A sign-flip involution") shows that flipping `λ(e)` exchanges
`{e regular}` with `{e good singular}` and fixes `{e bad singular}`, changing the face count by `−1`,
`+1`, `0`; hence `P[good] = P[reg]` exactly and `E[#bad] = m − 2E[#good]`, so parts (2) and (3) are
each equivalent to part (1). §4 states that neither Theorem 1 nor §3 uses this.

§5 ("What is, and is not, settled") adds a second measure: uniform on the `2^n` orientable (pure
rotation) systems, where on `K_4` the eight coboundary signatures give `E[#bad] = 9/4`,
`E[#good] = 0`, `E[#reg] = 15/4`, again all different from `m/3 = 2`, with the sign of the deviation
flipped. §5 also records that the conjecture's only consumer in the source, its Theorem 4.1
(`main.tex` lines 317–319), keeps its conclusion, which the source proves unconditionally in the
referee remark at lines 324–325; what is lost is a conjecture and one redundant conditional proof.

## 2. What the program checks

`verify.py` (Python 3.9+, standard library only, exact integer and `Fraction` arithmetic, no
floating-point decision) implements the face tracing of the paper's equation (1) from the definitions.
The recorded run `verify.output.txt` ends `VERDICT: ALL 133 CHECKS PASS` with
`program exited with status 0`, and prints one `PASS <name> [detail]` line per check. The blocks, and
the claims they bear on:

* **8 `carried_table_*` checks** — the 64-row `K_4` table (signature, faces, Euler genus, six-letter
  class string, `(bad,good,reg)` triple) that §3 says is carried in `verify.py` rather than printed:
  64 rows, 64 distinct signatures, class string agreeing with its triple, `bad+good+reg = 6` on every
  row, `f = 4 − Euler genus`, column sums `108/138/138`, total `384`, and
  `carried_table_reproduced_row_by_row_from_the_definitions` (0 mismatches against the program's own
  trace of the 24 arrival states). These sums are the input to the proof of Theorem 3.
* **21 `K4_*` checks** — Theorem 3, Corollary 4 and the orientable paragraph of §5: 1024 embeddings
  traced (100% of `Ω = 2^10`); `E[#bad] = 27/16`, `E[#good] = E[#reg] = 69/32` versus `m/3 = 2`;
  `m/3` not a member of the three expectations; per-edge `P[bad] = 9/32` on all six edges; the three
  expectations summing to `m = 6`; `E[#singular] = 123/32`; the eight coboundary signatures recomputed
  from the 16 vertex subsets, with good count `0` and even Euler genus on all eight; and
  `E[bad] = 9/4 > 2`, `E[good] = 0`, `E[reg] = 15/4` for the orientable model, tagged in the run as
  the point where "the sign of the deviation FLIPS".
* **15 or 16 checks for each of six further cells** — `theta` (the 3-dipole), `K_{3,3}` and the 3-prism
  at 100% of `Ω` (32, 32768, 32768 embeddings), and `Q_3`, Wagner `V_8`, Petersen at all `2^m`
  signatures at each of two fixed rotations (8192, 8192, 65536 embeddings), which the run marks
  "exact by Lemma 2, tested below". On every cell: `bad+good+reg = m` on every embedding, the Euler
  and orientability parity condition, the three expectations all different from `m/3`, `m/3` not a
  member of them, and `P[good] = P[reg]` on every edge. The closing
  `all_seven_cells_have_all_three_expectations_different_from_m_over_3` records 7 of 7 cells, 3 of 3
  expectations, 0 exceptions.
* **Theorem 1's arithmetic instance** — `..._dyadic_obstruction_arithmetic` on the five cells with a
  single per-edge value (`theta` 1/4, `K_4` 9/32, `K_{3,3}` 39/128, `Q_3` 313/1024, Petersen
  1263/4096), each reported with its power-of-two denominator; for the prism and `V_8` the run instead
  records two distinct per-edge values and that the dyadic theorem does not reach them.
* **Lemma 2** — `..._per_rotation_counters_identical` on each cell (e.g. 16 rotation systems of `K_4`,
  all equal), and `switching_invariance_lemma_2_of_the_paper`: 23616 single-vertex switchings over all
  seven cells, 0 class changes.
* **Theorem 5** — `theorem_4_involution_and_face_count_delta`: 6240 (embedding, edge) flips on `theta`
  and `K_4`, 0 violations of regular↔good, bad fixed, and `Δf = −1/+1/0`.
* **Controls and cross-checks** — the source's assertion at `main.tex` line 119 (good count `0` on
  every orientable embedding) on all seven cells; `forced_positive_control_planar_tetrahedron_is_all_regular`
  (`f=4`, Euler genus 0, classes `RRRRRR`, and the `m/3` detector fires); `opposite_polarity` (880 of
  the 1024 `K_4` systems carry a good singular edge, so no negative is vacuous); 16 of 1024 `K_4`
  systems and 0 of 32768 `K_{3,3}` systems of Euler genus 0; and the two conditions §5 says are
  transcribed from the literature and used nowhere in the paper — the Bender–Richmond bound on all 64
  rows, attained in all four surface classes, aggregate `246 ≤ 259`, and the Gross–Furst genus
  distribution of `K_4` (`1 : 7` at `π_0`, `16 : 112` over the full orientable space).

## 3. What the program does not check

The run states its own limits in a closing `NOTE SCOPE` / `NOT RE-RUN` block; they are carried over
here together with the limits §5 and §6 of the paper state.

* **Theorem 1 is a hand proof, and the program is a control on it.** The program verifies only the
  arithmetic instance on the five single-orbit cells above; the general statement over all finite
  connected edge-transitive cubic graphs is not verified by computation. **Lemma 2** and
  **Theorem 5** are likewise hand proofs (for all graphs, resp. all loopless graphs); the program
  tests them, exhaustively only on `theta` and `K_4` — the switching test covers 2 rotations × 256
  signatures × every vertex on the other five cells, and the involution test runs on `theta` and
  `K_4` only.
* **Coverage of graphs.** Any bridgeless cubic graph outside the seven cells is NOT RE-RUN. The paper
  does not claim the conjecture fails for every bridgeless cubic graph; §5 notes that, being
  universally quantified, it needs only one counterexample.
* **The fixed-rotation cells.** For `Q_3`, `V_8` and Petersen only two rotations are traced; that this
  is exact rests on Lemma 2, which the program tests rather than proves.
* **The reading of the conjecture.** The run says in the detail of
  `K4_all_three_expectations_differ_from_m_over_3` and of the seven-cell check that the identification
  of its three computed expectations with the three parts of the conjecture is assumed and NOT
  VERIFIED, and that simplicity and 3-connectivity of the hard-coded `K_4` edge list are NOT CHECKED.
  §5 says the same and leaves both to the reader.
* **Transcribed numbers.** The Bender–Richmond bound is transcribed from the restatement by the
  authors of the source (their companion e-print arXiv:2511.07285, `CDC.tex` line 133) and not
  recomputed from the 1990 paper; the Gross–Furst `2 : 14` distribution is transcribed and then
  reproduced, not proved. §5 states that neither the 1990 nor the 1987 paper was consulted beyond its
  published abstract, that these two checks test the tracer against restatements only, and that no
  claim is made about either.
* **Quoted, not verified.** The infinitude of the connected cubic arc-transitive family is cited to
  Biggs and Hoare and, per the run and §5, not verified. Likewise no file in this folder verifies the
  quotations from and line numbers of the source's `main.tex`; they are transcriptions, and the checks
  named `..._main_tex_L119_...` test the *content* of that line against the tracer, not the quotation.
* **Multigraph cell.** `theta` is the 3-dipole; §6 says it is a cell of this conjecture only if
  multigraphs count as bridgeless cubic, a reading the source neither grants nor denies, and that
  nothing in the paper depends on it.
* **Measure.** Only the two uniform models named in §5 are settled; no non-uniform or limiting measure
  is covered. §5 explicitly does not address whether `E[#bad]/m → 1/3` along some family, nor whether
  `E[#bad] < m/3` in the uniform signed-rotation model. The unproved Properties 6–9 of the source are
  used nowhere in the paper.

## 4. How to check it

From this folder:

```
shasum -a 256 verify.py
python3 verify.py
```

The digest must be

```
f1de8ebb6b5eb20a7f0fa6898a7de299ce9d207e62d83121c95dd3d15ab46044
```

which is the value recorded in the header of `verify.output.txt` beside the program name, so the
shipped transcript and the shipped program can be paired. The second command reprints the `PASS`
lines, the `NOTE SCOPE` block and `VERDICT: ALL 133 CHECKS PASS`; the recorded run was made with
Python 3.9.25 and exited 0. There is no input file and no network use, and §3 notes that any single
row of the 64-row `K_4` table can be redone by hand from equation (1).

## 5. One naming discrepancy

The involution check is named `theorem_4_involution_and_face_count_delta`, and the per-cell
`..._P_good_equals_P_reg_on_every_edge` checks are annotated "(Theorem 4)"; the statement they test is
numbered **Theorem 5** in `paper.tex` and `paper.pdf` (Theorem 1, Lemma 2, Theorem 3, Corollary 4,
Theorem 5). The mathematics tested — `P[good] = P[reg]` per edge and `Δf = −1/+1/0` — matches
Theorem 5 and equation (2) as shipped.
