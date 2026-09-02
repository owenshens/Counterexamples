# A Four-Vertex Counterexample to the Clow--van Bommel Eternal Distance-k Domination Conjecture

`a-four-vertex-counterexample-to-the-clow-van-bommel-eternal-distance-k-domination-conjecture`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

Clow and van Bommel (*Eternal distance-2 domination in trees*, Ars Combinatoria **163**
(2025) 29--49, doi `10.61091/ars163-03`, = `arXiv:2308.00054v1`) conjecture, as
**Conjecture 8.2** of their Conclusion:

> For all `k > 2` and trees `T`, if `gamma_k(T) = gamma_{all,k}^inf(T)`, then
> `gamma_{all,k}^inf(T) = gamma_{floor(k/2)}(T)`.

The paper refutes it with `T = P_4` at `k = 3`, where
`gamma_3 = gamma_{all,3}^inf = 1` while `gamma_{floor(3/2)} = gamma = 2`. Four vertices, four
distances, no program needed: `diam(P_4) = 3 = k` so a single guard distance-3 dominates and
can answer any attack, while no single vertex *dominates* `P_4`.

One further statement is proved, and it is a statement about infinitely many objects,
so a referee should read it as a separate claim:

- **Proposition 2.** For every odd `k >= 3`, *every* tree of diameter exactly `k` is a
 counterexample. The mechanism is parity:
 `2*floor(k/2)+1` is `k+1` for even `k` but only `k` for odd `k`. This is why `k = 2` is a
 theorem in the source and why extrapolating it is unsafe.


## What was checked, and how

The decisive argument is hand-sized and is carried out in the paper on the objects printed
there. `verify.py` re-derives it mechanically and independently:

- `P_4` is confirmed to be a tree with the printed
 order, edge set, diameter and radius (all-pairs BFS, integers only);
- every `gamma_j` the paper names is recomputed by an **exact minimum set cover** (iterative
 deepening on the budget, branching only on balls covering the lowest uncovered vertex);
 no bound is taken on faith, and each explicit dominating set the paper names is separately
 certified;
- every `gamma_{all,k}^inf` is recomputed **from the definition in Section 1**, as the
 greatest fixed point of "for every attack there is a live successor containing the attacked
 vertex" over the full multiset configuration space --- and in **both polarities**, so `g`
 guards are shown to win *and* `g-1` to lose;
- Proposition 2 is instantiated at `k = 3, 5, 7` and, as a negative control, the same
 construction is shown to yield nothing at `k = 4, 6, 8`;
- and, as a negative control, the same construction is shown to yield nothing at `k = 4, 6, 8`,
 which is bounded evidence for the surviving even-`k` case and nothing more.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only (`sys`, `itertools`): no third-party package and no
external data file. There is no floating-point value anywhere in the file and no decision is
taken on one; all arithmetic is exact integer arithmetic. Its inputs are exactly the objects
printed in the paper --- the adjacency list of `P_4`, the printed path labelling,
and the printed dominating set `{2,3}`. The program prints one line per check and a closing verdict, and exits 0 only if
every check passes. The recorded run reports **24 checks, all passing**, in about 0.1 s:

 VERDICT: ALL 24 CHECKS PASS

The recorded exit status is 0.

## Provenance

**Of this folder's program.** `verify.py` was written for this folder and run locally on this
machine. It is **not** one of the scripts that produced the result, and nothing in the paper
depends on it. `verify.output.txt` holds its output, preceded by a provenance header and
followed by an exit status, both written by the run harness. The header records the SHA-256 of
the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 94971758e92ba64cde6cc97a66346c4c45d8a6cbce8cadab5daf424c2254f4bc

`verify.py` reads nothing outside itself, so a referee needs neither this folder's other files
nor any of the material described next.

**Of the original computation, stated as the record states it.** The refutation used **no
compute at all**; a census was run *after* the verdict was decided, purely as corroboration.
That run's artifacts are recorded with per-file SHA-256 digests in an internal manifest.
**None of the files named in the rest of this section --- the manifest, `census2.py`,
`census2.stdout.txt` --- is shipped in this folder, and nothing in the paper or in `verify.py`
reads them;** they are described here only so that the corroborating run is on the record
rather than asserted. The four files in the table above are the whole folder. The manifest
records:

- one script, `census2.py` (9,303 B, sha256
 `8c37a827062b6a46ee7b7e7b5f725197b4ae48bc941791bc5e1103a2743391cb`), dispatched
 2026-08-31 as `aws/slot_run.sh --detach AUTO census2.py 900 angle1-t11158-evenk-census
 --force-decided` to one slot, EC2, SSM CommandId
 `(a dispatch id, redacted)`; STATUS Success, RC=0, elapsed 85.3 s, single
 process, pure Python 3 standard library, **no randomness and no seed**;
- its complete standard output, `census2.stdout.txt` (2,133 B, sha256
 `c14acc15c59fd8f9efa921169d277af8ebfc600a51bba1fd96aeb27242bf9dda`), recorded as complete
 (the wrapper reported 2,132 stdout bytes against a 2,133-byte file, an order of magnitude
 below the transport's ~24 KB cap, with no truncation marker), reporting five controls PASS,
 102 counterexamples at `k = 3` over all 2,288 trees of order at most 13, 416 at `k = 5`,
 and **0 counterexamples at each of `k = 4, 6, 8` over all 32,508 trees of order at most
 16**, with 0 undecided cells;
- the manifest's reproduction instruction, verbatim and in full: **`python3 census2.py`**.

**Honest gaps in that record, recorded rather than papered over.**

- The manifest states that `census2.py` was identified as the script **as run** from the
 recovered output itself, not from testimony: it is the only candidate in the attacking
 agent's scratch directory whose format strings reproduce that stdout verbatim, a sibling
 draft written two minutes earlier lacking three of the prints. So the code/output binding is
 an inference from the bytes, a good one, but an inference.
- The manifest's `reproduce` field says the script prints 27 lines; `wc -l` on the recorded
 stdout says 28. The integers in the output are what matters and they are unaffected.
- There is **no archived second copy** of that output: the wrapper recorded three
 `head_object 404` lines, so the internal artifacts directory is the only copy. That
 directory was filed by a later audit stage, not by the agent that dispatched the run.
- The recorded stdout prints only the three smallest witnesses per `k`, so `P_8` has **no
 counterpart in it**; the `P_8` claim rests on the hand argument of Proposition 3 and on
 `verify.py`, not on agreement with that run.
- No census of any kind is reproduced by `verify.py`; its own `k = 4` control covers order at
 most 9, not 16 (see Scope).

## Scope

What is refuted is Conjecture 8.2 **as written**, i.e. as a statement universally quantified
over `k > 2` and over trees; one odd `k` settles that, and `P_4` settles it at `k = 3`. What is
**not** settled:

- **the even-`k` restriction of Conjecture 8.2 remains open.** The mechanism of Proposition 2
 requires `k` odd: for even `k` a diameter-`k` tree has radius exactly `floor(k/2)` and the
 conclusion holds. The only evidence in the even case is a bounded census --- 0
 counterexamples among all 32,508 trees of order at most 16 at each of `k = 4, 6, 8`. That is
 bounded evidence, **not a proof**, and a referee should treat the even case as untouched.
- **no exhaustive classification at a fixed odd `k`.** Proposition 2 covers only the trees of
 diameter exactly `k`, which are not the full counterexample set: the internal census found
 102 counterexamples at `k = 3` among trees of order at most 13.
- **minimality only at `k = 3`.** `P_4` is of least order there because no tree on at most 3
 vertices has `gamma >= 2`; no minimality is claimed for other `k`, and the paper's own
 `P_{k+1}` is not asserted to be the least-order witness at odd `k > 3`.
- **no new mathematics.** The decisive ingredients are the source's own Lemma 3.3, the
 identity `rad = ceil(diam/2)` in a tree, and `gamma(P_4) = 2`. The contribution is the
 observation that these already contradict Conjecture 8.2, plus the parity account of why.
- **quoted results are not reproved:** the sandwich
 `gamma_k <= gamma_{all,k}^inf <= gamma_{floor(k/2)}` of Cox--Meger--Messinger and Lemma 3.3
 of the source. Both are quoted, and the paper says so at the point of use. `verify.py` does,
 however, recompute every `gamma_j` it needs by exact minimum set cover, so the numbers are
 corroborated on the finite ranges where they are used.
- **the other conjectures and questions of the source's Conclusion are untouched.**

The program's own closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE: Proposition 2 is stated for ALL odd k >= 3 and ALL trees of diameter exactly k; Step 3 instantiates it at k=3,5,7 on P_{k+1} and on one non-path tree. The general statement rests on the proof in the note, not on these instances. Nothing here bears on the even-k restriction of Conjecture 8.2, and no claim about the contents of the cited papers is checked by this program.

## One thing a referee should check against the source

Conjecture 8.2 is the **second** numbered conjecture of the Conclusion, immediately after the
polynomial-time one, in a `conjecture` environment sharing the `theorem` counter. In the
e-print it sits at lines 482--484 of `TheMainDocument.tex` (68,388 B) inside the tarball at
`arxiv.org/e-print/2308.00054`; the Conclusion is Section 8 in print and Section 7 in the
e-print, but the conjecture is numbered 8.2 in both and its sentence is verbatim unchanged in
the version of record, although the paper was otherwise revised. The paper reproduces the
sentence in full, so what is refuted does not depend on the labelling.

One transcription trap in the surrounding literature, which the paper avoids but a
referee checking by hand will hit: the e-print printing of the source's Theorem 2.3
**misprints its subscript** as `gamma_{all,2}^inf` against a `k` on the right-hand side; this
is corrected in print, and the primary source `arXiv:2104.03835` has subscript `k`.
