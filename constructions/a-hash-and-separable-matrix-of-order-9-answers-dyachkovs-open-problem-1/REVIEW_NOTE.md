# A Hash-and-Separable Matrix of Order 9: D'yachkov's Open Problem 1 Answers Yes

`a-hash-and-separable-matrix-of-order-9-answers-dyachkovs-open-problem-1`

Supporting material for this paper: the program that re-derives its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

That the `9 x 9` matrix `M` printed in Section 2 of the paper is a `C_HS(9,4)`-matrix in the sense
of A. G. D'yachkov, *Lectures on Designing Screening Experiments* — that is, it satisfies
Propositions 5, 6 and 7 of that text — and that therefore

 q_4^HS <= 9 < 13,

which answers **yes** to Open problem 1 of the same text: *"Is it possible to construct a
`C_HS(q,4)`-matrix if `q < 13`?"* The question is quoted verbatim in Section 1 of the paper with a
byte-level locator (`lect_dse.tex` inside the arXiv:1401.7505v1 e-print, lines 3790–3793, bytes
114,506–114,687), as are the three propositions the matrix has to satisfy.

The published value `13` that the paper improves on is a **reported construction, not a proved
minimum**: Open problem 1 sits one paragraph after the sentence that reports it, and nothing in
the source proves `q_4^HS = 13`. This note corrects a reported value; it identifies no error of
proof.

## What was checked and how

The decisive object is printed in full, so a referee who wants to check nothing else can check the
matrix directly. Two of the three conditions are genuinely hand-sized:

- **Proposition 5** is read off the grid: four symbols and five stars in every row and in every
 column, no repeat within a row or a column, each of the nine symbols in exactly four cells.
- **Proposition 6** reduces, by Lemma 2 of the paper, to nine `4 x 4` blocks — for each symbol
 `a`, the block on `a`'s four rows and four columns must be a permutation matrix over `{a, *}`.
 The row/column table needed to locate the nine blocks is printed in Section 2.
- **Proposition 7** is the one finite check that is not hand-sized. Lemma 3 of the paper removes
 the star conditions (they follow from Proposition 6) and leaves, for `M`, exactly **162**
 candidate configurations over the 84 triples of rows; each fails one of three column equations.

`verify.py` does all of that from the printed object, plus a second, independent route through the
source's raw Definitions 2, 4 and 5 (the induced code of 36 codewords: homogeneity, pairwise
Hamming distance at least 2, the 3-hash condition on all 7140 triples, and the 3-separability
condition on all 7806 subsets of size at most 3), and two brute-force encodings of Proposition 7
that use no reduction at all (84 x 84 increasing triples against all six printed forms, and
504 x 504 ordered triples against the single base form). All routes agree.

Both polarities of every checker are exercised on **published** objects before any of them is
believed, so a checker that could only say "yes" would be caught. The controls are the eight
labelled matrices of Examples 7–10 of the source, transcribed into `verify.py` with the source's
own labels. The load-bearing one is at the same `k = 4` as the witness and is printed in the
paper: D'yachkov's own `C_H(8,4)` passes Propositions 5 and 6 and **fails** Proposition 7, with 32
forbidden submatrices and, independently, 32 colliding unions under Definition 4. The published
`C_HS(13,4)`, `C_HS(7,3)` and `C_HS(4,2)` pass every check; `C_H(6,3)` and `C_H(3,2)` fail
Proposition 7; `C_S(3,2)` fails Proposition 6 and the 3-hash condition; the plain `C(6,4)` of
Example 7 fails Proposition 6 with 68 violations.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file. It
runs in under a second. The program prints one line per check and a closing verdict, and exits 0
only if every check passes. The recorded run reports **50 checks, all passing**:

 VERDICT: ALL 50 CHECKS PASS

Its only inputs are strings printed in the paper — the matrix `M`, the symbol/row/column table of
Section 2, and the `C_H(8,4)` of the controls remark — together with the other seven published
control matrices, which are transcribed into the program itself. Every quantity the paper asserts
is re-derived, including the counts `36 = 4*9`, `18/18/54`, `162`, `2304`, `42,336`, `7140` and
`7806`, the bounds `q >= k^2 - k + 1` and `q >= 2k - 1`, and the identity
`13 = 4^2 - 4 + 1 = |PG(2,3)|`. No floating point is used anywhere; every comparison is over
`int`.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an exit
status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 69d8030d84c9e134bc84e772b93ea06006c75c1fb697edaa3adb413a8c49e410

`verify.py` was written for this folder and run locally (Python 3.9.6, exit status 0). It is **not**
one of the programs that produced the result; it is an independent re-derivation of the paper's
printed claims from the printed object.

The result itself came out of four recorded jobs, transcribed here from the run's artifacts
manifest and not reconstructed. None of them is shipped in this folder, and nothing in the paper
depends on re-running any of them, because the witness is printed in full and Propositions 5–7 are
three paragraphs long.

- `sweep_q4_12.py`, sha256 `ad3022ac1f82c3ea786cac18d1a08920d85cf8b69b4a53374018c6b422391868`,
 9,976 B — the SAT encoding of Propositions 5/6/7 and the sweep over `q = 4..12` that returned
 `M`. Dispatched detached through the fleet dispatch script onto an automatically assigned slot,
 running the run record's `sweep_q4_12.py` under a 900 s timeout, a separate cloud instance,
 `RC=0`, about 20 s of solving in total.
- `controls.py`, sha256 `5a7dfff8c24e684f4af853c493b3cd8263e43c76f19eebeb6b9014a567325026`,
 13,674 B — the both-polarity control suite on the source's published matrices. a separate cloud instance, `RC=0`.
- `defcheck.py`, sha256 `7a8bf33813a19d68f16805e744fa2ae81ae196737828290f8d3916308d42cb38`,
 13,725 B — the raw Definition-2/4/5 route. a separate cloud instance, `RC=0`.
- `weak_break_q8.py`, 11,692 B — the re-run of the contested `q = 7` and `q = 8` cells under a
 weakened symmetry break.

Three limits of that record are stated by the manifest and are repeated here rather than smoothed
over. (i) All four scripts were dispatched from paths under a single scratch directory and copied
byte-for-byte into the run's artifacts directory afterwards; the invocation lines above are what
the manifest records. (ii) The exact cloud instance type, the runner image's `python3` version and the
`python-sat` build version **were not captured**; the solvers used were `Cadical153`, `Glucose4`
and `Minisat22` in default configuration with **no seed set**, so a different build may return a
*different* valid `C_HS(9,4)`-matrix — which is why `M` is recorded as a printed object rather
than left to be re-derived. (iii) The audit stage's independent re-implementation of Propositions
5, 6 and 7, written from the source's LaTeX without reading ours, confirmed all three on `M`, but
that program is **not on disk** and cannot be re-run from anywhere; its verdict is recorded in the
run's audit findings only. `verify.py` in this folder is a fourth, self-contained implementation
and does not depend on it.

The published control matrices in `verify.py` were transcribed by hand from `lect_dse.tex`. That
transcription is not independently attested inside this folder; a referee comparing them against
the source is comparing against the source, which is the point of printing the `C_H(8,4)` in the
paper.

## Scope

- **Only the upper bound is claimed.** `q_4^HS <= 9` is what the printed matrix proves. Whether
 `q_4^HS = 9` additionally needs non-existence at `q = 7` and `q = 8`. Lemma 4 of the paper
 disposes of `q <= 6` by hand (`q >= 2k - 1 = 7`), but `q = 7` and `q = 8` were decided in the
 run only by **SAT-solver UNSAT** — three solvers, two encodings, and at `q = 8` also under a
 deliberately weakened symmetry break. No unsatisfiability certificate was retained, an
 independent exhaustive census of those two cells was written and **timed out with zero bytes of
 output**, and `verify.py` does not attempt either cell. It says so in its own closing scope
 lines.
- **`q = 10, 11, 12` are not claimed.** The same sweep returned `C_HS(q,4)`-matrices at
 `q = 10, 11, 12`. They are not printed in the paper and not examined by `verify.py`. The
 `q = 10` object in particular was never re-checked by an independent implementation.
- **The `k = 3` improvement is not claimed.** The sweep also reported a satisfying assignment at
 `(q,k) = (6,3)`, which would give `q_3^HS <= 6` against the published `7`. **That matrix is not
 exhibited anywhere in this folder**, so nothing here verifies it and the paper claims it only as
 a reported by-product.
- **Open problem 2 is untouched.** The source's second question (`k >= 5`, `q < k^2`) is not
 addressed.
- **The source's published minima are not re-verified as minima.** `verify.py` re-checks the
 objects the source prints, never the claim that nothing smaller exists.
- **Quotations are from the arXiv v1 e-print source, not the printed lecture notes.** There is
 only one arXiv version (submitted 29 January 2014), so no revision withdrew the question; but
 the POSTECH Lecture Note Series 10 (February 2004) text was **not obtained**, so its numbering
 of the propositions, examples and open problems may differ. The paper reproduces (P5)–(P7) and
 the question in full so that nothing depends on the labelling. One locator correction found in
 review is worth recording: the "Open problems" box closes the subsubsection *Examples of hash,
 separable and hash&separable codes*, not the following *Existence of hash and hash&separable
 codes*, which opens at line 3796.
- **Prior-art search answered but is not exhaustive.** Semantic Scholar, OpenAlex, arXiv and
 OpenCitations queries were run and nothing asserting a `C_HS(q,4)` with `q < 13`, or any value
 of `q_4^HS` other than `13`, was found. The citing literature of the arXiv record was
 enumerated (OpenAlex 51 of 51, Semantic Scholar 55 of 55) and so was that of the journal home
 of the target section, D'yachkov–Rykov, *Optimal superimposed codes and designs for Rényi's
 search model*, JSPI **100** (2001) 281–302 (Semantic Scholar by DOI, 43 of 43); these are two
 citation nodes thirteen years apart and both were swept. Real gaps remain, and they are gaps,
 not clearances:
 - **This row must not be described as zbMATH-checked.** Of six zbMATH queries only **one**
 answered; the other **five returned HTTP 404** under a fetcher return code of 1, which is
 throttling or non-indexing of the phrasing, not zero hits. Three further phrase queries
 run independently at the audit stage returned HTTP 502, 404 and 404. A throttled channel
 is unread, not empty.
 - The full texts of Kim–Lebedev (JCD 2004) and Kim–Lebedev–Oh (JCD 2005) are **unread**
 (Wiley-paywalled; their zbMATH reviewer descriptions and complete reference lists were read
 instead). Füredi–Ruszinkó name these as the small-values literature for exactly this object.
 - Lebedev, *Some tables for (w,r) superimposed codes*, ACCT-8 (2002) 185–189, is **unread**;
 it is a table of small values for this family, and a targeted zbMATH search for it returned
 a different paper.
 - MathSciNet was **not reachable** from this harness.
 - The Semantic Scholar abstracts of D'yachkov's own later *Superimposed Codes and Designs for
 Group Testing Models* (2012) and *Superimposed Codes and Threshold Group Testing* (2014)
 both returned **HTTP 429** and were never read.

 So openness is **supported but not established**, and nothing in this folder should be read as
 a completed openness proof. The nearest published neighbour that was read is Füredi–Ruszinkó's
 grid-free characterisation, a counterpart of Proposition 7 that supplies no small explicit
 object; it is in the paper's reference list.
- **The reduction of Lemma 3 is proved, not assumed.** It is used to make Proposition 7 small, and
 `verify.py` also runs Proposition 7 with no reduction at all, and separately confirms that the
 reduction can return a **non-empty** answer on the negative controls `C_H(6,3)` and `C_H(8,4)`.
