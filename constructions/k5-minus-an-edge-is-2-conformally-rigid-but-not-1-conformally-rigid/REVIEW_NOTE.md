# K_5 Minus an Edge is 2-Conformally Rigid but not 1-Conformally Rigid

`k5-minus-an-edge-is-2-conformally-rigid-but-not-1-conformally-rigid`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

Assumpção, Coutinho and Godsil (arXiv:2605.08508v2) define *k-conformal rigidity* for
k in {1,...,n-1} and ask, in the final paragraph of their Conclusion, whether some graph is
2-conformally rigid but not 1-conformally rigid. The paper answers **yes** with `H = K_5 - e`
on five vertices, and both halves are certified by objects printed in full in the paper:

- **the answering half.** The paper makes no claim of priority for the upper-2 fact itself; whether it already follows from results or data in the two sources was not checked. The integer matrix `M` of equation (3) in Proposition 2 of the
 paper, with `X = M/45`, is a
 dual certificate: `tr X = 2`, `X.1 = 0`, `0 <= X <= I` exactly (spectrum
 `{0, 0, 5/9, 5/9, 8/9}`), and `L*(X) = 10/9` on **every** one of the nine edges, with edge
 sum `9 * 10/9 = 10 = S_2(1)`. By the Ky Fan maximum principle this gives
 `S_2(w) >= 10 = S_2(all-ones)` for every `w` in `Delta_E`, so `H` is upper-2 rigid; at
 `n = 5` and `k = 2` the source's own duality has `n-1-k = 2 = k`, so lower-2 rigidity is
 the same statement and `H` is 2-conformally rigid.
- **the negative half (already published as a fact).** The weight `w*` equal to `9/7` on the
 six edges meeting `{4,5}` and `3/7` on the three edges inside `{1,2,3}` lies in `Delta_E`
 and has `spec L(w*) = {0, 27/7, 27/7, 27/7, 45/7}`, so `s_1(w*) = 27/7 > 3 = s_1(all-ones)`
 and all-ones does not maximise `s_1`. That `K_5 - e` is not lower-1 rigid, hence not 1-conformally rigid, is
 **published numerical data**: row `id=450` of `hog/numerical/numerical_certs.csv` in Niu's
 dataset (arXiv:2605.15017v1, `github.com/andrewmniu/conformally-rigid-graphs`) reads
 `lcr=False, ucr=True`. The paper says so in Section 1. What is ours in this
 half is only the exact rational `w*`, which turns that table cell into a hand-checkable
 witness.

The paper's Section 1 also flags a definitional trap that a referee will hit:
Niu defines "conformally rigid" as lower-1 **or** upper-1, while Steinerberger–Thomas and
Assumpção–Coutinho–Godsil use **and**. Under the disjunctive reading `K_5 - e` *is*
conformally rigid.

## What was checked, and how

Every step of both proofs is 5x5 rational arithmetic and can be followed by hand from the
printed matrix `M` and the printed weight `w*`; the spectra are obtained from explicit
eigenvectors, not from a numerical eigensolver. `verify.py` re-derives all of it
mechanically, in exact arithmetic, from the objects printed in the paper.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only (`fractions`, `itertools`, `sys`): no third-party
package and no external data file. Every decision is an exact integer or `Fraction`
comparison; no floating-point value is ever compared. The program prints one line per check
and a closing verdict, and exits 0 only if every check passes. The recorded run reports
**68 checks, all passing**:

 VERDICT: ALL 68 CHECKS PASS

Its inputs are exactly the objects exhibited in the paper -- the graph6 string `D^{`, the
printed relabelling onto `{1,...,5}`, the printed nine-edge list, the integer matrix `M`, and
the weight `w*` -- and it re-derives, among other things: the decoding of the graph6 string
and its agreement with the printed edge list under the printed bijection; the characteristic
polynomial and exact spectrum `{0,3,5,5,5}` of the Laplacian and the values
`S_1, S_2, S_3, S_4, s_1` at all-ones; each of the four properties of `X` used in the proof,
including its characteristic polynomial, the exact containment `0 <= X <= I`, `L X = 5X`, the
constancy of `L*(X)` and its edge sum; the pairing identity
`tr(L(w) X) = <w, L*(X)>` on several exact rational weight vectors; the closed form
`X = (5/9) P + c c^T / 90` and that `P` is the projector onto the 5-eigenspace; that `w*` lies
in `Delta_E` and has the stated exact spectrum and eigenvectors; the two strict inequalities
`27/7 > 3` and `99/7 < 15`; the duality identity `S_k + s_{n-1-k} = 2|E|` for both weights and
all admissible `k`; the redundant upper-1 certificate `X/2`; the values of `L*(P)`; and the
resulting k-map at `n = 5`, namely that `H` is k-conformally rigid exactly for `k` in
`{2, 4}` -- these last three being recomputations the program supplies beyond the paper, which
states no k-map, no closed form for `X` and no upper-1 certificate.

## Provenance

**Of this folder's program.** `verify.py` was written for this folder and run locally on this
machine; it is not one of the scripts that produced the result. `verify.output.txt` holds its
output, preceded by a provenance header and followed by an exit status, both written by the
run harness. The header records the SHA-256 of the program that produced the output, so the
two files can be matched:

```sh
shasum -a 256 verify.py
```

 55722671189622c240a483b73241ad753f8d028d3a0b365ffbaf84efcc61853b

The recorded exit status is 0. `verify.py` reads nothing outside itself, so a referee needs
neither this folder's other files nor any of the material described next.

**Of the original computation.** `M` and `w*` were produced by an internal run whose
artifacts are recorded with per-file SHA-256 digests in that run's manifest. What that
manifest actually records, stated as it records it and not improved upon:

- the decisive script, which verified both halves in exact rational arithmetic (`sympy`
 `Rational`, no floating point), is archived with its digest, its invocation, and the
 identifiers of the machine it ran on; its complete standard output is archived beside it,
 and the numbers in that output are the numbers printed in this paper.
- the two halves were re-verified independently by a second and third script that share no
 code with the first; those are archived with their invocations and their outputs.
- **honest gaps, recorded in the manifest rather than papered over.** Four further scripts
 in that run -- including the two that swept the census of all 994 connected graphs on 3 to
 7 vertices -- had their standard output excerpted into an internal note and **never written
 to a file**, and the credentials needed to retrieve the archived copies had expired before
 the manifest was written, so those outputs are *unread, not absent*. One of those four also
 has **no command identifier captured**, so its run is not independently traceable. Elapsed
 runtimes were recorded for none of the seven scripts. **None of this bears on the paper**:
 the census is not a claim of this paper (see Scope), and every input the decisive script
 consumed is printed in the paper.
- the run used no randomness anywhere: the certificate arm is exact rational arithmetic, and
 the numerical arm that suggested the witness was a deterministic cutting-plane loop with a
 fixed iteration cap and no seed.

No reproduction instruction beyond `python3 verify.py` is offered here. In particular this
note does not tell a referee how to re-run the census, because the record does not support
such an instruction: its output was not preserved.

## Scope

What the paper claims is that this one graph is 2-conformally rigid and is not
1-conformally rigid, and that claim is certified in full. What is **not** settled:

- **no minimality.** The paper does not claim that `K_5 - e` is the smallest graph with this
 behaviour, nor that it is the only such graph on five vertices. Those questions were
 examined only by a float64 census at a bound tolerance near `1e-7`, whose output is not in
 the record (above); orders `n >= 8` were never examined at all. Nothing in the paper rests
 on any of it, and `verify.py` certifies none of it.
- **the other two open questions of the same paragraph of the source remain open.** In
 particular the question of whether, for every constant `C`, some graph is k-conformally
 rigid for all `k <= C` but not for `C+1` is untouched: our witness has the reverse pattern.
- **`n = 5` is degenerate in the result's favour,** and the paper says so in a remark: it is
 the smallest order at which a witness can exist, and at `k = 2` the two halves of
 2-conformal rigidity collapse into one condition. The witness therefore says nothing about
 how the upper and lower conditions interact where they genuinely differ.
- **two quoted theorems are not reproved** by `verify.py`: the Ky Fan maximum principle, and
 the upper-k / lower-(n-1-k) duality of the source (which the paper does prove, in two
 lines, as Lemma 3).

The program's own closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE -- what this program does NOT cover. It re-derives, in exact arithmetic, the quantities named in the PASS lines above for the single graph K_5 - e. Its inputs -- the graph6 string, the edge list, the matrix M and the weight w* -- are TRANSCRIBED BY HAND from the paper and hard-coded here: this program never reads paper.tex, so neither the faithfulness of that transcription nor the completeness of the list above against what the paper states is machine-checked, and "every quantity the paper states" is not something this program verifies. NOT RE-RUN here: (a) the Ky Fan maximum principle and the source's upper-k/lower-(n-1-k) duality, which are quoted theorems, not recomputations; (b) minimality and uniqueness -- the paper claims NEITHER, asserting neither that K_5 - e is the smallest graph with this behaviour nor that it is the only 5-vertex witness; the census of the 994 connected graphs on 3..7 vertices that looked at those questions was float64 at a bound tolerance near 1e-7 and is NOT reproduced here, so nothing here bears on them either; (c) every order n >= 8, which was never swept; (d) any infinite family. BEYOND THE PAPER, in the other direction: all of Step 5 (the upper-1 certificate X/2, the values of L*(P), the k-map) and the closed form for X in Step 3 are established here but are not stated or claimed by the paper. This program certifies the two witnesses and nothing else.

## One thing a referee should check against the source

The question answered here is unnumbered running prose, not a numbered Question or Problem:
it is one of the questions in the final paragraph of Section 9 (Conclusion,
label `sec:conclusion`) of `main.tex` in arXiv:2605.08508v2, at line 1072, bytes 63519-63597
of a 64127-byte file. It must not be cited as "Question 9.1" nor as "the first open
question", which names a different sentence. The paper reproduces the sentence verbatim so
that what is answered does not depend on how it is labelled.
