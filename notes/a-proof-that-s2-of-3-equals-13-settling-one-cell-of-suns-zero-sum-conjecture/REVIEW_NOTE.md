# A Proof that s_2(3) = 13, Settling One Cell of Sun's Conjecture 1.2

`a-proof-that-s2-of-3-equals-13-settling-one-cell-of-suns-zero-sum-conjecture`

Supporting material for this paper: the program that re-derives its finite claims from the two
objects printed in it, and a record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

For integers `n > 1` and `r >= 1`, Zhi-Wei Sun's `s_r(n)` (the environment labelled `Def1.2`,
part (i), of the e-print source of arXiv:2606.18234v5) is the least `k` such that any `k` vectors
in `Z^r`, none congruent to `0` mod `n`, admit an index set `I` with **`|I| = n` exactly** and
`sum_{i in I} a_i = 0 (mod n)` but `!= 0 (mod n^2)`. Conjecture 1.2 of that paper asserts
`s_1(n) = 2n+1` and `s_2(n) = 4n+1` for all `n > 2`. The paper proves the single cell

 s_2(3) = 13,

which is the `n = 3` instance of the `s_2` half.

**Attribution, which is most of what a referee needs to weigh this.**

* The **lower bound** `s_2(3) >= 13` is **Sun's**, the left half of his published bracket
 `3*2^r + 1 <= s_r(3) <= 2*3^r - 1`. The size-12 multiset `W` printed in Section 5 is the unit
 multiple `4` of Sun's own one-dimensional extremal set placed on each coordinate axis. It is
 reproduced only so the paper is self-contained; it is a control, not a contribution.
* **Fact 3** of the paper is **Sun's Lemma 2.1** and **Fact 4** is **Sun's own collapsing claim**,
 both quoted with source line numbers. Together with the elementary Fact 2 they give the master
 bound `|S| <= 3|T| - |A(T)|`, which is also Sun's.
* The **only new content** is Section 4: Sun double-counts the master bound to `2(3^r - 1) = 16`
 and stops, which is how his published upper bound reaches `17`; a case analysis on the mod-3
 support size `|T| in {0,...,8}` carries `16` down to `12`. Seven of the nine cases are
 immediate; `|T| = 7` and `|T| = 8` need a linear elimination mod 9.
* **Sun already asserts the value.** Remark 1.5 of the same paper (lines 287-288) says *"It is
 easy to see that `s_r(3) = 2 x 3^r + 1` for `r = 1, 2`"*, gives no proof, and is followed seven
 lines later by the same value being posed as an open conjecture. As printed the formula yields
 `19`, which that paper's own bracket refutes; the intended formula is `3 x 2^r + 1 = 13`. So the
 honest description of this result is **the first proof of an author-asserted value**, together
 with a correction of a printed misprint that is not offered as a separate result. The argument is
 short, and the paper says so rather than dressing it up.

## What is *not* settled

* **Conjecture 1.2 remains open.** One cell of one half is proved.
* **`n >= 4` is untouched and the method does not generalise in `n`**: it uses the coincidence that
 the required index-set size `n = 3` equals the number of points on a line of `AG(2,3)`, so
 Facts 2 and 4 have no analogue for `n != 3`.
* **`r >= 3` is untouched**, and the source itself records that the pattern breaks there
 (`s_3(3) > 25`).
* **The `s_1` half is not Sun's conjecture** but Gao-Jiang-Lei-Lin-Yang's, and its `n = 3` instance
 `s_1(3) = 7` is already pinched between the two sides of Sun's own bracket at `r = 1`. It is a
 control here, not a result.
* **No minimality or uniqueness.** The bound 12 is attained at `|T| = 4`, `5` and `6`, so extremal
 multisets are not unique even up to the `GL_2` action, and the paper claims nothing about them.
* **Nothing about `t_r(n)`**, the source's sibling family of invariants.
* **Statement numbering.** Our record of the *printed* numbering of Sun's definition is internally
 inconsistent (it appears as both "Definition 1.2" and "Definition 1.3"), so the paper cites the
 LaTeX labels `Def1.2`, `Th-value`, `Th-bound`, `Conj1.2` with line and byte offsets in
 `P95z.tex` instead. "Conjecture 1.2", "Lemma 2.1" and "Remark 1.5" are printed numbers as read.
 Both statements the paper touches are reproduced in full, so what is proved does not depend on
 the labelling.
* **Priority search, and its two holes.** The prior-art pass ran ten channels; the nearest source
 found anywhere is Sun's own Remark 1.5, described above. Two channels did not answer and a
 referee should know: zbMATH's review text and reference lists for the target and for three Gao
 et al. papers are **licence-blocked** (the API returns the literal "contents unavailable due to
 conflicting licenses"), and **MathSciNet was not attempted** (no institutional access), which for
 a result whose only near miss is an Acta Arithmetica paper is the likeliest remaining place for a
 prior observation. The full texts of the two paywalled Gao et al. papers were read as abstracts
 only. Semantic Scholar reports exactly one citer of the source; that citer's full LaTeX was read
 and contains no occurrence of `s_r`, `s_2` or `4n+1`.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package, no external data file, no
network. Its inputs are the two objects **printed in the paper** and nothing else: the multiset `W`
of equation (6) and Table 1 of Section 3 (the four antipodal pairs, the eight origin-missing lines
of `AG(2,3)`, and the eight equations `E1..E8` in the paper's own naming `A, A', B, B', C, C', D,
D'`). Both are transcribed into the program as text blocks and parsed, so a referee can compare the
program's input against the printed page character by character. All arithmetic is on integers mod 3
and mod 9; there is no floating point in the program and no decision depends on one. It prints one
line per check and a closing verdict, and exits 0 only if every check passes. The recorded run
reports **64 checks, all passing**:

 VERDICT: ALL 64 CHECKS PASS

What the 64 checks cover, in the program's seven steps:

1. **The ambient combinatorics and the transcription.** 72 admissible classes in `(Z/9)^2`; the
 reduction-triples of Fact 2 split as 8 all-equal and 8 all-distinct with **zero mixed**; the 8
 all-distinct ones are exactly the origin-missing lines; each of the 8 nonzero points lies on 3
 of them; no such line holds an antipodal pair and each meets exactly 3 of the 4 pairs; the
 printed pairs `P1..P4` really are the antipodal pairs; **the eight rows of Table 1 are, as a set
 of triples, exactly the eight computed lines**, each row's class symbols name the points printed
 beside them, and each symbol is used three times; and `3a != 0 (mod 9)` for all 72 admissible
 `a`, which is the mechanism of Fact 3(a). The Table 1 checks matter: the whole upper bound runs
 through `E1..E8`, and a single mis-transcribed row would change the eliminations.
2. **The printed witness `W`.** Admissible, BAD, size 12, class multiplicity 2, support exactly
 `P1 u P2` of size 4 with all four fibres of size 3, containing no full line so `|A| = 0`, and
 with exactly four candidates whose sums are `(9,0), (18,0), (0,9), (0,18)` -- all `0` mod 9.
 Plus **maximality**: all 72 admissible one-element extensions fail to be BAD.
3. **Directed negatives, so a "not BAD" answer means something.** Three copies of `(1,0)` are not
 BAD (the must-fire polarity of Fact 3(a)); and `W` is **not** BAD under the `|I| <= n`
 misreading of the definition, which is why the paper stresses that `|I| = n` exactly is
 load-bearing.
4. **The master bound over all 2^8 = 256 supports.** `max(3|T| - |A(T)|)` is
 `0, 3, 6, 9, 12, 12, 12, 14, 16` for `|T| = 0..8`; it settles 247 supports and leaves exactly
 the `C(8,7) + C(8,8) = 9` with `|T| >= 7`; and `|T| = 5, 6` are tight at 12.
5. **The three eliminations of Section 4, as identities over all `9^3 = 729` lifts of `(A,B,C)`**:
 `B + C + D' = 3C`, `B' + C' + D = -3(A+C)`, and `A' + C' + D' = -3(B+C)`, each mod 9, with every
 derived class checked to have the right reduction mod 3, plus the three reductions
 `pi(C) = (1,1)`, `pi(A)+pi(C) = (2,1)`, `pi(B)+pi(C) = (1,2)` that make them contradictions.
6. **Stronger than the paper argues.** The paper relabels by `GL_2` and treats one representative
 support of size 7. The program instead solves **every one of the nine remaining supports
 separately** -- all eight with `|T| = 7` and the one with `|T| = 8`. For each it verifies that
 the full lines cover `T` (so every fibre is frozen), then substitutes `x_v = v + 3 u_v` with
 `u_v in F_3^2`, which turns each full-line equation `x_p + x_q + x_r = 0 (mod 9)` into a linear
 equation `u_p + u_q + u_r = -w_L` over `F_3`, and solves those systems exhaustively: all nine
 are inconsistent. This does not rely on the relabelling used in Proposition 5.
7. **The source's own arithmetic**: `4n+1 = 13` at `n = 3`; the bracket `13..17` at `r = 2`; that
 Remark 1.5 as printed gives `19 > 17`; that the two candidate formulas agree at `r = 1` (both
 `7`) so `r = 2` is the only place they can be told apart; and that the double count gives 16.

The program runs in about 0.1 s.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an exit
status, both written by the run harness. The header records the SHA-256 of the program that produced
the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 4c48d631769d57b3339db7838ca0625a907f063330e5ca3b6331c34579544843

The run recorded there was made on the author's laptop with Python 3.9.6 and exited 0.

**`verify.py` is fresh code written for this folder, and it is not the program that produced any
object here.** `W` is Sun's construction, read off his own proof. The provenance of the machine work
that *echoed* this result, as recorded in the result's artifacts manifest
(`runs/wave23/artifacts/t5904/MANIFEST.json`, 7 files with a SHA-256 each; that manifest and those
files are **not** shipped inside this folder and nothing here depends on them), is:

* `fast.py` (sha256 `4eb3b5a4a50efa368b116031e0869afe653c66223bf627c71955b19f5bfbf2fd`, 9,453 B) ran
 an exact CP-SAT maximisation over the 72 classes with multiplicity capped at 2 and 6,288 enforced
 triple constraints, using none of the fibre structure of the paper. It returned
 `status=OPTIMAL alpha=12 bound=12` in 67.4 s on slot **S31**, instance
 **a separate cloud instance**, SSM CommandId **a recorded command id**,
 `STATUS Success`, `RC = 0`, `ortools 9.15.6755`, 32 vCPU. Its stdout is filed as `fast.out`
 (5,172 B) and the manifest marks it complete: the wrapper's own tail line reads
 `ARTIFACT_STDOUT_BYTES=4141`, well below the ~24 kB SSM ceiling, and the job's final `DONE` is
 present. It also reproduced nine published integers of the source to proved optimality as
 controls, and its decision-form run (does a BAD multiset of size 13 exist?) returned `UNKNOWN` at
 its 900 s cap -- **INCOMPLETE, not a negative**, and superseded by the maximisation.
* `verify_witness.py` (sha256 `28a6648620316fbadd9d6c2a35f6381ae81d8afb251b5dfdab185e1a27c4686f`)
 ran on the control plane in 0.09 s, stdlib only, and re-verified four extremal size-12 multisets
 against the definition; its output is filed as `verify_witness.out`. The manifest's `reproduce`
 field is that invocation, `python3 runs/wave23/artifacts/t5904/verify_witness.py`.

**Gaps in that record, stated rather than filled.** The manifest itself flags them and this note
does not improve on it.

* The two slot invocations in the manifest are marked **COMPOSED, NOT COPIED**: `slot_run.sh`
 printed `ARTIFACT_NOT_FILED` (the jobs were dispatched from `/tmp/aa1-t5904/` without
 `--artifact-wave`), so no `ARTIFACT_INVOCATION=` line exists to quote. The dispatched files are
 recorded as byte-identical to the indexed ones, differing only in path.
* **The stdout of `check_s2_3.py` was never written to a file.** That job (slot **S21**, instance
 **a separate cloud instance**, CommandId **a recorded command id**) returned
 `Success`/`RC = 0` and was an independently-encoded replication, but no S3 artifact was filed and
 its numbers survive only as prose. They are therefore **not** cited in the paper.
* Three further dispatches (`census.py` twice, `mini.py` once) ended `TimedOut` with 0 B stdout:
 **INCOMPLETE, never negatives**, and nothing in the chain rests on them.
* Instance types, slot Python versions and total wall clock are recorded as **NOT RECORDED** in the
 manifest, not estimated. `fast.out` has no S3 copy (`head_object` 404 for all three keys), so the
 local capture is the only one.
* One `UNCHECKABLE` claim in the artifact gate: `witness.json` is JSON rather than a graph, so the
 library's graph6/edge-list cross-decoder cannot read it. The label-equality check was run by
 `verify_witness.py` instead (`MULTISET-EQUAL = True` against the multiset printed in `fast.out`).

**None of that touches this folder.** The paper claims a hand proof plus one exhibited multiset, and
`verify.py` re-derives every finite claim about both from the printed page, with no dependence on
any slot output, on `witness.json`, or on any file outside this directory. The CP-SAT run is an echo
and is not shipped, because it needs a third-party solver.

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE. This program re-derives every quantity the paper claims about the objects printed in
> it, and it settles the two cases the master bound leaves open for EVERY support of size 7 or 8,
> not only for the relabelled representative the paper argues. NOT COVERED, and stated as such in
> the paper: (i) Sun's Lemma 2.1 (Fact 2) and Sun's collapsing claim (Fact 3) are quoted from the
> source and proved there and in Section 3 by hand -- only Fact 2(a)'s must-fire mechanism is
> machine checked here, the cancellation arguments of Fact 2(b) and Fact 3 are not; (ii) the master
> bound |S| <= 3|T| - |A| itself is a counting step, taken from Fact 2(b) and Fact 3 rather than
> re-derived; (iii) NO independent exhaustive search over the 2.25 x 10^34 multisets in {0,1,2}^72
> is performed here -- the certified CP-SAT maximisation that echoes this result (OPTIMAL
> alpha = 12, best bound 12) needs ortools and is NOT shipped in this folder, so the upper bound
> here rests on the hand proof plus the checks above; (iv) no minimality or uniqueness of W -- the
> paper claims none, and other extremal multisets exist at |T| = 5 and 6; (v) n >= 4 and r >= 3,
> which are untouched; (vi) the line and byte locators quoted from the e-print source, which are
> transcribed and not fetched by this program.

Two further limits of the checking, outside the program. First, the numbering in that quoted note
counts the paper's facts in the order they are stated ("Fact 2" there is the paper's Fact 3, Sun's
Lemma 2.1, and so on); the *content* referred to is unambiguous. Second, the byte and line locators
in Section 1 were read off a stored fetch of the v5 e-print tarball (14,100 B, `sha256` beginning
`9686627e`; one file `P95z.tex`, 48,435 B, `sha256` beginning `f4973376`) and are transcribed here.
A reader without that fetch must retrieve the e-print source to confirm them, and a different
arXiv version may move them.
