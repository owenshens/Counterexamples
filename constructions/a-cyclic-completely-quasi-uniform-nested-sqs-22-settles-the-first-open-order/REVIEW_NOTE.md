# A Cyclic Completely Quasi-Uniform Nested SQS(22)

`a-cyclic-completely-quasi-uniform-nested-sqs-22`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |
| `REVIEW_NOTE.md` | this note |

Those five files are the whole folder; nothing here refers to a file it does not ship, and every
external file named under *Provenance* below is flagged there as not shipped.

## What is claimed

One existence statement: **a completely quasi-uniform nested SQS(22) exists**, exhibited by an
explicit object printed in full in Section 2 of the paper as 19 base nested blocks over `Z_22`
together with the rule `x -> x+1`. No side results are proved: neither the design nor its
splitting is claimed to be new, and no claim is made that `Z_22` is the full automorphism group.

The target is a **statement of the state of knowledge**, not a numbered conjecture. Lu writes, in
the last paragraph of his section "More results for small non-Boolean orders" (arXiv:2509.06663v3,
e-print source `UnifNestBoolSQS_rev1109.tex`, bytes 64,824-65,106, lines 1037-1039), that for
`v = 4 (mod 6)` "apart from the Boolean cases and the smallest nontrivial case SQS(10) ... no
further examples are currently known". No question or conjecture is attached to that remark, so
**nothing here refutes anything in print and nothing posed is settled**. The smallest order not
covered by the two results Lu cites is 22 (10 is Lu's own example, 16 = 2^4 is Boolean), and 22
is the order treated.

Three things a referee should hold onto:

- **This is one order, not the class.** 28, 34, 40, 46, ... remain open. Problem 3 of Chee, Dau,
 Etzion, Kiah and Zhang, which asks for an *infinite* family of quasi-uniform nested SQSs, is
 untouched -- and note that the known set was **already infinite** through Lu's Boolean orders
 `2^m`, `m` even, so one more order does not bear on it. Two extension searches, at `v = 28` and
 `v = 34`, were run and returned INCOMPLETE: **zero** shards finished inside their time budgets, so
 they are not evidence in either direction about those orders.
- **No count and no isomorphism test.** Nothing bounds the number of such systems at order 22, and
 no isomorphism test was run -- not against the second (`Z_11`-invariant) example found during the
 search, and not against any published cyclic `S(3,4,22)`. No catalogue of cyclic quadruple
 systems was searched. The claim is that the printed object *is*
 a completely quasi-uniform nested SQS(22), not that it is new *as a design*.
- **Definitions.** "Quasi-uniform" is used in Lu's sense, which *excludes* uniformity: the
 multiplicities must not all be equal. Chee et al. drop that clause. The object here has
 multiplicities 3 and 4, so it qualifies under either reading; the class sizes `154, 77` are
 counted for this object, and nothing here shows the pair `{3,4}` is forced.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **22 checks, all passing**:

 VERDICT: ALL 22 CHECKS PASS

It reads the object exhibited in the paper as input -- the 19 lines of the Section 2 table, pasted
verbatim into the program and re-parsed from that text -- and derives every quantity it compares
against the paper's statements. All arithmetic is over `int`; no floating-point value is computed,
so no decision depends on one.

What the 22 checks cover, by step:

1. the printed table parses to 19 nested blocks with labels `B01..B19`, every entry in `0..21`,
 each a 4-set split into two disjoint pairs, and the printed orbit lengths summing to 385;
2. **Lemma 2**: the computed orbit length under `x -> x+1` equals the printed one for all 19 rows;
 every orbit has length 22 or 11; the three short orbits are exactly `B14, B17, B19`; a base block
 has a short orbit **iff** its quadruple is a union of two diameters `{a,a+11}` (the biconditional,
 both directions, on all 19 rows, plus the general fact that all three splittings of a
 diameter-union quadruple are fixed by `x -> x+11`); the printed diameter decompositions are
 correct; the 19 orbits are pairwise disjoint and give 385 nested blocks; the set is
 `Z_22`-invariant, verified element by element;
3. **Theorem 1, the design**: 385 distinct underlying quadruples, `C(22,3) = 1540`, and every one of
 the 1540 triples covered **exactly once**;
4. **Theorem 1, the nesting**: `C(22,2) = 231`, all 231 pairs occur (complete), the profile is
 exactly `{3: 154, 4: 77}`, two values differing by 1 and not all equal (quasi-uniform, not
 uniform), and multiplicity sum `770 = 2*385`.

Every value the paper prints is recomputed and compared, so each comparison can actually fail.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 d600f1c048b9fc09ee4f9dedd6bcbcf0f0be04b6af12a1935e3a77b10d5b57ca

`verify.py` was written for this folder and run here, on the machine that assembled it, under
Python 3.9.6; it exited 0.

**Provenance of the search that produced the object** (recorded because it is *not* what `verify.py`
re-runs). The files named below are **not shipped in this folder** -- see the *Contents* table for
what is. They live in the underlying run record, and the digests are
pointers for matching them there, not files a referee can hash here. Nothing in this folder depends
on them: the object is printed in full in Section 2 and `verify.py` checks the printed object.

Per the artifacts manifest of the underlying record:

- The 19 base nested blocks of Section 2 were produced by the `Z_22`-cyclic sub-cell of
 `job2.py` (SHA-256 `88865de3d950eb567f290f24c40ffbc22fe1fc0f7f21bc30737a06887ac17132`,
 24,509 bytes), invoked as

 aws/slot_run.sh S02 runs/wave24/artifacts/t10841/job2.py 2400 a2t10841verify

 on one 32-vCPU worker (), status Success, exit 0, 17,716 bytes of captured stdout.
 **No solver of any kind**: stdlib `itertools`, `collections`, `math.comb` and `multiprocessing`
 only -- no SAT, no ILP, no nauty. The manifest records **no randomness anywhere in that file** and
 a shard order fixed by the depth-0 candidate list, so that search is deterministic and would
 return this same object again. Its lane was **not exhausted**: 15 of 24 shards finished and 9 hit
 a 400 s budget, which is why nothing is claimed about the *number* of cyclic examples.
- A second, `Z_11`-invariant object, not mentioned in the paper as it now stands, came from a different script,
 `nsqs22.py` (SHA-256 `d92f39eba120819e9277d545db096ef5d482cb9ebeeec99044f7d41a5607aac5`,
 13,334 bytes), on, status Success, exit 0, 5,705 bytes of stdout. It is not
 printed in the paper and is not checked by `verify.py`.
- The machine-readable certificate in the run record, `witness.z22.txt` (SHA-256
 `700187a1979d7215f3bd368df433dbe39d3d1b247ad7adb7d3d7a89ddd377dff`, 4,075 bytes, 385 nested blocks
 one per line), was **derived** from the captured stdout by `witness_extract.py` (SHA-256
 `bb4bec2a4db6c27751291dd4ec86cff8500f8688d2869cf3fdf79e2c083264b2`, 7,986 bytes) rather than
 retyped from prose, and cross-checked against a JSON encoding of the same object: symmetric
 difference 0 over 385 labelled nested blocks. The table printed in Section 2 is the 19 base blocks
 of that certificate.

Four gaps in that record, stated rather than dressed:

1. **No elapsed time was recorded for either slot job.** `slot_run.sh` printed no wall-clock figure
 into the captured stdout and none exists on disk, so none is given here. Neither command timeout
 (2400 s, 900 s) was hit.
2. **No S3 copy of either slot output exists** (`head_object` 404 on all planned keys), so the
 captured `.out` files in the run record are the only copies of that stdout.
3. The `job2.py` dispatch was filed in the record with its SSM command id in two non-identical
 forms, `...f1b3...` and `...91b3...`. `slot_run.sh`'s own dispatch log reads `f1b3`, which is what
 is quoted above; `91b3` is a transcription error. Nothing in the result depends on it. That
 dispatch also ran as a post-decision verify over an already-terminal row and is logged in the
 run's SSM override log rather than its main dispatch ledger.
4. A third script in the record, a by-hand transcription re-check of the `Z_11` object, has **no
 captured slot stdout at all** -- that dispatch left no log -- so its reported slot runtime is a
 claim in the record and not a reading. It bears on the object *not* printed here.

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE: this program checks the single completely quasi-uniform nested SQS(22) printed in the note: its orbit lengths, 385 distinct underlying quadruples, unique coverage of all 1540 triples, completeness, multiplicity profile {3: 154, 4: 77}, and invariance under x -> x+1. It does not test isomorphism, enumerate other systems, or determine the full automorphism group.

Beyond that:

- **The locator has since been checked by hand, outside `verify.py`.** Item (d) of the program's
 scope note -- the byte locator, which the program does not fetch -- was checked at bundle review by
 downloading `https://arxiv.org/src/2509.06663v3`. The source file is
 `UnifNestBoolSQS_rev1109.tex` at **68,869 bytes**, matching the paper; 0-indexed bytes
 `[64824, 65106)` are **282 bytes** and are exactly the quoted sentence, character for character
 including its LaTeX macros; those bytes are lines **1037-1039**; and that paragraph is the last of
 the section `More results for small non-Boolean orders` (line 788, running to 1041, with
 `\section{Concluding remarks}` at 1042). Every `\label` the paper cites resolves in that file:
 `thm:all_quasi_uniform` (313), `prop:split_2designs_4` (420), `thm:main_boolean` (567),
 `ex:16` (269), `ex:cqu-SQS10` (824), `tab:cu_SQS_le_50` (1016), and the definition lines 108-116,
 119, 121, 124, 131, 132. In particular Lu's `ex:16` is confirmed to be a union of full `Z_15`
 orbits **together with four partial orbits over `t in {0,1,2,3,4}`**, so it is not
 `Z_15`-invariant. This
 check is a bundle-review addition and is **not** part of the 22 checks; `verify.py` was not
 changed and still does not fetch those bytes.
- **The statements attributed to Lu and to Chee et al. are transcribed, not re-read by
 `verify.py`.** The paper says so explicitly. Nothing load-bearing rests on those transcriptions:
 the paper uses no result of Lu or Chee et al. as input, and every quantity it states about the
 printed object is recomputed by the program. Lu's result
 labels are quoted by their e-print `\label` names, not by journal numbers, because the e-print
 source is the text that was read; the paper says this explicitly.
- **Bibliographic limits on "open".** Novelty is claimed only against the sentence quoted above.
 Searches of the arXiv metadata API (five queries), zbMATH's document search, Crossref from
 2025-11, and the single-record citation counts of both source papers at Semantic Scholar and
 OpenAlex all returned nothing at order 22, and both e-prints were fetched and grepped at source
 with zero occurrences of order 22. But: arXiv's `all:` field is metadata, not full text; zbMATH's
 editorial text for the older cyclic-SQS literature was licence-blocked; **OpenAlex's works-search
 and `cites:` endpoints, and Semantic Scholar's keyword-search endpoint, never answered** (HTTP
 429), so a citer indexed by one and not the other would have been invisible; and **MathSciNet was
 not consulted at all**, which for a 2026 journal paper is a plausible home for a prior
 small-order remark. No correspondence or erratum channel was checked. Nothing we did excludes
 that such an object has been recorded somewhere we did not look.
- **Two published results were checked for collision and do not collide.** Chee et al.'s own `v = 22`
 table row bars a *uniform* complete nesting at that order and leaves the remaining cell marked as
 having no known design, which corroborates rather than pre-empts this construction; and their
 doubling lemma runs the other way, taking an SQS(22) as *input*. Lu's own remark that a completely
 quasi-uniform nested SQS(10) cannot be obtained by cyclic translation of base blocks is about one
 order and the cyclic group `Z_10`, and says nothing about order 22.
