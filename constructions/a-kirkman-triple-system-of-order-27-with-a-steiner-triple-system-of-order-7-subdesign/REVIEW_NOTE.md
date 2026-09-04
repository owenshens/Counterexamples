# A Kirkman Triple System of Order 27 with a Steiner Triple System of Order 7 as a Subdesign

`a-kirkman-triple-system-of-order-27-with-a-steiner-triple-system-of-order-7-subdesign`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |
| `REVIEW_NOTE.md` | this note |

## What is claimed

One existence statement: **a KTS(27) with an STS(7) subdesign exists**, exhibited by the thirteen
parallel classes printed in full in Table 1 of the paper. Both halves of Theorem 1 are inspections
of that table -- each class is a partition of `{0,...,26}` into nine triples, the 117 triples cover
each of the 351 pairs exactly once, and the seven triples lying inside `U = {0,...,6}` are a Fano
plane -- so **the decisive argument needs no program**. A referee who does not want to run anything
can check the result from the printed table alone.

Two further statements are proved and are not the main claim: Lemma 2 (a proper subdesign of order
`w >= 3` in an STS(v) forces `v >= 2w+1`, elementary and standard) and Proposition 3 (the design has
exactly one subdesign of order 7 and none of order 9 or 13, so the number of Fano subdesigns through
a point is 1 on seven points and 0 on twenty, and the automorphism group is therefore not transitive
on points). Proposition 3 is proved by closing all 2925 point triples of the printed design; it appeals to nothing outside Table 1.

Three things a referee should hold onto:

- **The antecedent is a hedge, not a declared open problem, and the paper says so.** Dukes and
 Lamken write that a KTS(27) with an STS(7) subdesign "appears difficult to find", and reserve the
 word "unknown" for `v = 33` and `v = 39` in the very same sentence. So what is claimed is only the existence of a
 design for the cell `(v,u) = (27,7)`, exhibited here, **not** priority or novelty for it and **not**
 that Stinson's Open Problem 4 is settled; whether such a design has been exhibited before was not
 checked.
 The orders 33 and 39 are untouched.
- **The subdesign is not required to be resolution-compatible, and that is the source's own
 convention, quoted verbatim in Section 1.** Dukes and Lamken distinguish a *subdesign* from a
 *Kirkman subsystem* and say explicitly that no extra structure is assumed for subdesigns. At
 `u = 7` the stronger notion is vacuous anyway, since `3` does not divide `7`.
- **One object, no count, no group order.** No enumeration and no isomorph rejection was performed,
 so nothing here bounds how many such designs exist, and no minimality or canonicity is claimed.
 The order of `Aut(D)` is **undetermined**: all that is claimed is that it contains the order-3
 automorphism `sigma` of Section 4 and, by Proposition 3, that it is not transitive on points. An
 earlier lens in the underlying record printed `|Aut| = 1`, which contradicts its own verified fact
 that `sigma` is an automorphism; that figure is withdrawn and appears nowhere here.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **59 checks, all passing**:

 VERDICT: ALL 59 CHECKS PASS

It reads the object exhibited in the paper as input -- the 26 lines of Table 1, pasted verbatim into
the program and re-parsed from that text, label line and continuation line alike -- and derives
every quantity it compares against the paper's statements. All arithmetic is over `int`; no
floating-point value is computed anywhere, so no decision depends on one. It runs in about 0.2 s.

What the 59 checks cover, by step:

1. the printed table parses to 13 classes labelled `C00..C12`, nine triples each, every entry in
 `0..26`, every triple three distinct points printed in increasing order, 117 triples in all and
 pairwise distinct;
2. **Theorem 1(a)**: each of the 13 lists is a partition of the 27 points; the `117*3 = 351` pair
 slots are 351 *distinct* pairs and hence all of `C(27,2)`, each once; the design is an STS(27) and
 the 13 classes are a resolution of it; `13 = (27-1)/2`; every point is in exactly 13 blocks and
 in exactly one block of every class;
3. **Theorem 1(b)**: exactly seven blocks lie inside `U`, they are the seven the theorem prints,
 they cover all 21 pairs of `U` once, each point of `U` is on three of them and any two meet in
 exactly one point (the Fano plane); they lie one each in `C00..C06` in the printed order; and no
 block meets `U` in exactly two points;
4. **the Remark's census**: `C00..C06` have profile (1 inner, 4 one-point, 4 disjoint), `C07..C12`
 have (0, 7, 2), the totals are `7 + 70 + 40 = 117`, and the 70 is checked point by point (13
 blocks per point of `U`, three of them inner);
5. **sigma**: a permutation of the 27 points, order 3, fixing exactly `{0,7,8}` with orbit shape
 `1+1+1+3*8`, preserving the block set; the fixed set `{0,7,8}` is itself a block and lies in
 `C01`; `sigma` preserves `U` and acts on it as `x -> 2x mod 7`; it carries each class onto a
 single class, and the induced permutation of the 13 classes is **exactly the `pi` printed in
 Section 4**, of order 3, fixing only `C01`;
6. **the compact form**: `pi` has five orbits and the five classes named in Section 4 meet each
 exactly once; applying `sigma` to those five rebuilds all 13 classes, and the rebuilt design
 equals the printed one class by class;
7. **Proposition 3**, re-derived rather than quoted: `C(27,3) = 2925` point triples are closed under
 the line map, each closure is checked to be closed, every closed set of size `>= 3` is checked to
 be a subdesign (`k(k-1)/6` inner blocks covering its pairs once), and the sizes found are
 **117 of size 3, one of size 7, one of size 27 and nothing else** -- so exactly one order-7
 subdesign, namely `U`, and none of order 9 or 13. The Fano-count-through-a-point invariant is then
 computed (1 on `U`, 0 elsewhere) and Lemma 2's use is checked: the inequality excludes a proper
 subdesign of order 7, 9 or 13 inside one, which is why closing
 triples is *exhaustive* for those orders rather than merely a search;
8. admissibility arithmetic: `27 = 3 mod 6`, `7 = 1 mod 6`, `27 >= 2*7+1`, and `3` does not divide
 `7` -- of which the paper states only the last (Section 1) and the inequality (Lemma 2), the two
 congruences being recorded in the transcript as background.

Every value the paper prints is recomputed and compared, so each comparison can actually fail.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 671e216725054c27212e87a97284556077e14edc222ce6699d433cc2e7836f3d

`verify.py` was written for this folder and run here, on the machine that assembled it, under
Python 3.9.6; it exited 0.

**Provenance of the search that produced the design** (recorded because it is *not* what `verify.py`
re-runs, and because it has one real gap). The files named below are **not shipped in this folder**
-- the folder ships exactly the five files listed under *Contents*. They live in the underlying run
record, and the digests are pointers for matching them there, not
files a referee can hash here. Nothing in this folder depends on them: the design is printed in full
in Table 1 and `verify.py` checks the printed design.

Per the artifacts manifest of the underlying record:

- The design was produced by `pres3.py` (SHA-256
 `30bd72c1b09a7c381e12d8a6ebd509969fdb4e4bdc7c5002b298e744aa5543ff`, 7,069 bytes), a CP-SAT model
 with the automorphism `sigma` and the class action `pi` **prescribed**, invoked as

 the fleet dispatch script, detached, with automatic slot selection, on `pres3.py` in a scratch
 directory, with the timeout argument 2900 and the job tag kts27pres

 from the repository root with a 128-vCPU request, on a fleet slot, timeout 2900 s.
 Solver: OR-Tools CP-SAT, with `num_search_workers = min(96, cpu_count())`.
- Two independent stdlib-only programs on the object are also in that record and were re-run in the
 document stage to capture their output: `witness.py` (SHA-256
 `9a11882847677fe4f1e39a5962f88e217c7c723a4e3ff98b10c4b6416174fb81`, 4,266 bytes, output
 `witness.out`, SHA-256 `1ba42b42edbd2f48e4d94110a29bbeb497212069adebb7e1ddf456dba2a78a84`, ending
 `ALL CHECKS PASSED`) and `recon.py` (SHA-256
 `249e717548752d8a14c9e0138db5cc315db46b44807b325071116ddd8cbfaf95`, 778 bytes, output
 `recon.out`, SHA-256 `3defe2a1a4d8a6eef9e09c2f77dcba1a318e17b786512cd25b8293243a910edd`), which
 rebuilds all 13 classes from the compact certificate. `verify.py` in this folder **shares no code
 with either**; it was written fresh against the printed table, and it re-derives everything they
 do and more.
- The machine-readable forms of the design in that record are `witness.json` (SHA-256
 `c0c1daca70d855ff94ab2e14ee49d6805ccdf6c8fcf5b351e3e4cb1199e277b2`, 4,027 bytes), `witness.txt`
 (SHA-256 `926a0deb13dba37cade9388bfdf01ad50f7e573066a6b3b17482e2b1f093b25b`, 1,366 bytes) and
 `certificate.txt` (SHA-256 `4dc573836f699bb619b6395711726457c02634e527c1b1f9ebd3dd19993a2489`,
 965 bytes, the compact form: `sigma`, `pi` and the five base classes). Table 1 of the paper is
 that same design; the two encodings were checked label-equal class by class and block by block in
 the document stage, and a reader parsing `certificate.txt` alone rebuilt all 13 classes with label
 equality against `witness.json`.

Five gaps in that record, stated rather than dressed:

1. **The search job's own stdout is not in the record.** It was never filed at dispatch (the log
 printed `ARTIFACT_NOT_FILED=no wave`), the slot has since self-terminated, and the off-box copy
 `(an internal object-store copy)` could not be read
 from the control plane (the credential check returned `ExpiredToken`,
 `head_object` returned HTTP 400). The CP-SAT figures quoted in the underlying record -- SAT in
 **29.7 s** on the target with **12,681** orbit variables, and SAT in **18.2 s** with **12,753**
 for the positive control -- therefore survive only as quotations in that record, **not** as a
 captured log. They appear nowhere in the paper, and nothing in the paper or in `verify.py`
 depends on them.
2. **The instance type is not recorded** -- the dispatch log says only that the slot was chosen
 automatically `(... preferring >=128 vCPU)` and the box has self-terminated -- and **total job wall time was not captured**.
 The 2900 s timeout was not hit.
3. **The OR-Tools version is not recorded.** `pres3.py` pip-installs `ortools` unpinned if the
 import fails.
4. **No random seed is set** in `pres3.py`. SAT/UNSAT is deterministic, but *which* satisfying
 design the solver returns is not, so re-running the search need not return this design. That
 costs nothing: the design returned is banked verbatim, printed in Table 1, and independently
 re-verified here.
5. The artifacts gate reported `UNCHECKABLE` on two entries of that manifest, because its
 cross-checker decodes graph6 strings and edge lists only and cannot read a block-list encoding.
 That is "could not check", not "checked and wrong"; the substitutes were the two label-equality
 checks named above, and in this folder the substitute is `verify.py`, which parses the paper's own
 table.

The controls run alongside the search are recorded there too and are worth a referee's eye, since
the trustworthiness of a SAT answer rests on them: the same prescribed-symmetry code path returned
**INFEASIBLE** where the answer is provably NO (`sigma` of order 3 fixing `U` pointwise) and
**SAT** on a positive control whose answer is known (a KTS(27) invariant under one fixed-point-free
translation of order 3 of GF(3)^3, witnessed by AG(3,3)); the plain model reproduced the source's own
settled cases (KTS(15) with a sub-STS(7), SAT) and returned INFEASIBLE on three cases barred by
Doyen-Wilson or by divisibility. None of that is shipped here and none of it is needed: the object
is printed and checked.

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE: this program certifies exactly ONE object -- the KTS(27) of Table 1 of the paper, copied by hand into TABLE above -- and the following quantities, which were inventoried by hand from the paper when this program was written: the resolution, the 351 pairs, the seven inner blocks and their Fano structure, the 7+70+40 census, sigma and its induced class permutation, the reconstruction from five base classes, and the subdesign census behind Proposition 3. NOT ESTABLISHED HERE: no check below opens paper.tex, so two things are asserted rather than verified -- that the table above is the paper's Table 1, and that the list just given is COMPLETE for what the paper states about the object; a quantity the paper states and that list omits would go unnoticed. NOT COVERED: (a) any COUNT of KTS(27) with an STS(7) subdesign -- one object is exhibited and no enumeration or isomorph rejection is performed, so nothing here bounds how many exist or claims minimality or canonicity; (b) the ORDER of Aut(D) -- only that it contains sigma and, by the invariant above, that it is not transitive on points; (c) the orders v = 33 and 39 of the quoted sentence, and Open Problem 4 itself, which are untouched; (d) the BIBLIOGRAPHIC claims -- the byte locator of the quoted passage in arXiv:2110.07874v1, the sentence quoted from Colbourn-Magliveras-Mathon, and the question of priority, which the paper's Scope section expressly declines to claim -- none of which this program fetches or re-reads; (e) a second witness recorded in the underlying run but NOT printed in this paper, which no check here touches. NOT RE-RUN: the constraint search that produced the design -- this program checks the object, not the search, and the paper needs only the object.

Beyond that:

- **A second witness exists in the underlying record and is deliberately not exhibited here.** It is
 a starter-adder construction over `Z_13`, and its subdesign census is *different* (26 sub-STS(7)s
 and a sub-STS(13)), so **no invariant in this paper transfers to it** -- in particular
 Proposition 3 is about the printed design only. It is not printed, not checked, and not relied on.
- **Bibliographic limits, recorded because the paper claims no priority or novelty.** The quoted passage of Section 1 is byte-exact for
 the e-print arXiv:2110.07874v1 (single gzip member; the expanded source is 87,094 bytes, SHA-256
 beginning `dd556250`; the passage occupies bytes 61,204-61,733, lines 1139-1143). The **Wiley
 version of record was not obtained and was not diffed against it**, so the quotation is verified
 for the preprint only, and the two versions carry different titles. The nearest paper, Kokkala and
 Östergård, *Kirkman triple systems with subsystems* (Discrete Math. 343 (2020) 111960) -- which is
 reference 12 of Dukes-Lamken -- has been read in **abstract and metadata only**: three attempts at
 the full text returned HTTP 403 (the ScienceDirect PDF twice, the Aalto accepted manuscript once),
 and `export.arxiv.org` has no record of it, so a remark about order 27 in its concluding section
 would be **unseen rather than absent**. Residual risk is judged low, since Dukes and Lamken cite it
 and still wrote in October 2021 that (27,7) appears difficult. **Section 19.7 of Colbourn-Rosa**,
 the very section Dukes and Lamken name for the problem's status, is a printed book and was **not
 read**. **MathSciNet was not consulted at all** (no institutional access). One citing item of the
 source is **unidentified**: Crossref reports `is-referenced-by-count` 3 against OpenCitations' 2
 and Semantic Scholar's 1, and the OpenAlex citer query never answered (HTTP 429 then a fetch
 deadline, three times), so a citer indexed only there would have been invisible. What *was* run:
 Semantic Scholar on the DOI and on the arXiv id separately, on the nearest paper's record, and on
 a keyword search (772 hits, top 50 read); OpenCitations on the target and on the nearest paper's
 whole citation cone; Crossref on the target, both identified citers, and all five unidentified
 ancestor citers; five zbMATH searches with every hit enumerated, including author sweeps of both
 authors 2021-2026; two arXiv API queries; and the AMS free back file, from which
 Colbourn-Magliveras-Mathon was read directly off pages 1-4 of the
 1992 PDF.
- **The Stinson reference carries a title that the run record does not itself pin.** The record fixes
 the journal, volume, year and pages (Discrete Math. 92 (1991) 371-393) and the fact that Open
 Problem 4 lives there, quoting Dukes and Lamken's own sentence to that effect; the article title
 in the bibliography is standard bibliographic knowledge and was not re-fetched. Nothing in the
 result depends on it.
- **Five automorphism lemmas in the underlying record are not used here and their novelty is
 withdrawn.** They were written up as new, no literature search was run on them, automorphisms of
 STS(27) are classical, and Colbourn-Rosa §19.7 is unread; two of them were also stated too
 broadly (they concern *resolution-preserving* automorphisms). None of them appears in this paper,
 and the closure does not use them. Lemma 2 here is elementary and standard and is proved in place,
 not claimed as new.
- **One figure from an earlier stage is deliberately not carried.** The underlying record mentions
 587 parallel classes in the underlying STS(27); that number was not re-derived in the document
 stage and is not stated in the paper or checked by `verify.py`.
