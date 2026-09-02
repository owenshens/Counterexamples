# Under a simplicity hypothesis, the order of a [3,1;4]-mixed cage is twelve

`a-simple-3-1-4-mixed-cage-of-order-twelve-answering-araujo-pardo-et-al`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

Araujo-Pardo, Hernández-Cruz and Montellano-Ballesteros, *Mixed cages*, Graphs Combin. **35**
(2019) 989–999 (arXiv:1702.07255), prove `n[z,1;4] <= 3(z+1)` for odd `z` and `3z+2` for even
`z`, add "Moreover, the equality holds for `z ∈ {1,2}`", and in their Section 4 propose as
their **first problem** the search for a matching lower bound. The paper quotes both passages
verbatim with line and byte locators into the `/e-print` source, whose SHA-256 it also prints.

The paper computes the single case `z = 3` under the standing hypothesis that mixed graphs are
simple: **for simple mixed graphs, `n[3,1;4] = 12`**, so the odd-`z` bound is tight at `z = 3`
for simple mixed graphs. Whether that is the same assertion as an extension of the source's
equality clause from `{1,2}` to `{1,2,3}` depends on the unresolved simplicity convention, and
the paper makes no claim to settle or extend the source's problem. Two components:

* **Upper bound — not ours.** `n[3,1;4] <= 12` is the source's own Lemma 3.1(1) at `n = 3`,
 the mixed circulant on `Z_12` with arcs `i -> i+1, i+2, i+7` and edges `{i,i+6}`. The paper
 prints its 36 arcs and 6 edges in full.
* **Lower bound — the new content.** An elementary lemma (**Lemma L**): a digraph on `m`
 vertices with no loop, no digon, no directed triangle and minimum out-degree `δ >= 1` has
 `m >= 2δ + 2`. Partitioning the vertex set around one vertex `u` traps `N⁺(a)` inside
 `(A \ {a}) ∪ S` for every `a ∈ A = N⁺(u)`, so Lemma L applies to the digraph induced on `A`,
 which has only `z` vertices. At `z = 3`, order 10 forces `|S| = 2`, hence a 3-vertex
 digon-free triangle-free digraph with minimum out-degree `>= 1`, which Lemma L forbids. In
 general this gives `n[z,1;4] >= ⌈(5z+6)/2⌉` for `z >= 2`.

There is **no solver, no case split and no census anywhere in the proof**, and the paper claims
no novelty for Lemma L or for the count, recording that counts of that kind are standard
(Behzad–Chartrand–Wall).

## What was checked, and how

The decisive argument is short enough to follow by hand, and a referee should do so: two
displayed inequalities kill order 10, and six modular computations validate the order-12
witness. `verify.py` mechanises all of it anyway, and adds independent routes where one
exists:

* It **parses the 36-arc/6-edge link list transcribed by hand from Section 2 as text**, rebuilds the same
 object independently from the circulant description in that section, and checks the two
 label-equal.
* Every girth is computed **twice**: by definitional enumeration of closed walks on distinct
 vertices, and by a breadth-first search over links that forbids reusing one link twice in a
 row. The two routes are required to agree.
* **Lemma L is confirmed by complete enumeration**: a loopless digon-free digraph on `m`
 labelled vertices is exactly a choice of `{no arc, i->j, j->i}` per unordered pair, so all
 `3^C(m,2)` of them are visited for `m <= 5` (59,049 at `m = 5`). This covers the only
 instance the lower bound uses (`m = 3`, `δ = 1`, where no such digraph exists) and its tight
 case (`m = 4`, the directed 4-cycle).
* The exclusion of the orders 6, 8 and 10 at `z = 3` is then re-derived **from that
 exhaustive table**, not from the printed proof, and the program also checks that the
 argument does *not* exclude order 12 — a lower bound that refutes its own witness is a bug.
* Forced positives: the source's other constructions `C_6(1)`, `C_8(1,5)` and
 `C_14(1,2,8,9)` are validated as `[1,1;4]`-, `[2,1;4]`- and `[4,1;4]`-mixed graphs of girth
 exactly 4; and Theorem 3's formula at `z = 2` returns 8, the published exact value.
* Anti-control: the multi-arc order-6 object of the paper's scope section is confirmed to have
 girth 4 under a multigraph reading **and** to be rejected as a simple `z = 3` object (its
 out- and in-degrees collapse to 1). That is exactly why simplicity is printed as a
 hypothesis rather than assumed.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
All arithmetic is Python integers; the single non-integer constant in the paper (Shen's 2.885)
is carried as `fractions.Fraction(2885, 1000)`, and no floating-point value is compared,
rounded or thresholded anywhere. Runtime is about 0.2 s. The program prints one line per
check and a closing verdict, and exits 0 only if every check passes. The recorded run reports
**66 checks, all passing**:

 VERDICT: ALL 66 CHECKS PASS

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 2ff23e55cb62afb292974afb9c87fa80839438d17affd656f0b36f3a7d5f19fa

`verify.py` is **new to this folder**: it was written at the hand-over stage and is not one of
the programs that found the result. Its provenance is the command above, run on the laptop
control plane with the output captured by the harness; Python 3.9.6, exit status 0.

The programs that ran while the result was being found are **not** shipped here, because no
claim in the paper depends on them. Quoting the row's artifacts manifest rather than
reconstructing it, so a referee is not misled about how re-runnable they are:

* Three scripts are on disk in the row's artifacts directory with SHA-256 and byte counts —
 `mixedcheck.py` (an independently written checker, 7,663 B), a driver that runs only its
 cheap part, and a witness cross-checker — together with the captured stdout of the two
 local runs. Those two local runs reproduce the five-object control table and the
 witness-labelling equality check, and nothing more.
* `mixedcheck.py`'s exhaustive census (its "PART B") is indexed in that manifest as code
 **with no output**. It was dispatched twice and *never collected*: once during the original
 run (slot S14, still `Pending` — SSM never
 delivered it) and once at the document stage (slot S07, timeout 3000 s, still `InProgress` at
 hand-off). **No stdout of either dispatch exists in the record.**
* Two earlier confirmation jobs are cited in the row's own note by slot and command id with
 `RC=0` — including a census reporting 0 labelled `[3,1;4]`-mixed graphs of order 10 over
 3,136 shards — but **their scripts are not on disk and are named nowhere in the record**,
 and the boxes have self-terminated. They are reported, not shipped, and a referee cannot
 re-run them from this folder or from the artifacts directory.

None of that is load-bearing, which is the reason the folder is shippable in spite of it: the
lower bound is a hand proof that excludes order 10 in two lines, and `verify.py` re-derives
the same exclusion from an exhaustive enumeration small enough to run in a fraction of a
second on a laptop.

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE. This program checks the note's own claims and nothing else. It reads the
> 36-arc/6-edge link list transcribed by hand from Section 2, rebuilds it independently from
> the circulant description, and re-derives every difference, sum, girth and bound; Lemma L is
> confirmed by complete enumeration of all loopless digon-free digraphs on at most 5 vertices,
> which covers the m = 3 instance the lower bound uses. GAPS NOT COVERED: this program never
> reads the paper, so the agreement of the transcribed link list with the printed one is not
> checked here; and the published constants n[2,1;4] = 8, n[z,4] = 3z+1 and n_AHM[1,1;4] = 6
> are QUOTED, not derived here, so part 8 and part 9 are arithmetic consequences of the
> literature rather than independent verifications of it; nothing here decides the multigraph
> variant beyond exhibiting an order-6 object for it; nothing here claims the order-12 witness
> is unique; and nothing here touches girth g != 4 or edge-degree r != 1.

Four further limits, stated in the paper and repeated here because they are what a referee
should press on.

1. **Simplicity is a hypothesis, not a convenience.** The source contradicts itself: line 65
 of its `/e-print` source says "We allow multiple edges and arcs", line 67 defines a mixed
 regular graph as "a simple and finite graph `G`". Under the multigraph reading the value 12
 **fails** — the paper exhibits an order-6 object (arc `i -> i+1` tripled plus edges
 `{i,i+3}` on `Z_6`), so there `n[3,1;4] <= 6`. The paper does not resolve which reading the
 source's problem intends; the multigraph variant is undecided here.
2. **A must-cite whose printed hypothesis-free form fails.** The unrefereed preprint arXiv:2401.14768 states
 `n[z,r;g] >= n[z,g] + n_AHM[1,r;g] - g` with no hypothesis, which substitutes at `(3,1;4)`
 to exactly our 12. A referee holding it reaches our number in one step, so it must be
 cited — and the failure of that printed hypothesis-free form disclosed in the same breath: it gives `>= 9` at `z = 2` against
 the published exact 8, and `>= 15` at `z = 4` against the source's own order-14
 construction, which `verify.py` validates. The paper does this; it does not rely on the
 preprint for anything.
3. **The general bound is elementary, not numerically best.** Deleting the perfect matching
 reduces the problem to a triangle-free digraph, so Shen's Caccetta–Häggkvist bound applies
 directly and, rounded to an even order, *beats* Theorem 3 at `z = 10` and `z = 12`. The
 comparison table in the paper is reproduced check-by-check by the program. Theorem 3
 strictly improves what was available at exactly one value of `z` — the target cell `z = 3`.
4. **One new cell only.** `n[4,1;4] = 14` follows from the published `n[4,4] = 13`, parity,
 and the source's Lemma 3.1(2); it is new here only in being written down. Every `z >= 5`
 stays open, e.g. `n[5,1;4] ∈ {16,18}`.

Four attributional risks remain that no program can close, and Section 5 of the paper states
its novelty claim as a bounded search for exactly this reason.

1. The body of a 2025 AWM-series survey chapter by the source's first author is paywalled and
 unread. Its 59-entry reference list *was* read and holds no girth-4 mixed-cage item, which
 bounds the risk without eliminating it.
2. MathSciNet was not consulted, so a review of the 2019 paper carrying a remark on the
 `z = 3` cell would not have been seen.
3. One of the OpenCitations citers of the 2019 paper — a 2026 row — has an empty citing DOI
 and was never resolved to a document.
4. The OpenAlex citation query returned HTTP 429 and then timed out. It asks the same question
 as the Semantic Scholar citer channel, which answered in full with 12 citers, so the
 citation graph is not unread; but that specific channel is a hole.

The novelty claim is therefore "no correct published proof of `n[3,1;4] >= 12` was found, and
the value 12 is printed nowhere we could reach", not a claim about the whole literature. The
two channels that carry the most weight both answered and both came back empty: the zbMATH
search for "mixed cage" contains the string "girth 4" zero times across 26 titles, and the
survey's 59-entry bibliography holds only four mixed-cage items, none of them at girth 4.

Separately, the record is explicit about which numbers are taken on trust. The transcript's
scope note names the three — `n[2,1;4] = 8`, the
girth-4 Behzad–Chartrand–Wall value `n[z,4] = 3z+1`, and `n_AHM[1,1;4] = 6` — and records
that we hold no primary reference for the second beyond arXiv:2401.14768's assertion that it
is proved for degree `z ∈ {2,3,4}`. `n[3,1;4] = 12` itself uses none of them; the side claim
`n[4,1;4] = 14` and the `z = 4` comparison entries do.
