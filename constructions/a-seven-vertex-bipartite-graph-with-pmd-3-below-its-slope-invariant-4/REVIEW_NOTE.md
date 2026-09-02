# A Seven-Vertex Bipartite Graph Whose Positive Matching Decomposition Number Lies Below Its Slope Invariant

`a-seven-vertex-bipartite-graph-with-pmd-3-below-its-slope-invariant-4`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

Farrokhi D. G., Gharakhloo and Yazdan Pour (*Positive matching decompositions of graphs*,
Discrete Appl. Math. **320** (2022) 311-323, arXiv:2110.12168) define the positive matching
decomposition number `pmd(G)` and, for a bipartite graph with a labelled bipartition, an
invariant `K(G)` obtained by minimising the number of parts their Algorithm 1 returns over all
`m! n!` labellings. They prove `Delta <= pmd <= K` and ask, in the Question immediately after
that corollary, whether the right-hand inequality is ever strict. The paper answers **yes**:

- **the witness.** `Gamma = Theta(2,2,4)`, seven vertices, eight edges, `Delta = 3`: two
 degree-3 vertices joined by internally disjoint paths of lengths 2, 2, 4. `pmd(Gamma) = 3`
 and `K(Gamma) = 4`.
- **`pmd = 3` needs no computation.** The paper gives the three-part decomposition explicitly
 and proves it positive from a lemma proved there in five lines (a matching whose induced
 subgraph is acyclic is positive): `Gamma[V(M_1)]` is a 4-vertex path, and `Gamma - M_1` is
 itself a path on all seven vertices, so its two alternating classes are positive. The lower
 bound is `pmd >= Delta = 3`, which is one sentence.
- **`K = 4` is a finite exhaustion and nothing else.** Algorithm 1 is deterministic once a
 labelling is fixed, so `K` is a minimum over `3! 4! = 144` runs on an 8-edge graph (288 if
 the two sides may also be exchanged). The paper reports the full multiset of outputs
 (`{4,5}` on the printed orientation, constantly `4` on the exchanged one). There is **no
 human proof of `K >= 4`** here and the paper says so.

Both readings of `K` (sides fixed, sides exchangeable) are given, so the asymmetry of the slope
`j - i` does not leave a gap in the claim. The program goes further than the paper does: its
Steps 7 and 8 compute a census of all bipartite graphs of order at most 7, the minimality and
uniqueness of `Gamma` among them, and a second, cubic witness `M_5` on ten vertices. The paper
asserts none of that, and neither does this note.

`verify.py` also records the **near miss inside the source itself**: the source's own second
worked example is a 6-vertex bipartite graph for which it computes `K(tau) = 4` and `pmd = 3`,
so a strict gap at a *fixed* labelling has been in print since 2021. That graph is not a
witness for the Question, which is about the minimised `K(Gamma)`: over all 36 labellings (72
with the sides exchanged) its least output is 3, so `K = 3 = pmd`. The paper makes no claim
about it.

## What was checked, and how

Everything except the two exhaustions can be followed by hand from the objects printed in the
paper. `verify.py` re-derives all of it mechanically, in exact integer and bitmask arithmetic
-- there is no floating-point value anywhere in the program, so no decision depends on a
tolerance. In particular it checks the two things a referee is most likely to doubt:

1. **that our Algorithm 1 is the source's.** The source's pseudocode is now *reproduced in the
 paper*, statement for statement, beside the streamlined offer-and-accept form the rest of the
 paper uses, and the program runs **both**. Each reproduces the parts the source prints for
 each of its two worked examples, edge for edge (checks
 `example_i_parts_match_the_source_edge_for_edge`, `example_i_literal_pseudocode_also_...`,
 and the same pair for example (ii)); and the two forms agree part for part on every bipartite
 graph of order at most 6 under every labelling together with every labelling of `Gamma`, the
 two source examples and `M_5` -- 82,402 runs, 0 disagreements
 (`alg1_equals_the_literal_pseudocode_part_for_part`). The program also verifies that the
 algorithm is insensitive to the order in which equal-slope offers are scanned, so `K(tau)` is
 well defined, and that every one of 648 runs returns an ordered partition of the edge set into
 matchings. What is left is the one thing no program here can do: confirm that the displayed
 statements are the ones printed at lines 549-574 of the source file. The locator is exact and
 that comparison is a referee's.
2. **that our reading of "positive matching" is not doing hidden work.** The criterion is
 implemented twice, once literally over *all* non-empty subsets of the matching and once as a
 pendant-edge peeling, and the two are compared on 678 matchings across nine graphs (0
 disagreements). As anti-controls, no perfect matching of `K_{3,3}` comes out positive, and
 the source's published values `pmd(K_{m,n}) = m+n-1` are reproduced for all ten pairs
 `1 <= m <= n <= 4`.

The remaining checks re-derive: the decoding of both graph6 strings and their agreement with
the printed edge lists and bitmasks; the identification of `Gamma` as `Theta(2,2,4)` from the
degree sequence and the three internally disjoint paths; the positivity of `M_1, M_2, M_3` in
their successive residuals and the exact value `pmd(Gamma) = 3` by exhaustive search; that
`Gamma` has eight matchings of size 3 and that none is positive; `K(Gamma) = 4` over 144 and
over 288 labellings, together with the parts Algorithm 1 returns on the printed labelling; the
five anti-controls `pmd = K` (`K_{2,2}`, `C_6`, `K_{3,3}`, `K_{1,4}`, and the source's second
example). Beyond anything the paper claims, it also computes the whole census table in both
readings, with the witness classes identified by canonical form; the published bound
`Delta <= pmd <= K` on all 215 and 407 classes of order at most 7; and the properties of `M_5`,
including that the rotation by one is an automorphism exchanging its two sides.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only (`itertools`, `sys`): no third-party package and no
external data file. Runtime about 35 seconds on one core. The program prints one line per
check and a closing verdict, and exits 0 only if every check passes. The recorded run reports
**80 checks, all passing**:

 VERDICT: ALL 80 CHECKS PASS

Its inputs are literals in the program. Those bearing on the paper are exactly the objects
exhibited there -- the row bitmasks `(3,13,14)`, the printed 8-edge list, the graph6 string
`FEhb?`, the three parts `M_1, M_2, M_3` -- and it also carries, for its own further steps, the
two worked examples of the source with the parts they print, the census table, the graph6
string `IhEIHCPaG` and the edge list of `M_5`.

## Provenance

**Of this folder's program.** `verify.py` was written for this folder and run locally on this
machine; it is **not** one of the scripts that produced the result. `verify.output.txt` holds
its output, preceded by a provenance header and followed by an exit status, both written by the
run harness. The header records the SHA-256 of the program that produced the output, so the two
files can be matched:

```sh
shasum -a 256 verify.py
```

 2abf3f970fd8dcb86f2245db2d3489387223abf79329655c922b52052e9f7678

The recorded exit status is 0. `verify.py` reads nothing outside itself, so a referee needs
neither this folder's other files nor any of the material described next.

**Of the original computation.** The result came out of an internal run whose artifacts are
indexed with per-file SHA-256 digests in that run's manifest (19 files, all present, all
digests re-hashed). Stated as that manifest records it, and not improved upon:

- the decisive script is a from-scratch transcription of Algorithm 1, of the positivity
 criterion and of the `K` minimisation, written from the source's LaTeX; it computed the
 headline values, the source's Example (i)/(ii) controls, `M_5`, and the census up to
 isomorphism. It ran on one AWS slot instance under SSM, status `Success`, `RC=0`, python 3.9,
 one vCPU, no randomness; its complete standard output (3,390 bytes of program output inside a
 4,442-byte capture) is archived beside it, and the numbers in that output are the numbers in
 this paper. A second archived script re-ran the census under a per-orientation reading and
 tested the run's own auxiliary lemma exhaustively; its output is archived too.
- **honest gaps, recorded in that manifest rather than papered over.** Seven further scripts of
 that run -- including the 22-minute census job and a second, independently written census --
 are indexed as code with **no captured output**: the standard output was never written to a
 file, the slot's S3 upload is verified *absent* (HTTP 404 on all three keys), and the AWS
 credentials needed to retrieve anything by command id had **expired** before the manifest was
 written. Those outputs are therefore *unread, not absent*. Two of the indexed scripts were
 dispatched from a scratch directory rather than from the run tree, so the harness did not file
 them at dispatch; the copies indexed are the files at the exact paths that were dispatched,
 and they could not be re-verified against the payload the dispatch carried. Elapsed runtimes
 were recorded for only one of the scripts.
- one computation in that run *was* randomised -- a non-exhaustive probe of larger cubic graphs
 by 128,000 random labellings, which yields upper bounds on `K` only. **Nothing in this paper
 uses it**, it is not indexed as an artifact, and `verify.py` does not reproduce it.

**Independence, and one bookkeeping difference.** Because those outputs are partly unread,
`verify.py` was written from scratch for this folder rather than derived from the original code,
and it reproduces every number this paper prints, so the paper does not rest on the internal
record. Where the two can be compared they agree: `K(Gamma) = 4` over 144 and 288 labellings,
the unique order-7 witness, `pmd(M_5) = 4` with `K(M_5) = 5`, and the order-7 class counts
restricted to graphs with no isolated vertex (`(1,6):1`, `(2,5):11`, `(3,4):42`, total 54). Two
conventions differ and a referee comparing numbers should know it: our census allows isolated
vertices and so counts 128 classes at order 7 where the internal run, which excluded them,
counted 54; and in the per-orientation reading our table counts *ordered* bipartitions and
therefore lists 2 witnesses at order 7, which `verify.py` confirms are the two orientations of
the single graph `Theta(2,2,4)`, reported as one class internally.

No reproduction instruction beyond `python3 verify.py` is offered here. In particular this note
does not tell a referee how to re-run the original census, because the record does not support
such an instruction: its output was not preserved.

## Scope

What the paper settles is the Question as printed -- some bipartite graph has `pmd < K` -- by the
single witness `Gamma = Theta(2,2,4)`. What is **not** settled:

- **no structural lower bound on `K`.** Every `K` value here is an exhaustion over labellings.
 We do not prove, use, or endorse any lower bound on `K` for regular bipartite graphs; in
 particular the paper does **not** claim that every cubic bipartite graph with `pmd = 4` is a
 witness, and says nothing about the circular ladders, the Moebius ladders or the
 generalised Petersen graphs at all. An earlier form of this result did
 advance such an infinite family, on the strength of a hand-proved lemma that no independent
 check ever machine-verified; it is deliberately absent here.
- **orders 8 and above are not surveyed,** so nothing here counts witnesses of order 8 or 9, and
 nothing bounds `K - pmd`: the one witness in this paper has `K - pmd = 1`.
- **Question 2 of the source is untouched.** It asks for a labelling with `K(tau)` below the
 minimum number of slopes, which is a different statement.
- **the combinatorial characterisation of a positive matching is quoted, not proved.** It is a
 theorem of the source; this paper and this program use it as the definition and never touch
 the underlying algebraic one.
- **the *reading* of the source's pseudocode is a human step.** The paper now displays that
 pseudocode statement for statement and `verify.py` runs it as written, so the gap is no longer
 "is our algorithm theirs" but the narrower "are those the statements printed at lines 549-574".
 A referee with the e-print open settles that by eye; the locator is exact and nothing else in
 the folder can do it for them.

The program's own closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE -- what this program does NOT cover. It re-derives, in exact integer arithmetic and from the objects printed in the note, the quantities the note states: the witness Gamma and its pmd and K. Steps 7 and 8 go beyond the note's claims: the census of all bipartite graphs of order at most 7 under both readings of K, the minimality and uniqueness of Gamma among them, and the second witness M_5 are computations of this program, asserted nowhere in the note. NOT RE-RUN here: (a) the source's combinatorial characterisation of a positive matching, which is a theorem of the source and is taken as the definition throughout -- this program never touches the underlying algebraic definition; (b) the READING of the source's pseudocode. Step 2 goes as far as a program here can: the pseudocode is written out statement for statement as `alg1_source_literal`, it reproduces both of the source's worked examples part for part, and it agrees with the streamlined form used everywhere else in this program on every bipartite graph of order <= 6 under every labelling and on every labelling of Gamma, the two source examples and M_5. What no program here can settle is that those statements are the ones printed at lines 549-574 of the source file; that single comparison is left to the reader; (c) every order >= 8 except the single graph M_5, so nothing here says how many witnesses of order 8 or 9 exist; (d) any infinite family, and in particular no lower bound on K for r-regular bipartite graphs is proved or used anywhere in this program.

## Two things a referee should check against the source

**The statement being answered.** It is a numbered `question` environment in the e-print source
of arXiv:2110.12168v2, at lines 635-637 of the single 69,193-byte file
`2022-05-01_PMD__Revised_1_.tex`, byte-identical to lines 497-499 of v1. It carries the number 1
because the environment has its own counter with no per-section reset and the first of the three
occurrences of `\begin{question}` in the file is commented out. The paper reproduces the
sentence verbatim, so what is answered does not depend on the label. The sibling Question 2, at
line 639, is a different statement and is not addressed.

Every part of that locator has been re-checked against the e-print as retrieved: the file is
69,193 bytes and 783 lines; the three `\begin{question}` occurrences are at lines 457
(commented out), 635 and 639; the preamble declares `\newtheorem{question}{Question}`, so the
counter is the environment's own and is not reset per section; and lines 635-637 of v2 are byte
for byte lines 497-499 of v1. The one liberty the paper takes with the quoted sentence is
typographic and is disclosed in the paper: the source sets the invariant in fraktur,
`\mathfrak{K}`, and the paper sets it upright. Nothing else in the sentence is altered.
Algorithm 1 is quoted from the `algorithmic` block at lines 549-574 of the same file.

**The journal reference.** Volume, year and pages (Discrete Appl. Math. **320** (2022) 311-323,
doi `10.1016/j.dam.2022.05.012`) are taken from Crossref, not from the e-print, which carries no
journal reference; the printed journal text was not accessible here, so statement numbering is
quoted from the e-print alone.

**Prior art, for the record.** `K(Gamma)` is defined four lines above the Question, so it did
not exist before October 2021 and any prior observation of `pmd < K` would have to cite the
source. The citing set was enumerated through Semantic Scholar, OpenCitations, OpenAlex, zbMATH
and the arXiv API and is very small (Semantic Scholar returned three citing papers,
OpenCitations one, OpenAlex a cited-by count of two). Four e-prints in that neighbourhood were
read as LaTeX source rather than as PDF and grepped: none contains the word "slope" or the
symbol used for `K`. One channel is **unread rather than empty**: an OpenAlex `cites:` filter query failed on four attempts, and the authors' own GAP
data file on their institutional page could not be fetched (TLS version rejected). MathSciNet
was not consulted. The one thing found that comes close is the source's own second example,
which `verify.py` examines in full.
