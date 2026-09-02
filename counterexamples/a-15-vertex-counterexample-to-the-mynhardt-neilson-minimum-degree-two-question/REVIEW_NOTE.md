# A 15-Vertex Graph of Minimum Degree Two with alpha_bnr > Gamma_b

`a-15-vertex-graph-of-minimum-degree-two-with-alpha-bnr-greater-than-gamma-b`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

Those four files, plus this note, are the whole folder. Nothing in the paper or in this note requires a
file that is not listed above: where the paper refers to work done outside it --- the
1,026,505-graph 2-connected census, and the values for the family `G_k'` with `k >= 2` ---
it says so and marks the work as not shipped, and no claim of the paper rests on it. The
figures in question, recorded here so that the paper's disclaimer is not a pointer to
nothing, are `Gamma_b(G_k') = 2k+3` and `alpha_bnr(G_k') >= 3(k+1)`, giving an unbounded
`alpha_bnr - Gamma_b` on `delta >= 2`; the paper explicitly neither asserts nor verifies
them.

## What is claimed

The second sentence of Question 4 of Mynhardt and Neilson (*Comparing upper broadcast
domination and boundary independence broadcast numbers of graphs*, Trans. Combin. **13**
(2024) 105--126; e-print arXiv:2104.02257v2) asks whether `alpha_bnr(G) <= Gamma_b(G)` when
`delta(G) >= 2`. The paper answers **no**, with an explicit graph `G'` on 15 vertices and 26
edges: `delta(G') = 2` and `alpha_bnr(G') >= 6 > 5 = Gamma_b(G')`; the exact value `alpha_bnr(G') = 6` is not needed for the theorem and comes from the census.

Three things a referee should note before reading further, all of them stated in the paper
itself.

1. **Only the second sentence of Question 4 is answered.** The first sentence, for
 2-connected `G`, is untouched and remains open. `G'` has four cut vertices, and
 2-connectedness is strictly stronger than `delta >= 2`, so the half that falls here is the
 easier one. The bare phrase "Question 4 is refuted" would misreport this result.
2. **The printed sentence has a defect and the paper discloses it first.** It writes the
 subscripted `G_k` under a hypothesis on an unbound `G`, while the same paper proves
 `delta(G_k) = 1`; read literally the sentence is vacuous. What is answered is one repaired
 reading, "`alpha_bnr(G) <= Gamma_b(G)` for every graph `G` with `delta(G) >= 2`"; a second
 repair --- whether `alpha_bnr - Gamma_b` can be arbitrarily large when `delta >= 2` --- is
 admissible and is *not* settled by the paper, which treats no member of `G_k` with `k >= 2`. The
 published version carries the same subscript at journal page 125.
3. **The one new step is small, and its novelty is not cleared.** `alpha_bnr(G_1) = 6 > 5 =
 Gamma_b(G_1)` is already printed in the source for its own `G_1`; the obstruction the source
 names is that `G_1` has end-vertices. The new content is a lemma saying that doubling an
 end-vertex into a true twin leaves `Gamma_b` unchanged, plus the observation that the
 source's own broadcast survives the doubling verbatim. That lemma was **not** checked
 against the broadcast-domination literature, and a prior form of it there is the likeliest
 way this note loses its content.

## The verification program

What was checked, and how. The decisive inequality does not need a computer. The lower bound `alpha_bnr(G') >= 6` is one
exhibited broadcast whose balls, spheres and private boundaries are printed in the paper and
can be read off the 26-edge list in a few minutes. The upper bound `Gamma_b(G') <= 5` is the
end-vertex doubling lemma, proved in full in the paper, applied four times to the source's own
published `Gamma_b(G_k) = 2k+3` at `k = 1`.

Everything else --- the two exact values, the four counts, and the structural table --- is
checked by the accompanying program.

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
It reads the object exhibited in the paper (the 26 named edges and the graph6 string) and
derives every quantity it compares against the paper's statements. It prints one line per
check and a closing verdict, and exits 0 only if every check passes. The recorded run reports
**53 checks, all passing**, in about 16 seconds:

 VERDICT: ALL 53 CHECKS PASS

In outline it (A) rebuilds `G'` and confirms the printed graph6 string decodes to the same
*labelled* graph, then re-derives the degree and eccentricity tables, `delta = 2`, `diam = 4`,
the cut-vertex set `{y1,y2,z1,z2}` and the broadcast-space size 1,474,560,000; (B) confirms the
doubling structure, including that the doubling moves no distance and no eccentricity;
(C) re-derives every set printed in the proof, and that `v` is the unique unheard vertex;
(D) runs its own exhaustive census and returns `Gamma_b(G') = 5`, `alpha_bnr(G') = 6`,
`alpha_bn(G') = 7` with 283 / 1140 / 1416 broadcasts in the three families, plus
`Gamma_b(G_1) = 5`, `alpha_bnr(G_1) = 6`, `delta(G_1) = 1`; (E) validates its own census
enumerator; and (F) re-runs the doubling sweep over all 22,654 (H, end-vertex) pairs with H
connected on 3..6 vertices, 0 mismatches, graph counts 4, 38, 728, 26704 matching OEIS A001187.

Two of those deserve emphasis, because they are what makes the census evidence rather than an
assertion.

* The census does not walk all 1,474,560,000 broadcasts of `G'`; it prunes on irredundance and
 on bn-independence, both of which are monotone (once violated on a partial assignment they
 stay violated), so no member of either family is missed. To check the *implementation* rather
 than the argument, the same seven quantities are recomputed for the 11-vertex `G_1` by an
 **unpruned** walk of all 3,686,400 of its broadcasts; the two agree exactly. No unpruned walk
 of `G'` itself is performed.
* The definitions transcribed from the source --- in particular the convention that a
 broadcaster of strength 1 lies in its own private boundary, and the characterisation of a
 minimal dominating broadcast as a dominating irredundant one --- are not re-derived. The
 evidence that they were transcribed correctly is that the census reproduces the source's own
 published `Gamma_b(G_1) = 5` and `alpha_bnr(G_1) = 6`, and five further published control
 values (`P_6`, `P_7`, `C_6`, the 3x3 grid, the 3x4 grid). A referee who disagrees with the
 conventions should start there.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 1fcc3c995bdbb27e5573971d384becf3c1ae697426d6133b605b6bc6f8eb62cb

`verify.py` was written for this folder and is a fourth, independent implementation. The values
`Gamma_b(G') = 5`, `alpha_bnr(G') = 6`, `alpha_bn(G') = 7` and the counts 283 / 1140 / 1416
were first produced by three separately written census engines dispatched to remote hosts; the
row's artifacts manifest records, for two of those engines, the script as run with its SHA-256,
the invocation, the EC2 instance and the SSM CommandId, together with the captured stdout of
each (`bcast.out`, `lemma.out`, both recorded as complete and ending in their own job's
done-marker). `verify.py` reproduces all six of those numbers, and additionally the two counts
283 and 1140 that the row's own record flags as having come from a single engine.

Two provenance gaps are recorded here rather than papered over, because the artifacts manifest
records them as gaps:

* **Six of the ten indexed scripts have no captured stdout.** The manifest states that the jobs
 returned success with RC=0 but that their output was transcribed into the row's database note
 and never written to a file, and that the dispatcher's own `head_object` calls returned 404
 for every S3 copy, so no S3 artifact exists to point at either. The claims that rested only
 on those six --- the 1,026,505-graph 2-connected census, the full uncapped census of the
 22-vertex `G_2'`, and the values for `k = 3, 4` --- are therefore **corroboration, not
 certificate**. None of them is a claim of the paper: the first appears in the paper only as
 explicitly-labelled weak evidence for the *open* 2-connected question, and the others only in
 a bullet that says the unbounded-gap statement is not claimed here.
* **One indexed script never returned.** `lemmas.py` was still `InProgress` with
 `ResponseCode=-1` when the attack finished; its result is unread, not negative, and no number
 in the paper comes from it.

Nothing in `verify.py` re-runs any of those jobs, and nothing in the paper depends on them: the
program reads only the 15-vertex object printed in the paper and the 11-vertex `G_1` obtained
from it by deleting four vertices.

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE -- what this program does NOT establish. (1) It settles only the SECOND sentence
> of Question 4 of Mynhardt and Neilson, the delta(G) >= 2 one, in its repaired reading; the
> FIRST sentence, for 2-connected G, is untouched and remains open, and G' is deliberately not
> 2-connected (four cut vertices, checked above). (2) NOT RE-RUN: any claim of minimality -- no
> census of graphs of order below 15 with delta >= 2 is performed here, so nothing above says
> G' is a smallest witness. (3) NOT RE-RUN: the family G_k for k >= 2 -- this program reads only
> the 15-vertex object printed in the paper and the 11-vertex graph G_1 obtained from it, so the
> unbounded-gap statement is outside its scope. (4) GAPS NOT COVERED: the definitions of
> bn-independence, irredundance, alpha_bnr and Gamma_b, and the characterisation of a minimal
> dominating broadcast as a dominating irredundant one, are TRANSCRIBED from the source paper
> and not re-derived [...]. (5) The census of G' uses two monotone prunes rather than a full
> walk of its 1,474,560,000 broadcasts [...]. (6) The end-vertex doubling sweep is exhaustive
> only for H on at most 6 vertices; the general lemma rests on the hand proof in the paper, not
> on the sweep.

Beyond the program's own scope, the paper states two further limits. "Smallest" is **not**
claimed: no census of graphs of order below 15 with `delta >= 2` was run, and `G'` is not the
unique repair --- the same lemma applied to other end-vertexed graphs of the source yields more
witnesses. And the prior-art search behind item 3 of "What is claimed" reached the arXiv, zbMATH,
Semantic Scholar, OpenAlex, DBLP and the published PDF, but **not** MathSciNet (no institutional
access), and it did no bibliography walk of the broadcast-domination literature for the doubling
lemma specifically.
