# A 28-vertex circulant with integral fractional broadcast number strictly between its multipacking and broadcast numbers

`a-28-vertex-circulant-answers-teshimas-integral-fractional-broadcast-problem`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

Problem 4 of the open-problems section of L. E. Teshima, *Multipackings in Graphs*
(arXiv:1409.8057v1), quoted verbatim in the paper, asks whether some graph `G` has
`gamma_{b,f}(G)` an **integer** with `mp(G) < gamma_{b,f}(G) < gamma_b(G)`. The paper
answers **yes** and exhibits the circulant `G = C_28(1,4)`, for which

 mp(G) = 3 < gamma_{b,f}(G) = 4 < gamma_b(G) = 5,

with the middle term integral.

The claim is an existence claim, so one witness settles it. The paper makes **no priority
claim**: it names the closest work it could read (the source survey itself, and F. Yang's
2015 M.Sc. and 2019 Ph.D. theses, whose only fractional examples are integral with
`gamma_{b,f} = gamma_b`), and it records the items it could not read (Hartnell–Mynhardt
*Utilitas Math.* 94 (2014), and Sen–Kola *Discrete Appl. Math.* (2024)). The paper also
records that the dual it uses in certificate (B) is a published theorem *of the survey
posing the problem*, credited there to an unpublished 2013 Brewster–Duchesne manuscript that
we were unable to obtain; that is the largest residual risk and it is attributional, not
mathematical.

## What was checked, and how

Every value in the chain is pinned by a certificate a referee can check by hand from the
ball profile `(n_1,...,n_5) = (5,13,21,27,28)` printed in the paper:

* `gamma_{b,f} = 4` — an explicit primal of cost 4 and an explicit dual of value 4. **Weak
 duality alone** closes it; no LP solver is involved, in the paper or in the program.
* `gamma_b = 5` — an explicit cost-5 broadcast above, and below it a ball-counting bound
 whose five rows at cost 4 are printed in full (max 27 < 28).
* `mp = 3` — an exhibited 3-element multipacking below, and above it a tiling obstruction:
 a 4-element multipacking would force the four sets `{v : d(u,v) >= 4}` to partition
 `Z_28`, which a mod-7 argument rules out because `t(x) = 1 + x + 2x^3 + 2x^4 + x^6` is
 monic of degree 6 and is not `Phi_7`.

`verify.py` mechanises all of that in exact arithmetic and adds two independent routes to
the one step (`mp <= 3`) whose paper argument is algebraic rather than finite: an exhaustive
scan of all 210 residue distributions of a 4-element set, and an exhaustive scan of all
2,925 four-element subsets containing 0 (legitimate because `Z_28` acts transitively, which
the program also checks).

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
All arithmetic is Python integers and `fractions.Fraction`; no floating-point value is
compared, rounded or thresholded anywhere. The program prints one line per check and a
closing verdict, and exits 0 only if every check passes. The recorded run reports
**68 checks, all passing**:

 VERDICT: ALL 68 CHECKS PASS

It reads the objects exhibited in the paper as input — the `graph6` string for the primary
witness, and the arithmetic definition of that witness — and derives every quantity it
compares against the paper's statements. It also reads the graph two ways and checks them
label-equal: the `graph6` line decoded by its own stdlib codec, and the circulant built from
`i ~ i ± 1, i ~ i ± 4 (mod 28)`. Runtime is well under a second.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 97718d60b68df30d89ca93716b4921e5cf16eff9af2934e1d682b0bfb6c0223f

`verify.py` is **new to this folder**: it was written at the hand-over stage, and it is not
one of the programs that found the result. Its provenance is therefore simply the command
above, run on the laptop control plane, with the output captured by shell redirection.

The programs that *found* the result live in the row's artifacts directory
(`runs/wave24/artifacts/t10849/`, MANIFEST.json with SHA-256 and byte counts on 21 entries)
and are **not** shipped here, because nothing in this paper depends on them. What that
manifest records about them, quoted rather than reconstructed, so that a referee is not
misled about how re-runnable they are:

* Six scripts ran on AWS slots via SSM; all six report `ssm_status = Success`, `job_rc = 0`,
 region `us-east-1`, no randomness. Their `invocation` fields are marked **RECONSTRUCTED**
 in the manifest: no `ARTIFACT_INVOCATION` line was preserved in the row, so the
 `<timeout_s>` argument each was dispatched with **is not in the record**. Command ids and
 instance ids *are* in the record, read off each job's own stdout header.
* One of the eight retrieved stdout files, `construct.out`, is **truncated** by SSM's ~24 KB
 cap (its body begins mid-tuple), and no S3 copy exists — the manifest records that a
 listing of 2,046 objects under the slots prefix matched none of these jobs.
* Two review-stage confirmation jobs shipped stdout but **their scripts are not on disk**;
 the manifest marks them `code_on_disk = False` and corroborating rather than load-bearing.
* The manifest's `reproduce` field is a re-run recipe for the *witness encoder*, not for the
 slot jobs, and says so explicitly: "THE SLOT RUNS ARE NOT RE-RUNNABLE AS RECORDED".

We do not restate those runs as reproducible here, and the paper no longer relies on them for
any claim: the census and sweep figures they carried have been removed from the paper, which
now claims no minimality bound and no count of witnesses.

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE. This program re-derives the paper's claims about the one printed graph G = C_28(1,4)
> and nothing else. It examines no other graph, it claims and proves no minimality, and it reads
> no external manuscript -- in particular not the unpublished Brewster-Duchesne manuscript to
> which the source survey credits the vertex-transitive dual that certificate (B) instantiates.

That quoted note and the paper agree: the paper states no minimality bound and no
census/sweep figure, and discusses no second graph.

Three further limits, stated in the paper and repeated here because they are what a referee
should press on:

1. **Attribution, not correctness.** `Brewster–Duchesne (2013)` is unobtainable: no dblp,
 zbMATH or Crossref record, and absent from the deposited references of the 2021 survey
 chapter. We cannot exclude that it already contains a witness of this kind. Two further
 items were unreadable: Hartnell–Mynhardt, *Utilitas Math.* **94** (2014) 19–29, and
 Sen–Kola, *Discrete Appl. Math.* (2024), `doi:10.1016/j.dam.2024.02.010` — the latter one
 of only two zbMATH hits for "broadcast domination circulant", so nothing here claims that
 no circulant broadcast values are in print.
2. **No minimality.** The paper claims no lower bound on the least order of a witness.
3. **Convention.** The paper uses the source's multipacking definition, ranging `s` up to
 `diam(G)`; much of the literature ranges `s` up to `rad(G)`. On the witness
 `rad = diam = 5`, and the program checks the exhibited multipacking under **both** ranges,
 so the result does not turn on the convention.
