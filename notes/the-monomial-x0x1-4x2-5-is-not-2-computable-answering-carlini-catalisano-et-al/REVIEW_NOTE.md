# The Monomial x_0 x_1^4 x_2^5 is Not 2-Computable under Definition 3.5 as Printed

`the-monomial-x0x1-4x2-5-is-not-2-computable-answering-carlini-catalisano-et-al`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

Carlini, Catalisano, Chiantini, Geramita and Woo, *Symmetric tensors: rank, Strassen's
conjecture and e-computability*, Ann. Sc. Norm. Super. Pisa Cl. Sci. (5) **XVIII** (2018)
363-390, ask in print, inside numbered published **Remark 6.2** on journal page 387, whether
`G_1 = x_0 x_1^4 x_2^5` is 2-computable. **The paper answers no under Definition 3.5 as printed**, and the answer is a
three-step argument that fits on half a page:

- **the obstruction (Lemma 2).** For *every* ideal `I` generated in degree `e` and *every*
 `t`, `T_{>= d-e+1}` is contained in `F^perp : I`, because a form of degree `>= d-e+1` times a
 form of degree `e` has degree `> d` and so annihilates `F`. Hence the Hilbert-function sum
 appearing in Definition 3.5 is at most `sum_{i=0}^{d-e} HF(A_F, i) = N(F) - sum_{j<e} HF(A_F, j)`.
 No hypothesis on `t` is used anywhere, so the bound holds under either reading of
 "for general `t` in `I_e`".
- **the application (Theorem 1).** `G_1^perp = (X_0^2, X_1^5, X_2^6)`, `d = 10`,
 `HF(A) = 1,3,5,7,9,10,9,7,5,3,1`, `N(G_1) = 60`; `rk(G_1) = 30` by the monomial Waring-rank
 theorem of Carlini-Catalisano-Geramita, which is also the value the source prints. So
 2-computability would require `60 <= 60 - 1 - 3 = 56`. **Margin 4.**
- **a free remark.** The same cap forces `e <= a_0` for every monomial
 `x_0^{a_0} ... x_n^{a_n}` with `0 < a_0 <= ... <= a_n` whenever that monomial is
 `e`-computable, which narrows the window of the source's Remark 4.3 from unbounded to
 `e <= a_0`. No converse is claimed.

Nothing here is a counterexample and nothing here bears on Strassen's additivity conjecture:
see Scope.

## What was checked, and how

Every quantity in the paper is a small exact count over a monomial complete intersection, and
the proof of Theorem 1 can be redone by hand in a few minutes from the exponent vector
`(1,4,5)` printed in the paper. `verify.py` re-derives all of it mechanically and exactly,
reading only what the paper prints. Beyond the paper's own claims it adds two things a referee
will want:

- **controls that reproduce the source's own positive numbers.** With `I = (X_0^2)`, `t = X_0^2`
 the recipe of Definition 3.5 returns `60`, i.e. `(1/2)*60 = 30 = rk(G_2)` for
 `G_2 = x_0^3 x_1^4 x_2^5`; with `I = (X_0)`, `t = X_0` it returns `30 = rk(G_1)`. Those are
 exactly the two assertions the source makes on page 387, so the machinery is calibrated
 against published values before it is used to refuse one. The cap of Lemma 2 is checked to be
 **silent** in both of those cases (`116 >= 60` and `59 >= 30`), as it must be.
- **an exhaustive census that shows the cap is attained.** Over all 63 nonzero coordinate
 subspaces `W` of `T_2`, `dim_k T/(G_1^perp : (W))` reaches its maximum `56` --- precisely the
 cap --- at `W = span(X_0X_1, X_0X_2, X_1^2, X_1X_2, X_2^2)`, and is `0` at the single subspace
 `span(X_0^2)`, which is the ideal the source's Proposition 4.2 would use and which is why the
 cell was left open. Adjoining `(t)` can only shrink each quotient, so this census bounds the
 quantity of Definition 3.5 for all 63 coordinate ideals without any genericity argument.

One check in an earlier draft of `verify.py` was an over-claim of mine and was corrected rather
than deleted: I had asserted that the cap excludes `e` *exactly* when `e > a_0`. It does not ---
the cap is sometimes strictly sharper. Five (monomial, `e`) pairs with `e <= a_0` are also
excluded, all in two variables, and both the paper (its Remark) and the program now record them.
The paper's remark claims only the implication.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only (`sys`, `fractions`, `itertools`): no third-party
package and no external data file, and the program checks that fact about itself by scanning its
own import lines. Every decision is an exact comparison of Python integers or of
`fractions.Fraction`; no floating-point value is created. It prints one line per check and a
closing verdict, and exits 0 only if every check passes. The recorded run reports **55 checks,
all passing**:

 VERDICT: ALL 55 CHECKS PASS

The recorded exit status is 0, and the run takes about two seconds in a single process. Its
inputs are exactly the objects printed in the paper: the exponent vector `(1,4,5)`, the printed
apolar ideal `(X_0^2, X_1^5, X_2^6)`, the printed Hilbert function, and the generators
`(W^2, YZ, XZ, XY, Y^3-Z^3, X^3-Z^3)` of the apolar ideal of the source's example
`w(x^3+y^3+z^3)`. Among other things
it re-derives: that the printed apolar ideal is the annihilator, on every one of the 560
monomials of degree at most 13; that the apolarity matrix of a monomial is diagonal, so the
annihilator really is a monomial ideal; the Hilbert function by two independent routes (counting
standard monomials, and ranks of the apolarity pairings); Gorenstein symmetry and the socle
degree; Step 1 of Lemma 2 on 2370 explicit products; the cap 56 as an honest `k`-dimension, by
three routes that must agree; the margin 4 and the derived bound `28 < 30`; the census and the
controls described above; the `e <= a_0` remark over 1485 (monomial, `e`) pairs; and, for the source's example
`w(x^3+y^3+z^3)`, the Hilbert function `1,4,6,4,1` computed from `F` itself rather than from the
printed generators, together with a separate confirmation that those generators do generate
`F^perp` in every degree.

Two integers are **quoted, not proved**: `rk(G_1) = 30` (Carlini-Catalisano-Geramita) and
`rk(w(x^3+y^3+z^3)) = 9` (Proposition 7.2 of the source). The program says so on the relevant
check lines and verifies only that the paper uses them correctly.

## Provenance

**Of this folder's program.** `verify.py` was written for this folder and run locally on this
machine. It is *not* one of the scripts that produced the result, and it shares no code with
them. `verify.output.txt` holds its output, preceded by a provenance header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 3798adbcf01382c0225b0930a1b43d324bc6007f3fae03dca6d28d0531f1a90c

`verify.py` reads nothing outside itself, so a referee needs neither the rest of this folder nor
any of the material described next. `paper.pdf` was compiled from `paper.tex` with `tectonic`.

**Of the original computation, stated as its manifest records it and not improved upon.** The
result was reached by hand; the internal run's manifest marks it computer-assisted because two
scripts ran as cross-checks and controls, and it records explicitly that **neither is
load-bearing** --- the claim is the pen-and-paper certificate, which the internal audit stage
re-derived independently from the source, number by number. What that manifest records:

- two scripts, each archived with its SHA-256: a 9,522-byte census over subspaces of `T_2`, and
 a 5,924-byte independent second path computing the same quantity by kernels and ranks in the
 60-dimensional algebra. Both are standard-library-only exact arithmetic over `GF(p)`; no
 Groebner engine, no solver. Randomness is seeded (`20260831`) and used only to draw the
 general `t` in each subspace.
- both jobs are traceable: instance ids `(a fleet slot)` and `(a fleet slot)`,
 command ids `(a dispatch id, redacted)` and
 `(a dispatch id, redacted)`, dispatched 2026-09-01T05:48:49Z and 05:52:03Z, both
 Success with RC=0.
- **honest gaps, recorded in the manifest rather than papered over, and not repaired here.**
 (i) *Neither job's raw standard output was ever captured to a file* --- the dispatch ran
 without the artifact flags, both compute boxes have self-terminated, and no output file exists
 on the control plane. What the record ships instead is an 11,561-byte **transcription** that
 the two attack agents wrote into the row's note at run time: measured integers, not the jobs'
 stdout bytes. (ii) The two planned object-store copies of those outputs are **unread** ---
 neither confirmed present nor confirmed absent --- because the credentials had expired before
 the manifest was written. (iii) Instance types and wall times were recorded for neither job.
 (iv) The invocation line recorded for the census script is a **reconstruction** of the
 dispatch line with the archive path substituted for the scratch path, not a byte-copy of an
 invocation the harness printed; the manifest says so.

Because of (i) and (ii), **no reproduction instruction for the original jobs is offered here.**
The manifest's own `reproduce` field describes what re-running them should print, but the record
does not contain a captured run to compare against, and this note will not present a
reconstruction as a record. Nothing in the paper depends on either job: the census in section 3
of the paper is `verify.py`'s own, is exact over the rationals rather than over `GF(p)`, and
bounds `dim T/(G_1^perp : I)` rather than the slightly smaller `dim T/(G_1^perp : I + (t))` that
the original census measured. The two are consistent (the original reports a maximum of 29,
which is below the 56 established here) but they are not the same quantity, and neither figure
is load-bearing.

## Scope

Theorem 1, together with the corollary and the remark that follow it, is settled under
Definition 3.5 as printed in the source.
What is **not** settled:

- **The definition is load-bearing, and this is the largest risk to the result.** Definition 3.5
 requires `I` to be *generated in degree e*. That clause is what makes Step 1 of Lemma 2 work.
 Restatements in the literature omit it: the definition in the authors' earlier, withdrawn
 preprint arXiv:1502.01107, the restatement in the survey arXiv:1812.10267, and the form in
 Bhat-Carlini-Dubey-Masuti arXiv:2511.23035, where it appears as a hypothesis of their
 Theorem 2.2. **Under those
 readings `I` may contain forms of degree below `e`, Step 1 fails, and Remark 6.2 is not
 answered by anything here.** A referee should check which definition they mean to use before
 reading the theorem.
- **Nothing bears on Strassen's additivity conjecture.** `rk(F + G_1)` is not determined. The
 result closes off the route of Remark 6.2 rather than completing it. The
 value `rk(F + G_2) = 55` quoted in the paper is the source's number, not ours. (An internal
 row key for this result carried the string `strassen-rank-55`; that key oversold it and is
 deliberately not part of this folder's name.)
- **Characteristic zero.** Lemma 2 and Corollary 3 are characteristic-free; `rk(G_1) = 30` is
 not, so Theorem 1 is stated in characteristic 0. An earlier internal draft called the whole
 proof characteristic-free; that claim was withdrawn.
- **The bound is necessary, not sufficient, so it settles nothing where it is silent.** It says
 nothing about the two negative `e = 1` results of the source, and the smallest monomial cell
 it leaves open is `x_0^2 x_1^2 x_2^2` at `e = 2`, where `HF(A) = 1,3,6,7,6,3,1` caps the sum at
 `23` against the required `18`. That cell is genuinely open. There is also **no witness
 object** here of any kind --- the result is a non-existence proof --- so there is nothing to
 encode and no two encodings to cross-check.
- **Precedent is disclosed in the paper and should be read before novelty is assessed.** The
 *shape* of the argument is the source's own: its remark at preprint lines 1387-1424 already
 bounds the same Hilbert-function sum from above in order to defeat 1-computability of
 `w(x^3+y^3+z^3)`, and line 1405 already asserts the vanishing that is our threshold at
 `e = 1`, `d = 4`. What is new is only the uniformity over every `I` generated in degree `e`
 and every `t`, which turns a case check into an obstruction; preprint line 373 says merely
 that the summands vanish "for `s` big enough" and never names `d - e + 1`.
- **Gaps in the prior-art search, recorded because they bound how confidently this can be called
 new.** The arXiv full-text metadata sweep for `"e-computable"` (20 hits) and the zbMATH term
 sweep (29 hits) each found exactly two Waring-sense hits, both by this author family; the 33
 citing papers enumerated by Semantic Scholar contain OpenCitations' 7, and the four
 subject-closest were cleared by grep over their sources. **But:** the OpenAlex `cites:` channel
 is **unread** (rate-limited, then a fetch deadline, twice), so a citer that OpenAlex indexes
 and Semantic Scholar does not would be invisible; zbMATH's deposited reference list for the
 record is empty, so that channel was unavailable rather than clean; Teitler arXiv:1604.07691
 was checked at identifier level only, never at source; MathSciNet was not consulted at all;
 and no erratum or correspondence channel was checked, which matters here because both halves
 of the target are author hedges and the likeliest way this is not new is that the same authors
 settled it somewhere unlooked. Duplication against other rows of the same internal batch was
 checked by no stage, with no evidence either way.

The program's own closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE -- what this program does NOT cover. rk(G_1) = 30 and rk(w(x^3+y^3+z^3)) = 9 are QUOTED theorems ([CCG] and Proposition 7.2 of the source); this program checks that the paper uses them correctly, and does not reprove them. Corollary 4 is executed here only as a FINITE CENSUS: 120 monomials in 2-4 variables with exponents 1..5, 1485 (monomial, e) pairs. Arbitrary n and arbitrary 0 < a_0 <= ... <= a_n are NOT executed anywhere in this program; that generality is the paper's proof of Corollary 4, imported. Proposition 4.2 of the source, which supplies e = 1 in the a_0 = 1 slice, is likewise quoted, not reproved. the obstruction is a NECESSARY condition, so it settles no cell on which it is silent: x_0^2x_1^2x_2^2 at e = 2 (cap 23 >= 18) is the smallest such monomial cell and remains OPEN. the e = 1 half of Corollary 5 is NOT reproved here: it is Remark 7.3 of the source, quoted. For the other half this program supplies the finite range e = 2..11 only; it does not reduce the unbounded range e >= 2 to finitely many cases, which is the paper's step, imported. NOT RE-RUN: the original run's GF(p) census of dim T/(G_1^perp:I+(t)) over 1,263 subspaces (whose stdout was never captured to a file); the 63-subspace census printed above is this program's own, is exact over Q, and bounds dim T/(G_1^perp:I) rather than dim T/(G_1^perp:I+(t)). ALSO NOT RE-RUN: anything about Strassen's conjecture, and anything under the looser readings of e-computability that omit the clause "generated in degree e" -- under those readings step 1 of Lemma 2 fails and Remark 6.2 is not settled.

## One thing a referee should check against the source

What is answered is an open question stated inside numbered, published **Remark 6.2** --- not a
numbered Conjecture, Problem or Question, and it must not be cited as one. Its locator: journal
page 387, which is page 25 of the 28-page article file; in the preprint arXiv:1506.03176v2, a
single inner file `2015-06-11-CCCGW.tex` of 65,505 bytes and 1,517 lines, the enclosing `rem`
environment spans lines 1293-1307, the passage the paper quotes is lines 1300-1306 (the remark's
last three paragraphs), and the two sentences about `G_1` are lines 1300-1301. The paper
reproduces that passage verbatim, so what is answered does not depend on the labelling.
All statement numbers used (Definition 3.5, Corollary 3.4, Proposition 4.2, Remark 4.3,
Remark 6.2, Proposition 7.2, Remark 7.3) are those of the version of record and were each read
there; Corollary 3.4 prints correctly in the journal, while the preprint has a typo in it, which
is why the journal text is cited.

The preprint-side locators above were re-checked at bundle-review time against the e-print source
fetched fresh from `arxiv.org/e-print/1506.03176v2`: it is the single file
`2015-06-11-CCCGW.tex`, 65,505 bytes and 1,517 lines exactly as stated, the quoted passage is
byte-for-byte the text at lines 1300-1306 (with the source's `\ref{monomi}` resolved to
"Proposition 4.2" and its display `$$\rk(F+G_2)=25+30=55.$$` inlined), the definition's wording is
at lines 379-381, the two precedent claims are at lines 373 and 1405 as stated, and the preprint's
Corollary 3.4 does carry the typo the paper describes. Two locator defects were found and fixed in
that pass and are listed under Corrections below. What was **not** re-checked is everything on the
journal side: page 387, page 25 of 28, and the statement numbers themselves rest on the earlier
stages' reading of the version of record, not on this pass.

## Corrections made at bundle review

- The preprint line range for Remark 6.2 was printed as **1291-1305** in both `paper.tex` and this
 note. The `rem` environment actually spans **1293-1307**; the stated range was shifted by two and
 did not even contain line 1306, which carries the `rk(F+G_2)=55` display that the paper quotes.
 Corrected in both files, and the quoted passage's own range (1300-1306) is now given.
- The definition environment in `paper.tex` was titled "Definition 3.5 of the source, **verbatim**"
 while its body was the source's definition restated in this paper's notation (the source also
 fixes `e>0`, writes the sum as `\sum_{i=0}^\infty`, orders the "for general t" clause before the
 display, and continues with two sentences about "the rank of F is computed by I and t"). The
 title now claims only "Definition 3.5 of the source", and the source's own wording of the
 load-bearing clause is quoted in the following sentence with its line locator. No mathematical
 content changed: the restatement was faithful, only the word "verbatim" was wrong.
- References [3] (Carlini-Catalisano-Oneto) and [4] (Teitler) appeared in the bibliography with no
 citation anywhere in the text. They are now cited once, in the Statement-numbering paragraph, as
 adjacent additivity literature that no argument here uses.
- The preprint's Corollary 3.4 typo, previously asserted without content, is now stated: the
 preprint's display reads `HF(T/(F^perp:(t)+(t)),i)` where `HF(T/(F^perp:I+(t)),i)` is meant.
