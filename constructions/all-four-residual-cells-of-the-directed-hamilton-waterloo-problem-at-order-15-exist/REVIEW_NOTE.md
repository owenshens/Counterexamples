# Review note

Paper: *All Four Residual Cells of the Directed Hamilton--Waterloo Problem at Order 15 Exist*

Files a referee has here, and nothing outside this list is referred to below:

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper: 4 sections, Theorem 1, Lemma 2 |
| `verify.py` | the checking program |
| `verify.output.txt` | the recorded run of that program |
| `REVIEW_NOTE.md` | this note |

## 1. What the paper claims

HWP\*(v; m^r, n^s) asks for a decomposition of the arc set of the complete symmetric digraph
K\*_v into r spanning factors that are disjoint unions of directed m-cycles and s spanning
factors that are disjoint unions of directed n-cycles; Section 1 notes that counting arcs and
vertices forces r+s = v-1, m | v and n | v, and that K\*_15 has 15·14 = 210 arcs. **Theorem 1**
claims that HWP\*(15; 3^11, 5^3), HWP\*(15; 3^12, 5^2), HWP\*(15; 3^13, 5^1) and
HWP\*(15; 3^13, 15^1) all have solutions — a witness for each is printed in Section 3 — and hence
that the exception clause "except possibly for r ∈ {11,12,13} when (m,n) = (3,5) and for r = 13
when (m,n) = (3,15)" may be deleted from printed Lemma 4.2 of Yetgin, Odabaşı and Özkan (*On the
directed Hamilton--Waterloo problem with two cycle sizes*, Contrib. Discrete Math. **20** (2025),
no. 1, 74--94; e-print arXiv:2209.14588).

So the statement settled is that Lemma 4.2, which Section 1 quotes verbatim. Section 1 is explicit
about the logical situation: with r+s = 14 the clause names exactly the four cells
(m,n,r,s) = (3,5,11,3), (3,5,12,2), (3,5,13,1), (3,15,13,1); each satisfies the necessary
conditions, so what is undetermined in each is sufficiency; and the source "poses no problem about
them", it leaves those four values undetermined.

The evidence is on the page. **Lemma 2** (Section 2, proved there) says that for m,n | 15 and
r+s = 14 a solution is the same thing as a list σ_1,…,σ_14 of fixed-point-free permutations of
**Z**_15 with (i) each σ_i of a single uniform cycle length, that length being m for r of the
indices and n for the remaining s, and (ii) the 210 arcs (x, σ_i(x)) pairwise distinct;
equivalently, the array whose row 0 is the identity and whose row i is σ_i is a Latin square of
order 15. Section 2 reduces the check on each witness to three lines (a), (b), (c); Section 3
prints the four witnesses in that form, as 14 labelled permutations of **Z**_15 in cycle notation;
and the proof of Theorem 1 is that hand check against the printed cycles. The abstract states that
a reader can check it by hand and needs no code.

## 2. What the program checks

`verify.output.txt` records `verify.py` printing one line per check and closing with

    VERDICT: ALL 53 CHECKS PASS
    === program exited with status 0 ===

The 53 checks fall into four blocks.

* **"Step 0: the ground set and the exception clause" — 8 checks.** K\*_15 has 210 arcs; a
  solution needs exactly 14 spanning factors (210/15 = 14 = v-1); the listing parses into four
  labelled cell headers; the four cells named by the exception clause are exactly those four, the
  run printing the re-derived list `[(3, 5, 11, 3), (3, 5, 12, 2), (3, 5, 13, 1), (3, 15, 13, 1)]`
  under the note "the clause is re-derived from its own wording, not copied from the listing";
  then one check per cell for the necessary conditions (m | 15, n | 15, r+s = 14, 15·14 = 210
  arcs). This block corresponds to the identification of the four cells in Section 1.
* **One block per cell — 9 checks each, 36 in all.** Per cell: 14 factors printed against 14
  required; every factor a permutation of **Z**_15; every factor fixed-point-free ("no cycle of
  length 1 and no loop"); every cycle type equal to its own printed label; the census of labels
  equal to the cell's (e.g. "labels seen {3: 13, 15: 1}, cell wants {3: 13, 15: 1}"); the 210 arcs
  exactly the arcs of K\*_15 *as a set identity*, none missing and none repeated; the
  identity-adjoined array a Latin square, all 15 columns permutations of **Z**_15; i ↦ σ_i(0) a
  bijection onto {1,…,14}; and a summary line that every checked condition of the criterion holds.
  These four blocks are the mechanical counterpart of Theorem 1 via conditions (i) and (ii) of
  Lemma 2.
* **"Step 2: the four witnesses are distinct objects" — 1 check.** The four are pairwise distinct
  sets of 14 permutations. The paper claims no such distinctness; this check is extra to it.
* **"Step 3: controls, both polarities" — 8 checks.** Two positive, on an object independent of the
  witnesses: the 14 translations x ↦ x+c of **Z**_15 partition the arcs of K\*_15, with the
  predicted cycle types (3^5 for c ∈ {5,10}, 5^3 for c ∈ {3,6,9,12}, 15^1 when gcd(c,15) = 1). Six
  negative, over five deliberate damages to the (3^11, 5^3) witness — a factor reversed, a factor
  duplicated, a factor deleted, a cycle type mislabelled, a printed opening parenthesis dropped —
  each rejected by the same checker with the reason printed; the reversal is convicted twice, by
  the arc-set identity and again by the σ_i(0) diagnostic.

## 3. What the program does not check

* **Theorem 1 is a hand proof and the program is a control.** The proof of Theorem 1 verifies (a),
  (b), (c) on the printed cycles and appeals to Lemma 2; that proof itself says `verify.py` "redoes
  the same check mechanically on the same printed text". It is not an independent derivation of the
  result, and Lemma 2 — the equivalence the whole check rests on — is proved by hand only.
* **Its only input is a listing carried in its own source.** The program's header asserts that this
  listing is character-for-character the one typeset in Section 3 of `paper.tex`; that identity is
  not one of the 53 checks. Both texts are in this folder, and a referee who wants the point tied
  down should compare them rather than take the header's word.
* **No non-existence, census or exhaustion claim is checked or made.** The run's closing note says
  so in those words and adds "every verdict above is the existence of an exhibited object";
  Section 4 says the same. Nothing is claimed about uniqueness or about how many solutions a cell
  admits.
* **The search that produced the four objects is not re-run and is not part of this folder.**
  Section 4 states that the witnesses were found by a constraint solver, that nothing in Theorem 1
  depends on the solver, that a re-run would return different objects, and that the printed
  objects are the claim. The run's closing NOT RE-RUN note carries the same, adding that the
  searches' logs were not preserved. No claim in the paper rests on that search.
* **The other cases of printed Lemma 4.2 are not re-verified.** The NOT RE-RUN note names "the
  cases of printed Lemma 4.2 outside the four cells, which are the source authors' own". Section 4
  identifies what it combines the four witnesses with: all of (m,n) = (5,15), and every (3,5) and
  (3,15) case outside r ∈ {11,12,13} and r = 13 respectively; those "were not re-verified here".
  The deletion of the exception clause therefore rests on the source's cases together with the
  four witnesses printed here. Section 4 also records that Lemma 4.1 of the source, which lifts a
  complete list of solutions at one odd order to larger odd orders, is the authors' own and is
  neither verified nor used.
* **Nothing outside v = 15.** The run's closing SCOPE note says the program checks the four
  exhibited factorizations of K\*_15 "and nothing else"; Section 4 says nothing outside v = 15 is
  constructed or claimed.
* **The quotation of printed Lemma 4.2 is transcribed, not derived.** Section 1 states that it is
  taken from the arXiv e-print (arXiv:2209.14588v1), that the statement stands unchanged in the
  version of record, and that only the HWP\* macro and the line breaks are the paper's. No check
  compares the quotation with the source, and the source is not in this folder.

## 4. How to check it

```sh
python3 verify.py
shasum -a 256 verify.py
```

`verify.py` prints one line per check and, per its header, exits 0 if and only if every check
passed; the recorded run exited 0 with all 53 checks passing. Its header states that the program
uses the standard library only (Python 3.9+), with no external file, no third-party module, no
solver and no floating point. The first lines of `verify.output.txt` carry, beside the interpreter
version used for that run (Python 3.9.25), the SHA-256 of the program that produced it, so
transcript and program can be paired. Computed here from the shipped `verify.py`:

    b7cd3e61f087d06527b978575eb40e48241bfff85e09dcf92545f9c9748fbf48

which is the digest printed in that header. Independently of the program, the four witnesses can be
checked by hand from Section 3 against conditions (a), (b), (c) of Section 2.
