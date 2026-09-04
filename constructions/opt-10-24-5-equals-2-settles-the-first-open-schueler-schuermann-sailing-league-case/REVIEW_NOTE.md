# Referee note: opt(10,24,5) = 2

Files in this folder, and nothing else is needed to check the paper: `paper.tex` and `paper.pdf` (the
paper), `verify.py` (the verification program), `verify.output.txt` (the recorded run of it).

## 1. What is claimed

Schüler and Schürmann (*Sailing league problems*, J. Combin. Des. **32** (2024), no. 4, 171–189,
doi:10.1002/jcd.21929; arXiv:2212.02865 — reference [6] of the paper) close their case study of the
Asian Pacific Sailing Champions League with unnumbered prose recording that the parameters
N_teams = 10, N_inrace = 5, N_flights = 8n with n ≥ 3 "still appears to be open". Section 2 quotes
that sentence verbatim, with the source's definitions of a pairing list and of the sailing league
problem, and the paper settles exactly its n = 3 member, N_flights = 24:

> **Theorem 1.** opt(10,24,5) = 2.

The upper bound is one explicit 24 × 10 tournament plan over {1,2}, printed in Section 3 twice — as
the grid of Table 1 and as the same 24 rows in a character-exact row-string list — whose 45 pair
counts are 30 tens and 15 twelves, so λ_min = 10, λ_max = 12 and the utility λ_max − λ_min is 2. The
Remark of Section 3 discloses that the plan uses 21 distinct flights, flight 7 = flight 10, 8 = 11,
9 = 12, which the source's definition permits ("it is allowed to repeat races and even flights").

The matching lower bound is Theorem 5 in Section 4, which opens "No computer is needed for any step
of this section." Lemma 2, the authors' own averaging lemma at this cell, gives λ = 96/9 = 32/3 ∉ ℤ
and kills utility 0. Lemma 3 (triangle parity) states λ(t,u) + λ(t,v) + λ(u,v) ≡ F (mod 2) for any
(10,F,5) pairing list, so at F = 24 every triangle sum is even. Lemma 4 shows that a graph on the ten
teams all of whose 120 triples span an even number of edges is a complete bipartite cut, of edge count
|U|(10−|U|) ∈ {0,9,16,21,24,25}. Utility 1 would force the window {10,11} and a 6-regular graph on the
pairs at 11, which by Lemmas 3 and 4 would have to be a cut K_{k,10−k} with 10 − k = 6 and k = 6 at
once.

## 2. What the program checks

`verify.output.txt` closes with `VERDICT: ALL 61 CHECKS PASS` and
`=== program exited with status 0 ===`. Its header records the program name, its SHA-256 and
`python: Python 3.9.25`; the arithmetic is exact integers and `fractions.Fraction` only, with no
third-party package. The 61 checks fall into twelve labelled steps.

Blocks matched to the paper's claims:

* **Upper bound — Steps 1–2, 13 checks.** Table 1 and the row-string list agree on all 24 of 24 rows;
  every row is a five-five split; all 45 entries of the printed λ table are reproduced; λ_min = 10,
  λ_max = 12, utility 2, profile [(10, 30), (12, 15)], no pair at 11; grand total 480 over the 45
  pairs; every per-team sum 96.
* **The Remark of Section 3 — Step 4, 4 checks.** 21 distinct rows of 24, repeat classes exactly
  (7,10), (8,11), (9,12); rows 1–12 differ from rows 13–24 and rows 1–8, 9–16, 17–24 are not all
  equal, so the plan is neither a doubled 12-flight plan nor three copies of one base schedule.
* **Lemma 2 — Step 5, 3 checks.** λ = 24·4/9 = 32/3, denominator 3, so no perfect pairing list at 24
  flights; per-team sums and grand total agree with the average (960/90 = 32/3).
* **Lemma 3 — Step 6, 3 checks.** All 126 five-five splits against all 120 triples: 15120 checks, 0
  violations; and all 120 triangle sums even on the exhibited plan at F = 24.
* **Lemma 4 — Step 7, 6 checks.** Exactly 512 graphs on ten labelled vertices have all 120 triples
  even; every one is a complete bipartite cut, each bipartition exhibited and verified against all 45
  pairs; realised edge counts exactly 0, 9, 16, 21, 24, 25; realised degree sequences exactly 0^10,
  1^9 9, 2^8 8^2, 3^7 7^3, 4^6 6^4, 5^10; none 6-regular; the only regular degrees realised are 0
  and 5.
* **Theorem 5, and its combination with Table 1 — Step 8, 5 checks.** 9m ≤ 96 ≤ 9m+9 has the unique
  solution m = 10; 90 + a(t) = 96 forces a(t) = 6 for every team; 30 edges is not of the form k(10−k)
  and 6 is not a realised regular degree, so utility 1 is impossible; hence the lower bound, and with
  Step 2, opt(10,24,5) = 2. (The transcript heads this block "Step 8: Theorem 1" and names its check
  `theorem_1_opt_10_24_5_is_at_least_2`; the statement re-played is the paper's Theorem 5, which with
  Table 1 gives the paper's Theorem 1.)

Four further blocks are controls and by-products beyond what the paper asserts:

* **Step 3, 5 checks.** The fifteen λ = 12 pairs computed from the plan are exactly the printed ones;
  that graph is cubic and strongly regular with parameters (10,3,0,1) and is identified with the
  Kneser graph K(5,2) by an explicit bijection; design parameters v = 10, k = 5, r = 24, b = 48, with
  λ_1 = 10 on intersecting 2-subsets and λ_2 = 12 on disjoint ones. Of these the paper claims only
  v = 10, k = 5, r = 24 and the concurrences 10 and 12, in the "Adjacent work" paragraph of Section 5.
* **Step 9, 9 checks.** Census of arithmetically admissible utility-2 profiles at 24 flights: the
  window {9,10,11} is impossible; 211 odd graphs survive, each pinning its per-team counts uniquely;
  exactly two profiles remain up to relabelling, P1 = ((10,30),(12,15)) and
  P2 = ((10,18),(11,24),(12,3)), both closing at 480; the surviving odd graphs are the empty graph
  (P1) and the 210 cuts K_{4,6} (P2); the exhibited plan realises P1; and k = 5 (odd graph K_{5,5})
  is excluded by d(t) = 6 − 2c(t).
* **Step 10, 4 checks.** The same route gives opt(10,12,5) ≥ 2, opt(10,16,5) ≥ 2, opt(10,32,5) ≥ 2,
  plus an 8-flight control confirming the route never contradicts the published opt(10,8,5) = 3.
* **Steps 11–12, 3 + 6 checks.** The source's own published 16-flight table is a legal pairing list
  and reproduces the 16-flight theorem (λ_min = 6, λ_max = 8, utility 2, profile [(6,20),(8,25)],
  total 320, per-team sums 64) and, on its first 8 rows, the 8-flight theorem (λ_min = 2, λ_max = 5,
  utility 3). Truncations of it to 7, 8, 9, 11, 13 flights carry 20, 24, 36, 20, 20 odd-λ pairs, so
  the parity prediction is falsifiable there rather than vacuous; in each the odd-λ graph has the
  predicted shape, a cut for even F and the complement of a cut for odd F.

## 3. What the program does **not** check

* **The lower bound is a hand proof and the program is a control.** Section 4 needs no computer. The
  program re-plays the arithmetic of Theorem 5 and exhausts the two finite censuses behind Lemmas 3
  and 4; it does not supply the proof. Section 6 states that neither the program nor the transcript
  is needed in order to redo the work, and that a reader who trusts no program can verify Theorem 1
  by hand.
* **The upper bound is a checked object, not a reproduced search.** Section 6 records that the plan
  of Table 1 was found by a constraint solver, that this is provenance only, that the captured output
  of the search is not part of this record and that no program re-runs it. The transcript repeats it:
  "NOT RE-RUN: the search that FOUND the plan ... Nothing above depends on the search." That search
  is not in this folder and no claim of the paper rests on it.
* **The quantifier is not closed.** The source sentence ranges over all integers n ≥ 3; only n = 3 is
  settled. Per the transcript's SCOPE note, n = 4 (32 flights), n = 5 (40 flights) and every larger n
  remain open, only the bound opt(10,32,5) ≥ 2 is established, and no 32-flight witness exists in
  this record; Section 5 says the same.
* **Lemma 3 holds for all F, and only its per-flight step is exhausted.** The 15120 checks cover every
  five-five split against every triple, i.e. the single-flight contribution; the accumulation over F
  flights is the hand argument, and the F = 24 statement is confirmed only on the exhibited plan, plus
  the truncations of Step 12 at F = 7, 8, 9, 11, 13, 16. No other plan is tested. The paper notes
  that Lemma 3 uses that each flight has exactly two races, so it is claimed for none of the 18-team
  or 32-team case studies of the source.
* **Bounds, not values, at the other cells.** The transcript states that whether any (10,24,5) plan
  realises P2 is unknown and irrelevant to opt(10,24,5) = 2, and that the exact values of
  opt(10,12,5) and opt(10,16,5) are not claimed — only ≥ 2 in each case, opt(10,16,5) = 2 being the
  authors' own theorem, quoted and not reproved. Nothing in Steps 3 and 9 says every utility-2 plan
  at 24 flights has a cubic λ = 12 graph, and the paper claims no such thing.
* **The distinct-flights variant is untouched.** The Remark of Section 3 and the transcript's second
  SCOPE note both say that question, and its optimum, are answered in neither direction.
* **One object used by the controls is not printed in the paper.** The source's 16-flight table
  driving Steps 11–12 is transcribed inside `verify.py` (constant `PAPER16`) from the cited e-print's
  LaTeX source; that transcription cannot be audited from this folder alone, and no claim of the paper
  depends on it.
* **Nothing bibliographic or prior-art is computed.** Per the transcript's last closing note, the
  journal reference, the DOI, the citer registries and the design catalogues are reads, not
  computations, and are not checked. That covers the paper's quotations of the source, its quoted line
  numbers, its statement that `Nflights = 24` does not occur in the source's LaTeX file, and the two
  quoted theorems opt(10,8,5) = 3 and opt(10,16,5) = 2. Section 5 records its own limits: two
  catalogues not consulted page by page because both are paywalled, two search services not
  reachable, citer registries reporting no citing work, and no priority claim made.

## 4. How to check it

By hand, from the paper alone: check that each of the 24 rows of Table 1 carries five 1s and five 2s,
count same-symbol columns for each of the 45 pairs against the printed λ table, and read Section 4.

With the program:

```sh
shasum -a 256 verify.py
python3 verify.py
```

Python 3.9 or later, standard library only, no external data file and no argument. It prints one
`PASS <name> [detail]` line per check and exits 0 iff every check passed. The header of
`verify.output.txt` carries the SHA-256 of the program beside its output, so the two files can be
paired; recomputed from the shipped `verify.py` it is

    847f82d47a34b2757a5e032de95f04027b7e6d4c417b4fd8416e2d54932d7582

which is the value in that header.
