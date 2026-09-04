# Review note

**Paper.** *Non-Monotonicity in k at n = 3: A Census Contradicting a k-Monotonicity Clause
Transcribed from D'Angeli and Rodaro* (8 pages, Sections 1–6).

**Files in this folder.** `paper.tex`, `paper.pdf`, `verify.py`, `verify.output.txt`. Nothing
outside the folder is needed to check the paper, except the cited e-print itself, which the
paper says it cannot supply (see §3 below).

## 1. What the paper claims

For a uniformly random *k*-out digraph on *n* labelled vertices in the ordered
(with-replacement) model fixed in Section 2 — probability space `V^{kn}` with the uniform
measure, states weighted by the multinomial `W(c)` of display (2) — put
`P(n,k) = Pr(G totally simple | G strongly connected) = A(n,k)/B(n,k)`, displays (4) and (5).
**Theorem 1** asserts

    A(3,2) = 114     B(3,2) = 296       P(3,2) = 57/148    = 0.385135...
    A(3,3) = 5132    B(3,3) = 13754     P(3,3) = 2566/6877 = 0.373127...

hence `P(3,2) > P(3,3)`, decided by the integer inequality
`114 * 13754 = 1 567 956 > 1 519 072 = 5132 * 296` (display (6)); the gap is `12221/1017796`,
about 1.2 percentage points (display (7)).

What that settles is display **(1)** of Section 1: `P(n,k) <= P(n,k+1)` for every `n >= 1` and
every `k >= 2`. This is the second of three clauses of an unnumbered conjecture that Section 1
transcribes **by hand** from tex line 1463 of the LaTeX source `totally_main.tex` of the
e-print cited as [1]; the paper reports the clause is unnumbered there, carries no `\cite`, and
is attributed to nobody. Both `k = 2` and `k = 3` lie inside the transcribed clause's own
"Fix `k >= 2`"; `k = 1` is explicitly not used (`P(3,1) = 1 > P(3,2)` is also a decrease but is
outside the source's range).

The proof is a hand proof and needs no program. **Proposition 2** reduces total simplicity at
`n = 3` to the three column inequalities `a_10 != a_20`, `a_01 != a_21`, `a_02 != a_12`;
**Proposition 3** gives `B(3,k) = (3^k-1)^3 - 3(2^k-1)^2(3^k-1)`; **Proposition 4** gives
`A(3,k) = tr((W_k D)^3) - 3[(sum f)^2 - sum f^2]`; Step 4 of Section 3 evaluates both at
`k = 2, 3, 4` from a 3×3 and a 4×4 integer matrix. **Proposition 5** (Section 4) is secondary
but load-bearing for the numerator: the 2-out graph `0 -> {0,0}`, `1 -> {0,2}`, `2 -> {1,1}` is
totally simple and not strongly connected, so `P(n,k)` is not `Pr(TS)/Pr(SCC)`, and two further
sentences the paper transcribes from [1] (tex line 1455, and a remark at tex 1516–1521) are
false.

## 2. What the program checks

`verify.output.txt` records **`VERDICT: ALL 56 CHECKS PASS`** and `program exited with status 0`,
in eight labelled blocks. Block by block, with the claim each supports:

| block | checks | supports |
|---|---|---|
| Step 1 — the predicate: the partition lattice, Δ/∇ excluded | 5 | display (3) and Proposition 2. The column criterion agrees with a generic decider reading lumpability from the definition over all set partitions (4591 states at `k = 2,3,4`, 0 mismatches); Δ and ∇ are lumpable on all 216 states at `n=3, k=2`, so excluding both is forced |
| Step 2 — exhaustive census in the ordered model, `n = 3` | 12 | Theorem 1 and Step 4: `B(3,2)=296`, `A(3,2)=114`, `B(3,3)=13754`, `A(3,3)=5132`, `B(3,4)=458000`, `A(3,4)=180890`, with the total-weight identity `n^{kn}` and a state-count identity in each cell |
| Step 3 — a second enumeration over raw ordered edge tuples, weight 1 | 4 | display (2): 729 and 19683 raw tuples at `n=3, k=2` and `k=3` reproduce the same `A` and `B` with no multinomial weight anywhere |
| Step 4 — the closed forms against the census | 7 | Propositions 3 and 4 and Step 4: traces 150, 5798, 189518; corrections 36, 666, 8628; and two independent closed-form routes (transfer matrix, inclusion–exclusion) agree with each other over the 29 values `2 <= k <= 30` |
| Step 5 — the exact comparison | 11 | Theorem 1 and displays (6), (7), (9): `57/148`, `2566/6877`, `18089/45800`, the cross-multiplication, the drop `12221/1017796`, `P(3,4) > P(3,3)`; and over `2 <= k <= 30` the violating `k`-steps are exactly `[2]`, with `P(3,k)` strictly increasing from `P(3,3)` to `P(3,30) = 0.708531594` |
| Step 6 — the `n = 4`, `k = 2` census | 5 | Section 5: its own census yields `B(4,2) = 20958` and `B(3,2) = 296`, matching the two hand-transcribed `A027834` terms; also `A(4,2) = 6720` and `P(4,2) = 160/499` |
| Step 7 — the model-sensitivity disclosure, itself checked | 6 | the Remark of Section 2: rival-model values `34/65`, `133/243`, `87/148` with no violation at `n=3` for `k = 2,3,4`; and `Pr(no edge lands on v) = 64/729` in the ordered model against `1/8` in the rival one, the discriminating identity the paper quotes from tex 1452 |
| Step 8 — exhibited states and the TS-but-not-SCC weights | 6 | Proposition 5 and the closing sentence of Section 4: an exhibited totally simple state that is not strongly connected; TS-and-not-SCC weight 36 at `k=2`, 666 at `k=3`, 8628 at `n=3,k=4`, 1584 at `n=4,k=2`; plus controls not printed in the paper (the graph `0->{1,1}, 1->{1,2}, 2->{2,0}`, a state that is strongly connected but not totally simple, and `n = 2` where `P(2,2) = P(2,3) = 1`) |

Every decision is exact integer or `Fraction` arithmetic; the run records that the decrease is
decided without any floating point.

## 3. What the program does *not* check

The transcript states its own limits in four closing `NOTE SCOPE` lines; they are carried over
here, with the paper's matching disclaimers.

* **The main theorem is a hand proof and the program is a control.** Theorem 1 and
  Propositions 2–5 are proved in closed form in the paper, which says a reader needs no
  program; `verify.py` re-derives the printed integers independently by enumeration, but it does
  not verify any proof.
* **It does not read or parse the source.** The quoted clause, all tex line numbers, the
  definitions of Section 2 and the identification of the probability model are hand
  transcription. Section 1 states that the paper can supply neither an immutable copy of that
  e-print nor a cryptographic hash of it, and has not confirmed that the arXiv identifier
  resolves for a third party. That transcription is unchecked here in full and is what a
  referee must check against [1] directly. The transcript's banner nevertheless speaks of a
  counterexample and calls the clause refuted; Section 6 flags that this is the reading of the
  source, not something the program establishes.
* **The model is a reading, and it is load-bearing.** Both readings are computed. The
  counterexample exists in the ordered/with-replacement model only and vanishes under the rival
  uniform measure on labelled *k*-out multigraphs, where `P_u(3,2) = 34/65 < 133/243 = P_u(3,3)`.
  Since the drop is only ~1.2 points, the residual risk is this reading, not the arithmetic.
* **The asymptotic clauses are untouched.** `P(n,k) -> 1` as `n -> infinity` (tex 1461) and the
  expansion `1 - exp(-c_k n + beta_k) + o(1)` with `c_k` increasing (tex 1465–1470) are not
  addressed; no finite computation can settle them and the paper claims nothing about them.
  They remain open.
* **Cells not enumerated.** Exhaustively enumerated: `n=3` at `k=1,2,3,4` (and `n=3` at
  `k=2,3` a second time over raw ordered tuples), `n=4` at `k=2`, `n=2` at `k=2,3`. The closed
  forms are exercised at `n=3` only, for `k = 2..30`. Cells with `n >= 5`, and `n = 4` with
  `k >= 3`, are not recomputed, so the repair "for fixed `n >= 4`" floated in Section 5 is, in
  the run's words, **consistent with, not verified by**, this program. Section 3 likewise
  disclaims anything about `k > 30`, and Section 5 disclaims that `k = 2 -> 3` is the only
  violation along `n = 3` beyond `k = 30`.
* **Quantities transcribed rather than recomputed.** `296` and `20958` sit in the program's
  `CLAIMED` block as hand-typed literals and are only *compared* against its own census; no
  database is consulted. The `A006691`/`A006692` normalisation by which Section 5 places
  `B(3,3) = 13754` on record is not touched, and Section 5 says the paper has not verified that
  normalisation independently. The sampled table at `n = 8, 9, 10` attributed to [1]
  (tex 1484–1506) is not recomputed, nor is any other attribution in Section 5.

## 4. How to check it

```sh
python3 verify.py            # one PASS line per check; exits 0 iff all pass
shasum -a 256 verify.py
```

Python 3.9 or later, standard library only (`itertools`, `sys`, `fractions`, `math`): no
third-party package, no data file, no network. `verify.output.txt` opens with a provenance
header giving the program name, its SHA-256 and the interpreter version (Python 3.9.25), so the
transcript can be paired with the program shipped beside it. The digest, computed here from the
shipped file, is

    75888c5fc3348f349ffd9104d7cdc87f436a589f574022292fee55c29ba47d24  verify.py

which is the value that header carries.

## 5. One typographic note

The Remark of Section 2 is printed unnumbered, headed *Remark (the model is load-bearing)*, but
the cross-references to it on pages 7 and 8 render as "Remark 2"; there is no separately
numbered Remark 2 in the paper, and "2" there is the section number. Mathematically nothing
turns on it.
