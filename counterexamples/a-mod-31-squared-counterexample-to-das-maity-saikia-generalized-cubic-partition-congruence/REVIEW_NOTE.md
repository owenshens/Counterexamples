# Referee's note: a mod $31^2$ counterexample to a generalized cubic partition congruence of Das, Maity and Saikia

Files in this folder, and nothing else is needed:

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper (5 sections, one theorem, one remark) |
| `verify.py` | the verification program (Python 3.9+, standard library only, exact integers) |
| `verify.output.txt` | the recorded run of `verify.py`, whose header carries the SHA-256 of the program |

## 1. What the paper claims

For $c\ge1$, $a_c(n)$ is defined by equation (1) of the paper,
$\sum_{n\ge0}a_c(n)q^n=1/(f_1f_2^{\,c-1})$ with $f_k=(q^k;q^k)_\infty$ (even parts in $c$ colours).
Section 1 quotes, from the Concluding Remarks of arXiv:2503.19399v2, a conjecture of Das, Maity and
Saikia consisting of three displayed congruences to a prime-square modulus, numbered (2), (3), (4) in
the paper; the first is
$$a_{25}(31^2n+644)\equiv0\pmod{31^2}\qquad\text{for all }n\ge0 .$$
**Theorem 2 asserts that (2) is false, and fails at $n=2$:** with $N=31^2\cdot2+644=2566$,
$a_{25}(2566)\equiv682=22\cdot31\pmod{31^2}$, so $v_{31}(a_{25}(2566))=1$ exactly. The proof prints
$a_{25}(2566)$ as a 182-digit integer $W$ and $Q=W/31$ as a 180-digit integer, and notes
$Q\equiv22\not\equiv0\pmod{31}$. So the mathematical content is one coefficient of a 2567-term
truncation of $1/\bigl((q;q)_\infty(q^2;q^2)_\infty^{24}\bigr)$ plus one division by 31.

The paper identifies the refuted line by its position in the source's Concluding Remarks and says in
Section 1 that it makes no claim about how the source numbers the conjecture (noting that another
preprint calls the *third* line "Conjecture 7.3"). Section 3 records that the first-power statement
(5), $a_{25}(31^2n+644)\equiv0\pmod{31}$, is the $p=31$ case of Theorem 3 of Guadalupe, *Integers*
**25** (2025), Article A20, and is quoted, not proved; only the lift to $31^2$ is refuted. Remark 1
notes that the source's *prose* ("each part may appear in $c\ge1$ different colors") would give
$1/f_1^{\,c}$, a different object, under which the congruence at issue already fails at $n=0$ with
residue 550 mod 961; equation (1) is taken as authoritative.

## 2. What the program checks

`verify.output.txt` ends with `VERDICT: ALL 35 CHECKS PASS` and `program exited with status 0`. The
35 checks are grouped into the transcript's seven sections:

* **Section 1, 5 checks** — Theorem 2's arithmetic on the printed digits alone: 182 and 180 digits,
  $31Q=W$, $W\bmod961=682$, $W\bmod31=0$, $Q\bmod31=22$, $31^2\cdot2+644=2566$, and
  $682=22\cdot31\ne0\bmod961$.
* **Section 2, 7 checks** — the one computational step of the proof, the identity $a_{25}(2566)=W$.
  $a_{25}(0..2566)$ is recomputed over $\mathbb{Z}$ by the two routes described in the paper's
  Section 4 (pentagonal $q$-product then inversion; the logarithmic-derivative recurrence with each
  division by $n$ asserted exact); the run reports exact agreement at all 2567 coefficients,
  agreement with $W$ and with $Q$ **as decimal strings**, $a_{25}(2566)\bmod961=682$, $v_{31}=1$,
  and hence that (2) fails at $n=2$. One further check records that the cells $n=0,1$ *do* vanish
  mod 961 with $v_{31}$ exactly 2, as the paper's Section 5 states.
* **Section 3, 5 checks** — definition (1) pinned independently: a purely combinatorial
  colour-counting route agrees at $n=0..6$; $c=1$ reproduces OEIS A000041 at $p(50)$, $p(100)$,
  $p(200)$; $c=2$ reproduces OEIS A002513 up to $n=20$; and the prose reading $1/f_1^{25}$ is run as
  an anti-control, giving residues 550 and 492 mod 961 (the 550 is Remark 1's number).
* **Section 4, 7 checks** — a third, modular engine agrees with the exact routes on $0..2566$; the
  row $a_{25}(961n+644)\bmod961$ for $n=0..12$ is recorded, its vanishing cells are exactly
  $n=0,1,10$ (the window the paper's Section 5 refers to), and all ten nonvanishing cells are nonzero
  multiples of 31 — consistent with (5) while (2) fails. Also a neighbour row
  $a_{25}(641..647)\bmod961$; no class mod 961 identically zero over the swept range; and exactly the
  class $644\equiv24\pmod{31}$ identically zero, labelled in the output as a reproduction of
  Guadalupe's published theorem, not a discovery.
* **Section 5, 8 checks** — opposite polarity: the eight congruences the source *proves* must stay
  silent, and do, with 0 violations ($n=0..130$ for the mod-49 statement, $n=0..40$ for the seven
  mod-$p$ lines).
* **Section 6, 2 checks** — the other two displayed lines: $a_{41}(2465)\equiv1927\pmod{47^2}$ and
  $a_{61}(5044)\equiv2747\pmod{67^2}$, each a nonzero multiple of $p$ with the $n=0$ cell vanishing.
* **Section 7, 1 check** — from the source's own printed formula, $P(644)$ has fifteen residues
  rather than the single $\{644\}$ printed there, so the certificate window would be $15\cdot3814$
  cells.

The paper's Section 4 states that the run reports more than the paper uses and claims none of it.

## 3. What the program does not check

**The paper's theorem is a hand proof and the program is a control**, apart from Section 2: Theorem 2
is the division of one printed integer by another, and the only computation it needs is the single
identity $a_{25}(2566)=W$, which Section 2 of the run re-derives. Everything else in `verify.py` is a
control on the definition and on the decider.

The run states its own limits, verbatim from its closing block:

> NOTE SCOPE: this program re-derives every integer printed in the paper and nothing beyond it. It
> does NOT prove the conjecture on any other progression, it does NOT discharge Radu certificate
> hypothesis for the mod 31 companion congruence (that congruence is a published theorem of
> Guadalupe and is quoted, not proved, here), and it does NOT verify the source printed Radu
> parameter table against Radu original paper -- Section 7 re-derives P(644) from the formula AS
> PRINTED IN THE SOURCE only. The residue-class sweeps of Section 4 are bounded-range statements over
> the swept window, not uniqueness theorems; Sections 5 and 6 run to n = 40 and n = 1 respectively.

Additional points a referee should hold onto:

* Nothing is proved about $v_{31}$ elsewhere on the progression. Only the universally quantified
  statement (2) is refuted; the paper's Section 5 says (2) *does* hold at $n=0,1,10$ of the window
  checked, with $v_{31}=2$ at $n=0,1$, and makes no claim about how often that happens.
* No Ramanujan–Kolberg certificate is offered or completed in either direction, and (as the paper's
  Section 5 says) nothing checked would have proved (2) had the examined cells all vanished. The
  bound $\lfloor\nu\rfloor=3813$ quoted in Section 1 is transcribed from the source, not recomputed;
  Radu, *Ramanujan J.* **20** (2009), 215–251 is cited in the paper as not obtained.
* Equation (5) mod 31 is quoted from Guadalupe; the class sweeps only confirm it over the range
  swept.
* The two sibling lines are touched only mod $p^2$ at $n=1$; no exact integer behind either is
  computed, and the paper claims nothing from them.
* All quotations and parameter values are scoped by the paper to arXiv:2503.19399v2 (28 April 2026).
  The paper's Section 5 states that the version of record (*Eur. J. Math.* **12** (2026), no. 3,
  Paper No. 34; revised 26 May 2026, accepted 28 May 2026) is paywalled and was not read, so whether
  the accepted text prints the conjecture in this wording is not established.
* Two label mismatches to expect between transcript and paper, neither affecting a claim. The
  transcript's header and several check names call the statement "Conjecture 7.1", whereas the paper
  deliberately identifies it by position only; and some check names say "printed in the paper" of
  quantities the paper does not print — $a_{25}(0..12)$, the thirteen-cell residue row, the
  neighbour row $a_{25}(641..647)$ — which are controls belonging to the program, not tables in the paper. The paper prints
  the Radu bound as $\lfloor\nu\rfloor=3813$, while the Section 7 check line computes with 3814
  cells per residue.

## 4. How to check it

```sh
python3 verify.py            # one line per check, then the verdict; exits 0 iff all pass
shasum -a 256 verify.py      # 445acc462178b8d3932c0825d4f08406b34ef5bf3db7be6aab5d2a28c6b2fcd1
```

The header of `verify.output.txt` carries that same SHA-256 beside the recorded output, so the
transcript and the program in this folder can be paired; it also records the interpreter of the
recorded run, Python 3.9.25. No third-party package and no external data file is used. The digits $W$
and $Q$ are pasted into `verify.py` from the paper and compared against the fresh computation as
decimal strings, so a misprint in the paper would make the program fail.

Independently of the program, the proof can be checked in any computer algebra system, as the paper's
Section 2 says: take the coefficient of $q^{2566}$ in
$1/\bigl((q;q)_\infty(q^2;q^2)_\infty^{24}\bigr)$, a 2567-term truncation, and divide once by 31.
