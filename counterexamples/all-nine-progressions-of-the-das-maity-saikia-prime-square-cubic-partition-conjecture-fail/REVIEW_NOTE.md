# Referee note

**All Nine Progressions of the Conjecture of Das, Maity and Saikia on Generalized Cubic
Partitions Modulo Prime Squares Printed in arXiv:2503.19399v2 Fail**

The folder holds four files besides this note:

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper: 7 sections, Tables 1 and 2, one theorem, one lemma, two corollaries |
| `verify.py` | the verification program |
| `verify.output.txt` | its recorded run |

## 1. What the paper claims

Let $a_c(n)$ be given by $\sum_{n\ge0}a_c(n)q^n = 1/(f_1 f_2^{\,c-1})$, $f_k=\prod_{m\ge1}(1-q^{km})$
— equation (1) of §1: one uncoloured copy of every part, and $c-1$ further colours for the
*even* parts only. The unnumbered Remark of §1 records that the source's accompanying *prose*
instead describes $1/f_1^{\,c}$, and gives three reasons why (1) is the reading the authors
computed with.

The statement settled is quoted in full in §1: item 2 of the enumerated list in the concluding
remarks (Section 7) of Das, Maity and Saikia, printing there as **Conjecture 7.2**, as it stands
in the e-print arXiv:2503.19399v2. It asserts, for all $n\ge0$, nine progressions
$a_{p-6}(p^2n+r)\equiv0\pmod{p^2}$ over $p\in\{31,43,59,71,83\}$.

**Theorem 1.** All nine are false.
*(i)* The fourth line fails at $n=0$, the smallest index its own quantifier admits: equation (2),
$a_{65}(41)=451303321502143296879\equiv2170\pmod{71^2}$, and indeed $\equiv40\pmod{71}$, so it
fails already modulo the first power of the prime; the two residues are exhibited by the single
divisions (3) and (4).
*(ii)* Each of the other eight lines holds at $n=0$ and fails at $n=1$, seven of those eight
already modulo $p$. The eighteen values, with the $n=2$ column and the $n=1$ column reduced
mod $p$, are Table 2.

§3 replaces the large integer of part (i) by a small one. **Lemma 2** reduces $a_{p-6}$ mod $p$
to $b(m)=[q^m]\,f_2^{7}/f_1$ (its first step is Equation (14) of the cited Guadalupe paper, as §3
says), **Corollary 3** gives $a_{p-6}(N)\equiv b(N)\pmod p$ for $N<2p$, and the 21-term dot
product (5) of the two published rows A000730 and A000041 — all 21 products printed in Table 1,
positive ones summing to $1016099$ and negative ones to $1016130$ — gives $b(41)=-31$. Since
$71\nmid31$, the fourth line is false with no computer.

**Corollary 4** is the wider negative: for each of the five primes, with $c=p-6$, no residue $r$
modulo $p^2$ whatsoever makes $a_c(p^2n+r)\equiv0\pmod{p^2}$ hold at $n=0,1,2$, so no re-choice
of a printed residue rescues any of the five lines. §5 states that the criterion
$24r+13\equiv0\pmod p$ of (6) is Guadalupe's Theorem 3.2 and not the paper's, that no converse of
it is claimed or used, and that none of the nine printed residues satisfies it.

## 2. What the program checks

`verify.output.txt` records 55 `PASS` lines in ten labelled blocks, no failure, closing
`VERDICT: ALL 55 CHECKS PASS` with exit status 0. Standard library only, exact integer
arithmetic, no floating point, no randomness, as §7 states. Block by block, with the count of
checks and the claim each block bears on:

1. **The two printed tables** (3) — the rows $p(41-2j)$ and $e(k)=[x^k]\prod(1-x^n)^7$,
   $j,k=0..20$, of Table 1, re-derived from the products, plus an auxiliary row $p_{64}(j)$,
   $j=0..20$, from $\prod(1-x^m)^{-64}$, used by the independent exact route in block 3.
2. **Controls on the decider** (5) — the coefficient routine at $c=1$ against A000041 and
   $a_1(100)=190569292$, at $c=2$ against all 32 quoted terms of A002513, $a_c(0)=1$ for
   $c=1,2,3,25,65$, and the anti-control that the prose reading $1/f_1^{2}$ gives
   $1,2,5,10,20,36,\dots$, which is *not* A002513 (the Remark of §1).
3. **The refuting coefficient** (6) — the 21 printed products summing to
   $451303321502143296879$; that same integer again by an independent exact route; the residues
   $40\pmod{71}$ and $2170\pmod{71^2}$, both non-zero; and the divisions (3), (4). This is
   Theorem 1(i) and equation (2).
4. **The hand certificate** (6) — $\sum_{k=0}^{20}e(k)\,p(41-2k)=-31$; $71\nmid31$, agreeing
   with the exact integer; $b(41)=-31$ once more from the $b$-series; Lemma 2 numerically at
   $p=71$ for $N=0,\dots,1050$ with 0 mismatches; non-vacuity, 978 of those 1051 values being
   non-zero mod 71; and a control that must fire, the deliberately wrong exponent $f_2^{6}/f_1$
   mismatching at 1032 of the 1051 indices. This is §3.
5. **Companions** (2) — the exact integers $a_{25}(41)$ and $a_{65}(47)$, both genuinely
   $\equiv0$ modulo $31^2$ and $71^2$: the same computation that fails at $r=41$ passes here.
6. **The exhaustive census** (13) — for each of the five primes, the complete set of
   $r\in[0,p^2)$ with $a_c(r)\equiv0\pmod{p^2}$ and the survivor counts at depths 1, 2, 3; that
   $r=41$ is not a zero of $a_{65}$ mod $71^2$ at all; that no residue among the 18 221 classes
   gives three consecutive vanishing terms; and that the lone class surviving $n=0,1$ is
   $(p,r)=(31,644)$, which dies at $n=2$. This is Corollary 4.
7. **The nine printed cells** (11) — one check per cell, printing its residues at $n=0,1,2$,
   plus the mod-$p$ column $[25,22,26,38,0,57,16,39,14]$ and the summary that every one of the
   nine is false, one at $n=0$ and eight at $n=1$. These are the 27 entries of Table 2, its
   mod-$p$ column, and Theorem 1(ii).
8. **Positive controls on the same functions and primes** (5) — the source's proved mod-$p$
   congruences at $p=43,59,71,83$ over 129, 177, 213 and 249 consecutive values of $n$ with 0
   exceptions, and the identical census at $(c,p)=(3,7)$ returning exactly $[39]$ over all 49
   classes to depth 52, i.e. the source's proved $a_3(7^2n+39)\equiv0\pmod{7^2}$: a census that
   finds nothing here is not a census incapable of finding anything.
9. **The diagnostic** (3) — $(24r+13)\bmod p=[5,4,10,40,8,53,3,5,46]$ for the nine printed
   residues, none of them zero; the criterion holding for all seven residues of the proved mod-$p$
   theorem; and that it is the class $13(p^2-1)/24$, which at $p=71$ is 32, not the witness class
   41. This is §5.
10. **The prose anti-control on the nine cells** (1) — under $1/f_1^{\,c}$ the nine cells give
    $[514,234,558,561,1004,2860,4686,355,2764]$ at $n=0$, not one of them vanishing.

## 3. What the program does not check

For Theorem 1(i) and §3 **the program is a control, not the proof.** The refutation there is a
hand argument about one coefficient of a power series with 42 terms, checkable from Table 1 and
the two divisions (3), (4); the program re-derives it rather than carrying it. **Corollary 4 is
the one claim resting on computation alone**, and that computation is `verify.py`. Lemma 2 and
Corollary 3 are proved in §3 for every prime $p\ge7$ and every $N\ge0$; the program **samples**
those quantifiers only, at $p=71$ and $N\le1050$. Likewise the source's proved mod-$p$
congruences of block 8 hold for all $n\ge0$ but are checked over finitely many $n$.

The recorded run closes with its own `NOTE SCOPE` block; it is carried over here in full.

* **Nothing bibliographic is re-derived.** The locator of the conjecture (arXiv:2503.19399v2),
  the statement of Guadalupe's Theorem 3.2 and of Ahlgren's identity, and the OEIS entry numbers
  A000041 / A000730 / A002513 are **transcribed** from those sources, and are neither fetched nor
  checked. What *is* checked is the published **terms** of those three sequences, against the
  generating function (blocks 1 and 2).
* **The body of the version of record (Eur. J. Math. 12 (2026), Paper No. 34) has not been
  read.** No open-access copy exists; the refuted text is the e-print, and nothing here can
  certify that the published revision still prints the same residues. The paper says this itself,
  in the "Locator and scope" paragraph of §1 and as item (W1) of §6, and bounds it two ways: the
  block carrying the conjecture is byte-identical in v1 and v2, and by Corollary 4 every
  residue-level variant of the statement is false as well.
* **The census stops at depth 3 and at these five primes.** All 18 221 classes for
  $p\in\{31,43,59,71,83\}$ with $c=p-6$ are tested, but at $n=0,1,2$ only, so the honest form of
  the negative is exactly Corollary 4; it says nothing about other pairs $(c,p)$ — indeed
  $a_3(7^2n+39)\equiv0\pmod{7^2}$ is a true instance of the same shape, which block 8
  rediscovers. This is items (W2) and (W3) of §6. Depth is not a limitation of the refutation
  itself: Theorem 1(i) and §3 carry no depth dependence.
* **The primes 47, 67, 79 of the source's other statements are covered only by the arithmetic
  diagnostic, not by a census**, and Conjecture 7.1 — a separate statement, whose residue 644 is
  the lone two-term survivor found in block 6 — is not adjudicated at all.
* **How the nine residues came to be chosen is not computed.** The run lists that inference (a
  search cutoff) as not re-run; it is an inference about the authors, and nothing above depends
  on it.

§6 adds that Theorem 1 is a finite verification of explicitly exhibited coefficients, not a
structural result, and that nothing in the paper bears on the proved theorems of the source.

## 4. How to check it

```sh
python3 verify.py            # one PASS line per check; exits 0 only if all 55 pass
shasum -a 256 verify.py      # 1f17d405b18ce24e6412d182b0b4faa667c37ec2a28c6db029ff99a53ca37b8e
```

The digest above was computed from the shipped `verify.py`. The header of `verify.output.txt`
names the program and prints the SHA-256 beside it, so transcript and program can be paired
before the output is read; the program's own output begins at the `interpreter:` line, and the
recorded run was made with Python 3.9.25 and exited with status 0. The program's docstring notes
that the runtime is dominated by the exhaustive residue census and that a couple of minutes and a
few hundred MB of memory should be expected.

Theorem 1(i) needs no program at all: Table 1 prints the 21 products, their positive sum
$1016099$ and negative sum $1016130$, hence $b(41)=-31$, and $71\nmid31$ — checkable by hand from
the two published rows.
