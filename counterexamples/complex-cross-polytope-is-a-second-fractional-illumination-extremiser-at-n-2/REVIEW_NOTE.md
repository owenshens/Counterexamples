# Referee note

Paper: *The Complex Cross-Polytope $B_1(\mathbb{C}^2)$: Real Fractional Illumination
Number $2^2$, and Not a Complex-Linear Image of the Polydisc*.

Files here, and nothing else is needed: `paper.tex` / `paper.pdf` (the paper), `verify.py`
(the accompanying program), `verify.output.txt` (its recorded run).

## 1. What the paper claims

**Theorem 1.** For $K=B_1(\mathbb{C}^2)=\{z\in\mathbb{C}^2:|z_1|+|z_2|\le1\}$,
$\mathrm{ill}^*(K)=4=2^2$, and $K$ is not a complex-linear image of the polydisc $D^2$.

The statement addressed is **Conjecture 1.2** of Rotem, Schejter and Slomka, *The complex
illumination problem* (Combinatorica **46** (2026), no. 1, Paper No. 3; arXiv:2410.12021),
quoted verbatim in Section 1 ("The conjecture"): every complex convex body
$K\subseteq\mathbb{C}^n$ has $\mathrm{ill}^*(K)\le2^n$, with equality if and only if $K$ is
a linear image of $D^n$. Theorem 1 attacks the **equality conjunct at $n=2$** only: the
bodies in $\mathbb{C}^2$ attaining $2^2$ include at least the two $GL_2(\mathbb{C})$-orbits
of $D^2$ and of $B_1(\mathbb{C}^2)$.

Two readings, fixed in Section 1, are load-bearing. $\mathrm{ill}^*$ is read as the ordinary
**real** fractional illumination number of $K$ regarded as a convex body in $\mathbb{R}^4$ —
real directions, no circle invariance imposed on the measure — and "linear image" is read
**complex**-linearly. The paper states that it does not reproduce the source's own wording
of $\mathrm{ill}^*$, that every value is asserted for this reading only, and that it has not
established that the two agree; under a reading restricting $\mu$ to circle-invariant
measures the numbers could differ. The exhibited measure is supported on unit vectors, so
admissibility does not turn on whether non-unit directions are allowed.

Section 2 ("The result") lists what is *not* claimed: the bound $\mathrm{ill}^*(K)\le2^n$ is
not contradicted (the witness attains it with equality, and that conjunct remains open);
whether $B_1(\mathbb{C}^2)$ is a *real*-linear image of $D^2$ is not decided; nothing is
claimed for $n\ge3$; whether $\{D^2,B_1(\mathbb{C}^2)\}$ exhausts the equality set at $n=2$
is open; and Conjecture 1.1 of the same source, on the classical illumination number
$\mathrm{ill}$, is untouched.

All proofs are closed form. Section 3 ("The witness") records, in items (a)–(f), that $K$ is
a complex convex body with $K=\mathrm{conv}(De_1\cup De_2)=(D^2)^\circ$, its support
function $h_K(u)=\max(|u_1|,|u_2|)$, the exhaustive three-class boundary, the extreme set
$\mathrm{ext}(K)=C_1\cup C_2$ (two disjoint circles), the mass-$4$ measure (1) and the
first-order illumination criteria (2). Section 4 proves $\mathrm{ill}^*(K)\ge4$ (the extreme
circles are illumination-disjoint; an illuminating direction sees an arc of length
$2\arccos(c)\le\pi$; averaging over $\alpha$ with Tonelli gives $\mu(V_1)\ge2$ and
$\mu(V_2)\ge2$ on disjoint sets) and $\mathrm{ill}^*(K)\le4$ (the exhibited measure, checked
on each boundary class and the interior). **Proposition 2** (Section 5) proves $K$ is not a
complex zonoid by a mean-value obstruction; a complex-linear image of $D^2$ is a
$2$-generator complex zonotope, so Theorem 1's second assertion follows. **Remark 3**
records that no theorem of the source is contradicted. **Proposition 4** (Section 6) gives
$\mathrm{ill}^*(B_1(\mathbb{C}^n))=2n$ for every $n\ge1$, so the witness lies strictly below
$2^n$ for $n\ge3$; Section 6 adds that no priority check was made for that formula and
claims no novelty for it.

## 2. What the program checks

`verify.py` uses exact rational arithmetic and the standard library only; no floating-point
decision is made. Sums of two square roots are compared with rationals by exact
elimination, and $\pi$ enters only as a rational enclosure from Machin's formula. The
samples are 36 rational unit-circle points from the Pythagorean triples 3-4-5, 5-12-13,
8-15-17, 7-24-25. The run ends `VERDICT: ALL 59 CHECKS PASS` and
`=== program exited with status 0 ===`, in twelve blocks:

* **Step 0** (3 checks): the object block parses ($2+2=4$; claimed value $4=2^2$) and the
  square-root eliminations reproduce known values.
* **Step 1** (5) — Section 3(a): the defining norm, circle invariance on 576 (rotation,
  point) pairs, central symmetry, real dimension $4$ (four points of $K$ with determinant
  $1/16$), and $K=\mathrm{conv}(De_1\cup De_2)$ on samples.
* **Step 2** (3) — Section 3(b): $h_K=\max(|u_1|,|u_2|)$ on 64 test points of $K$, attained
  at an axis extreme point, and $h_{D^2}=|u_1|+|u_2|$, the polar situation.
* **Step 3** (5) — Section 3(c),(d): the three boundary classes are exhaustive on samples
  (axis1 216, axis2 216, mixed 432, unclassified 0); points of $C_1$ are exposed; sampled
  mixed points are proper combinations, hence not extreme; the mixed face is a
  $1$-dimensional segment, while the face of $D^2$ at a torus normal is a point.
* **Step 4** (5) — Section 3(f): the three criteria (2) are tested against the **definition**
  of illumination (membership of $x+tv$ in the interior) on $1136+848+576$ (boundary point,
  direction) pairs, and the first-order reduction is checked not to drop the term $t|v_2|$.
* **Step 5** (7) — Section 4, lower bound: 600 directions illuminating a point of $C_1$ all
  have $|v_1|>|v_2|$, 532 illuminating $C_2$ all have $|v_2|>|v_1|$, no sampled direction
  serves both circles, the arc shrinks strictly as $c$ grows with endpoints at
  $\cos\alpha=-c$; one check records that the cone mass $2$ is hard-coded; then $2+2=4$.
* **Step 6** (8) — Section 4, upper bound: an open half-circle has normalised measure $1/2$
  by the antipodal involution; a point of $C_1$ receives exactly $1$, a point of $C_2$
  exactly $1$, a mixed point $1+1=2$, an interior point $4$; total mass $4$.
* **Step 7** (8) — Proposition 2: the rational enclosure of $\pi$,
  $|a+e^{it}a|^2=(2a|\cos(t/2)|)^2$, $\int_0^{2\pi}|\cos(t/2)|\,dt=4$,
  $M(r,r)/r=4/\pi\in(1.2732,1.2733)>1$, $M(a,0)=a>a/2$, and that $T(D^2)$ is a
  $2$-generator complex zonotope.
* **Step 8** (6): controls of both polarities — the polydisc criterion on 204 pairs of
  $\partial D^2$; the same averaging identity reproducing the published
  $\mathrm{ill}^*(D^2)=2^2$ from both sides; the published $2=2^1$ at $n=1$; and the complex
  Euclidean ball, whose criterion matches interior membership on 910 pairs of $S^3$ and
  whose value $2<4$ shows the method separates.
* **Step 9** (4) — Proposition 4: the $n$ extreme circles of $B_1(\mathbb{C}^n)$ are
  pairwise illumination-disjoint on samples for $n=2,\dots,8$ (176 single-condition
  witnesses); $2n=2^n$ exactly at $n=1,2$ and $2n<2^n$ for $3\le n\le30$; and the real
  anti-control, that $B_1(\mathbb{R}^2)$ *is* a linear image of the square ($\det=-1/2$).
* **Step 10** (3): three directions at $120^\circ$ cover $C_1$; none of 25 sampled centre
  pairs does. (Outside the paper's claims — see §3, last bullet.)
* **Step 11** (2) — Theorem 1: $4\le4$, so the bound conjunct is satisfied with equality,
  not violated; and the equality clause fails at $n=2$.

## 3. What the program does not check

**The main theorem is a hand proof and the program is a control, not an independent
verification.** Section 7 of the paper says exactly that, and the run's closing SCOPE note
repeats it. In detail:

* Every check whose name ends `_on_samples` tests a closed-form criterion against the
  definition on finite sets of rational unit vectors: a test of the criterion, **not** a
  proof of it. None of the universal statements the theorem rests on is established — the
  criteria at *every* boundary point; $\mathrm{ext}(K)=C_1\cup C_2$ (only sampled mixed
  points are shown non-extreme, while the inclusion $\mathrm{ext}(K)\subseteq C_1\cup C_2$
  and the extremality of every point of the two circles stay closed-form arguments in the
  paper); $\mu(A_K(x))\ge1$ for *every* $x$; and $\mathrm{ill}^*(B_1(\mathbb{C}^n))=2n$.
* The cone mass $\mu(V_1)\ge2$ is **asserted, not derived**: the Step 5 check says so, and
  the averaging derivation
  $2\pi\le\int_0^{2\pi}\mu(A_K(x_\alpha))\,d\alpha=\int_{V_1}2\arccos(|v_2|/|v_1|)\,d\mu\le
  \pi\,\mu(V_1)$ is imported from the paper; only $2\cdot2=4$ is checked. The closing check
  of Step 6 says the same on the other side: the two bounds meet arithmetically at $4$
  **given** that cone-mass bound and the paper's closed-form universality of
  $\mu(A_K(x))\ge1$.
* In Step 7 the mean-value route is closed only *given* imported ingredients: the
  representation $h_K=\int M(|x_1|,|x_2|)\,d\nu$, the inequality
  $\int M\,d\nu\ge\int(|x_1|+|x_2|)/2\,d\nu$, the equality-almost-everywhere step, and the
  strict inequality $M>\max$ at every pair of positive moduli. Verified there are only
  $1=(1+1)/2$, $M(r,r)/r=4/\pi\in(1.2732,1.2733)>1$ and $M(a,0)=a>a/2$. The **second** route
  to the same conclusion imports the source's zonoid face lemma (a complex zonoid has no
  $1$-dimensional proper face) and is conditional on it.
* The topological step "a $2$-torus is not a disjoint union of two circles" is a
  homeomorphism statement and is **not** machine-checked; the SCOPE note observes that the
  zonoid route of Step 7 reaches the same conclusion and is the one the paper leans on.
* Nothing enumerates complex convex bodies, so the program does not claim that
  $\{D^2,B_1(\mathbb{C}^2)\}$ exhausts the equality set at $n=2$; and nothing decides any
  case $n\ge3$ beyond the arithmetic $2n<2^n$. Step 9 states its own conclusion as
  *assuming* Proposition 4's formula and records that $\mathrm{ill}^*(B_1(\mathbb{C}^3))$
  and $\mathrm{ill}^*(B_1(\mathbb{C}^4))$ are not computed.
* Transcribed from the cited source rather than recomputed: in the Euclidean-ball control of
  Step 8, the hemisphere measure $1/2$ and the matching lower bound. By contrast
  $\mathrm{ill}^*(D^2)=2^2$ *is* re-derived there by the same averaging identity and
  labelled as reproducing the published value.
* No program can settle the reading of $\mathrm{ill}^*$, or the intended sense of "linear
  image", on which everything above is conditional. That is what a referee must check
  against the source itself; the paper flags it in Section 1 and in the first bullet of
  Section 2.
* Step 10 concerns the **classical** illumination number, about which the paper claims
  nothing. `verify.py`'s object block nevertheless carries
  `claimed_classical: ill(K) = 6 < 7`, and that block's third check states that "this
  program does not establish ill(K) = 6", the necessity of three directions per circle being
  an imported arc argument.

## 4. How to check it

```sh
shasum -a 256 verify.py
python3 verify.py            # Python 3.9 or later, standard library only, no data file
```

One line per check is printed, then the verdict; the program exits $0$ if and only if every
check passes. The header of `verify.output.txt` carries the SHA-256 of the program beside
the output, so the two can be paired: it reads
`sha256:  a7e21243d1309b6f1fe3eb84fabfed6976946348306429830dfecb486f838dd8` and
`python:  Python 3.9.25`. Digests of the shipped files, computed with `shasum -a 256`:

| file | SHA-256 |
|---|---|
| `verify.py` | `a7e21243d1309b6f1fe3eb84fabfed6976946348306429830dfecb486f838dd8` |
| `verify.output.txt` | `24701b7c0046c4349c0ba0da7c59da040ae2286857e807cf049ea5dfc5219ecd` |
| `paper.tex` | `fb1a0a8908f8e9d6f5d1584e5aabf75f635c91bc031916c32ac526abfd990210` |
| `paper.pdf` | `80c504e84be51de6194afcbd39e4a734087dbedc95555fc7d4c66e0b7f8ce210` |
