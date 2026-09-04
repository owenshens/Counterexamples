# Review note

Paper: *A $3\times3$ Counterexample Family to Conjecture 4.5 of Chen--Cao's `arXiv:1801.00225v1`*

Files in this folder: `paper.tex` and `paper.pdf` (the paper: five sections, one numbered
statement, Theorem 1, one table, Table 1, and eleven displayed equations), `verify.py`
(the accompanying program), `verify.output.txt` (the recorded run of it). Nothing outside
this folder is needed to check the paper.

## 1. What the paper claims

Let $\omega_3$ be the $3\times3$ doubly substochastic matrices (nonnegative, every row sum
and every column sum at most $1$), $\sigma(A)$ the sum of all entries, and
$\omega_3^{\,s}=\{A\in\omega_3:\sigma(A)=s\}$. Conjecture 4.5 of Chen and Cao,
*On the maximum of the permanent of $(I-A)$*, asserts a two-branch closed
form for $\max\{\mathrm{per}(I-A):A\in\omega_3^{\,s}\}$ on $2<\sigma\le3$
— the single cell their Theorem 3.1, equation (1) of the paper, leaves open. The paper
refutes it. **Theorem 1:** for every $s$ with $2<s<3$ the matrix

$$A(s)=\begin{pmatrix}0&1-u&u\\1-u&0&u\\u&u&0\end{pmatrix},\qquad u=\tfrac{s-2}{2}$$

of equation (7) lies in $\omega_3^{\,s}$ and
$\max\{\mathrm{per}(I-A):A\in\omega_3^{\,s}\}\ \ge\ f(s)\ >\ V(s)$, where
$f(s)=(s^{3}-5s^{2}+4s+12)/4$ is $\mathrm{per}(I-A(s))$ (equation (9)) and $V$ is the
conjectured value, equation (6). Equivalently, Chen and Cao's own Lemma 4.4 — of which the
conjecture is exactly the assertion of tightness — is not tight anywhere on $(2,3)$. The
excess over the governing first branch is $s(s-3)^{2}/4$; at $s=5/2$ the conjecture gives
$23/16$ and the witness $51/32$.

Two qualifications belong to the claim; both are argued in Section 1 of the paper and are
not added here.

* **Which reading is refuted.** Equation (3), the conjecture as printed, is a maximum over
  all of $\omega_3$ with a right-hand side depending on $\sigma(A)$, which read literally is
  the single number $2$. The paper refutes the fixed-mass reading, with the printed second
  branch $6-4\sigma$ read as the $6-2\sigma$ that Lemma 4.4 proves, i.e.
  $V(s)=\max\{(s^{2}-5s+12)/4,\ 6-2s\}$ — the strongest reading. The paper states that
  refuting the printed branch would be worth nothing and does not claim it as content; the
  witness happens to exceed $6-4s$ as well, and at $s=5/2$ the governing branch is the first,
  printed identically in the proved Lemma 4.4, so there the readings coincide.
* **Which text is refuted.** Only the e-print `arXiv:1801.00225v1` was read. The journal
  version of record was not obtained, and the paper says so in the abstract, in Section 1 and
  in Section 4; every statement about what the source prints is a statement about v1.

The paper does **not** claim the true maximum. Section 4 gives the bracket
$51/32\le\max\le2$ at $s=5/2$ and claims no optimality of $A(s)$ on any slice or over the
whole polytope.

## 2. What the program checks

`verify.output.txt` records **97 checks, all passing** (`VERDICT: ALL 97 CHECKS PASS`),
in fourteen numbered steps, exact rational arithmetic throughout, permanents by the full
$n!$ expansion. By block:

* **Steps 1–3 (8+8+8).** At the eight masses of Table 1 ($21/10$, $9/4$, $5/2$, $14/5$,
  $11/4$, $29/10$, $2$, $3$): $A(s)\in\omega_3$ with row *and* column sums printed, mass
  exactly $s$, and $\mathrm{per}(I-A(s))=f(s)$ by the $3!=6$ expansion. This is the
  membership-and-value half of Theorem 1.
* **Step 4 (1+8).** Table 1 reproduced row by row — permanent, $V(s)$, governing branch,
  excess — together with the rational bracket
  $2274917/1000000<(-3+\sqrt{57})/2<1137459/500000$, enough to decide the branch at every
  tabulated mass.
* **Step 5 (6).** The strict inequality of Theorem 1 at the six *interior* masses: the
  permanent beats the first branch, the corrected $6-2s$, and the printed $6-4s$.
* **Step 6 (2).** Negative-polarity control at the two masses where the truth is published
  ($s=2$: equation (1) at $n=3$, $e=2$; $s=3$: Chen–Cao's Lemma 4.2, the value over
  $\Omega_3$) — the family returns the published value, excess exactly $0$, no violation.
  These are the last two rows of Table 1.
* **Steps 7–8 (4+5).** Controls on the source's own arithmetic: their $A_0$ of equation (5)
  reproduces their printed first branch exactly; their $A_1$ gives $6-2\sigma$ at every mass
  and $6-4\sigma$ at none; $6-4s<0$ for all $s>3/2$; and the printed threshold is the
  crossover of the first branch with $6-2s$, the crossover with $6-4s$ being $s=1$. These
  support the misprint reading of Section 1.
* **Step 9 (3).** The two identities of the proof of Theorem 1 — equation (10),
  $f(s)-(s^{2}-5s+12)/4=s(s-3)^{2}/4$, and equation (11), $f(s)-(6-2s)=t+t^{2}/4+t^{3}/4$
  with $t=s-2$ — by exact coefficient comparison rather than by sampling; the third check
  records that positivity of the two gaps is evaluated only at the six interior masses.
* **Step 10 (1+8).** The consistency paragraph after Theorem 1: the identity
  $2-f(s)=(s-2)(2+3s-s^{2})/4$, and $f(s)\le2$ at all eight masses, against the two proved
  upper bounds giving $\max\le2$ at $n=3$.
* **Steps 11–13 (13+10+10).** The symmetric zero-diagonal slice (its vertices are $A(s)$,
  $A_0$ is an edge midpoint, gap $s(s-3)^{2}/4$, centroid the slice minimum, restricted
  Hessian positive definite, its closed form verified entry-by-entry at 120 exact
  point-direction pairs), and a direct-sum lift with
  the full $n!$ permanent at $n=3,5,7$ plus lifted-endpoint consistency at odd $n=3,5,7,9,11$.
  Section 5 of the paper flags this material as supporting no claim made in the paper, and
  the transcript's own Step 13 heading calls its endpoint agreement non-discriminating.
* **Step 14 (2).** The headline bracket $51/32\le\max\le2$ at $s=5/2$, with excess exactly
  $5/32$ and ratio $51/46$.

## 3. What the program does not check

**The load-bearing verification is by hand; the program is a re-derivation, not evidence**
— Section 5 of the paper says exactly that, and three of the fourteen steps (6, 7, 8) are
labelled CONTROL in the transcript rather than checks of a claim. Theorem 1's universal
statement over $(2,3)$ rests on the two polynomial identities of Section 3, whose strict
signs on the interval are read off by hand.

Carried over from the closing SCOPE and NOT RE-RUN block of `verify.output.txt`, and from
Sections 4 and 5 of the paper:

* Membership in $\omega_3$, the exact mass and the permanent are verified at **eight sampled
  masses only**, and the strict refutation $f(s)>V(s)$ at the **six interior masses only**.
  Several check names and step headings are phrased over the whole open interval; what is
  evaluated at those places is exactly those samples. The paper states this too.
* At $s=2$ and $s=3$ the run establishes $f(s)=V(s)$ with excess $0$ — a control, not a
  violation.
* The two identities are verified universally by coefficient comparison, but the strict
  **sign** of the two gaps for all $2<s<3$ is **not** machine-verified.
* The program does **not** establish $\max\{\mathrm{per}(I-A):A\in\omega_3^{\,s}\}=f(s)$;
  that upper bound is open (Section 4). Step 11 establishes optimality of $A(s)$ only on the
  symmetric zero-diagonal slice, not over the full nine-entry polytope, and the paper claims
  no optimality at all.
* Not re-run: any global search, and any floating-point multistart search. No decision
  anywhere in the program is taken on a float.
* The program reads no source text. What Chen and Cao print — including the variant branch
  $6-4s$ used in Step 5, transcribed by hand — is quoted in Section 1, and the journal text,
  *Linear Algebra Appl.* **555** (2018) 412–431, was never obtained. The source-file digest,
  the citer count and the bibliography comparison discussed in Sections 1 and 4 are likewise
  not recomputed by anything here.
* Malek's bound $\mathrm{per}(I-A)\le2^{\lfloor n/2\rfloor}$, used only to give $\max\le2$ at
  $n=3$, is not reproved. Section 4 and the bibliography state that it was read only through
  its zbMATH review and that the article itself was not obtained.

## 4. How to check it

By hand: Section 2 exhibits a $3\times3$ matrix with entries in $\{0,\tfrac14,\tfrac34\}$,
six row and column sums, and a six-term permanent; Section 3 is two polynomial identities in
one variable. No computation is needed for any claim in the paper.

The program requires Python 3.9 or later and the standard library only (`fractions`,
`itertools`), with no external data file:

```sh
shasum -a 256 verify.py
# 14ca606c8089ce9bc3f6e1b37eb3d1fe1b3a23ec5b876fa862197521194354ef
python3 verify.py
```

That digest was computed from `verify.py` as shipped here, and the header of
`verify.output.txt` prints the same SHA-256 beside the program name, so transcript and
program can be paired before either is trusted. The recorded run used Python 3.9.25 and
exited with status 0.
