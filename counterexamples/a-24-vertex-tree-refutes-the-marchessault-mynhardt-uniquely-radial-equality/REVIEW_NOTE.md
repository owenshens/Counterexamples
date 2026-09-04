# Referee note: *A 24-Vertex Tree Refuting the Marchessault--Mynhardt Uniquely-Radial Equality*

Besides this note the folder holds four files, and nothing outside it is referred to below:

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper (7 pages) |
| `verify.py` | the verification program |
| `verify.output.txt` | the recorded output of one run of `verify.py` |

## 1. What the paper claims

Problem 1 of the paper quotes Problem 19 (p. 96, Section 5, Open Problems) of
E. M. Marchessault and C. M. Mynhardt, *Lower boundary independent broadcasts in trees*,
Discuss. Math. Graph Theory **44** (2024), 75--99, stated there to be Problem 2 of the
e-print `arXiv:2105.04611v2`: if the radial subtrees of a tree $T$ are uniquely radial
with radii at least $2$ for any maximum split-set $M$, is
$i_{\mathrm{bn}}(T)=\gamma_b(T)+\lceil(|M|+1)/3\rceil$?

**Theorem 2 answers no.** The counterexample, built in Section 2 (*The tree*), is the tree
$T$ on $24$ vertices obtained by chaining four copies of the spider $S(2,2,1)$: spine
$v_0v_1\cdots v_{19}$ with pendants $u_1,\dots,u_4$ at $v_2,v_7,v_{12},v_{17}$, with the
$23$ edges printed in full as (3), and with $\operatorname{diam}T=19$,
$\operatorname{rad}T=10$. Section 3 (*$T$ satisfies the hypothesis*) places $T$ inside the
hypothesis: Lemma 3 that the diametrical path is unique, Lemma 4 that any split-set has
cardinality congruent to $\operatorname{diam}T$ mod $2$, Lemma 5 that the maximum
split-set is unique and equals $\{v_4v_5,\,v_9v_{10},\,v_{14}v_{15}\}$ so $|M|=3$, and
Lemma 6 that each of the four radial subtrees is a copy of $S(2,2,1)$, uniquely radial of
radius exactly $2$. Equation (5) then gives
$\gamma_b(T)=(\operatorname{diam}T-|M|)/2=(19-3)/2=8$ from the Herke--Mynhardt formula (1)
as quoted from the source, so by (6) the problem demands the value
$8+\lceil 4/3\rceil=10$.

Section 4 (*The certificate*) supplies the witness: Lemma 8 shows that $f\equiv 1$ on
$S=\{v_0,v_2,v_4,v_7,v_9,v_{12},v_{14},v_{17},v_{19}\}$ is a *maximal* bn-independent
broadcast of cost $9$ --- maximality via the criterion quoted as Proposition 7 --- whence
$i_{\mathrm{bn}}(T)\le 9$. With the quoted inequality $\gamma_b\le i_{\mathrm{bn}}$ this
gives $i_{\mathrm{bn}}(T)\in\{8,9\}$, and both values differ from $10$.

No theorem of the source falls, and the paper says so after Theorem 2: the published
bound (2) asserts "$\le$" only and $9\le 10$; the source's Theorem 6, which gives equality
for $|M|\in\{1,2\}$, is untouched because $|M|=3$ here. Only the conjectured equality
fails. The paper claims **no** exact value of $i_{\mathrm{bn}}(T)$ and proposes no repair
of the problem.

## 2. What the program checks

`verify.output.txt` records one run of `verify.py`: one `PASS` line per check, then
`VERDICT: ALL 45 CHECKS PASS`, then a closing `NOT RE-RUN HERE` paragraph, and a final line
recording exit status 0. By block:

* **6 checks --- Section 2 and Lemma 3.** $n=24$, $m=23$, connected; the vertex set is
  exactly the printed spine plus the four pendants; the degrees are as printed (leaves
  `u1 u2 u3 u4 v0 v19`, degree 3 at `v2 v7 v12 v17`, the rest degree 2);
  $\operatorname{diam}T=19$ realised by $(v_0,v_{19})$; $\operatorname{rad}T=10$ attained
  exactly at `v9 v10`; and $(v_0,v_{19})$ is the only diametrical pair, which is Lemma 3
  and is what confines split-sets to the spine.
* **2 checks --- audit of the split-set enumerator.** It prunes on the
  even-positive-length condition, so it is compared against exhaustive enumeration of
  *all* edge subsets of the diametrical path on two small trees: `P12` (pruned 5 subsets,
  brute force 5) and the two-block chain (1 and 1). The paper's Exact verification
  paragraph announces exactly this cross-check.
* **3 checks --- Lemma 5.** Maximum split-set cardinality $|M|=3$ with exactly one maximum
  split-set (cardinalities present `[1, 3]`, five split-sets in all); that it is the
  printed `v4v5 v9v10 v14v15`; and that $T-M$ has four components of six vertices each.
* **1 check --- Lemma 6.** Each of the four radial subtrees has $\operatorname{rad}=2$,
  $\operatorname{diam}=4$, $\gamma_b=2$ and exactly one $\gamma_b$-broadcast.
* **4 checks --- (1), (5) and (6).** The printed strength-$2$ broadcast at the four block
  centres is dominating of cost $8$, so $\gamma_b(T)\le 8$; formula (1) gives $(19-3)/2=8$,
  agreeing with that witness; $T$ is non-radial (it has a split-set, and
  $\gamma_b=8<10=\operatorname{rad}T$); and the demanded value is $8+\lceil 4/3\rceil=10$.
* **6 checks --- Lemma 8, step for step.** $f$ is a legal broadcast of cost $9$; its
  support is independent; it is dominating (no unheard vertex among the $24$); it is
  bn-independent and the multiply-heard vertices are exactly `v1 v3 v8 v13 v18`; $U_f^E$ is
  exactly `v5v6 v10v11 v15v16`; and $T-U_f^E$ has four components carrying $3,2,2,2$
  broadcasters, which is the hypothesis of Proposition 7.
* **3 checks --- the conclusion of Theorem 2.** $i_{\mathrm{bn}}(T)\in\{8,9\}$ against the
  predicted $10$; both admissible values differ from $10$; and (2) is not contradicted,
  since $9\le 10$.

That is **25** checks on the paper's own claims. Two further checks re-establish the
maximality of $f$ by routes the paper does not use: a second criterion of the source
($B_f(v)-PB_f(v)$ non-empty at all nine broadcasters, sizes `[1, 2, 1, 1, 1, 1, 1, 1, 1]`),
and a direct test over all $24$ vertices that no unit increment of $f$ keeps it
bn-independent. The direct test does not depend on either quoted criterion.

The remaining **18** checks are on objects the paper's Scope section says it does not
discuss and does not rely on --- and a referee may ignore all eighteen:

* an anti-control (1 check): `P5` also has radius $2$ and $\gamma_b=2$ but four
  $\gamma_b$-broadcasts, so is not uniquely radial;
* a second cost-$9$ maximal bn-independent broadcast on $T$ of a different shape
  (5 checks): $g(v_4)=4$, $g(v_9)=1$, $g(v_{13})=2$, $g(v_{17})=2$, with
  $U_g^E=\{v_{10}v_{11}\}$ and multiply-heard set `v8 v15`;
* nine pairwise disjoint closed neighbourhoods (2 checks), giving
  $\gamma(T)=i(T)=9=\gamma_b(T)+1$;
* two controls (4 checks): the two-block chain, where $|M|=1$, $\gamma_b=4$ and no maximal
  bn-independent broadcast of cost $\le 4$ exists, so $i_{\mathrm{bn}}=5=\gamma_b+1$ and
  the equality *holds*; and `P12`, where $\gamma_b=4$, $|M|=3$ and $i_{\mathrm{bn}}=5<6$,
  outside the hypothesis because its radial subtrees are copies of `P3` of radius $1$;
* `T5` (6 checks), the $30$-vertex chain of five blocks, with $\operatorname{diam}=24$ and a
  unique maximum split-set of the **even** cardinality $|M|=4$: there $\gamma_b=10$, the
  predicted value is $12$, and a printed cost-$11$ maximal bn-independent broadcast gives
  $i_{\mathrm{bn}}(\texttt{T5})\in\{10,11\}$, so the equality fails there too.

## 3. What the program does *not* check

The run states most of this itself, in the closing `NOT RE-RUN HERE` paragraph of six
numbered items, which agrees with the paper's Scope section. Carried over below in seven
items; the framing of item 1 and the closing remarks of items 3 and 4 are added here.

1. **Theorem 2 is a hand proof; the program is a control on its finite ingredients.**
   Every quantity the proof uses is recomputed from the printed edge list (3) and the
   printed broadcasts, but no step of the argument is replaced by a search.
2. **No exhaustive search over the broadcasts of $T$ (or of `T5`) is performed, so no exact
   value of $i_{\mathrm{bn}}$ is established** --- only $i_{\mathrm{bn}}(T)\in\{8,9\}$ and
   $i_{\mathrm{bn}}(\texttt{T5})\in\{10,11\}$, which the run and the Scope section both say
   is all the refutation uses, since $8\ne 10$ and $9\ne 10$ alike. (Exhaustive search over
   cost $\le 4$ *is* performed, but only on the two control trees `P12` and the two-block
   chain.)
3. **The two results carrying the arithmetic are transcribed from the source, not proved.**
   The lower bound $\gamma_b(G)\le i_{\mathrm{bn}}(G)$ (attributed to the source's
   Proposition 3) and the formula $\gamma_b(T)=(\operatorname{diam}T-|M|)/2$ (the source's
   Theorem 1) are cited results; only their arithmetic on these trees is checked. The value
   $8$ is corroborated only *from above*, by the exhibited cost-$8$ dominating broadcast; no
   independent lower bound $\gamma_b(T)\ge 8$ is computed. Since the demanded value is
   $\gamma_b(T)+2$, this transcription is load-bearing and a referee should confirm
   formula (1) against the source.
4. **The maximality criteria are transcribed as criteria, not derived.** They are applied
   and cross-checked against each other and, on $f$, against the direct
   no-single-increment test. Note a mismatch a referee should be aware of: the paper quotes
   **one** criterion, Proposition 7 in Section 4, whereas the run's closing paragraph
   speaks of "the two maximality criteria the paper quotes"; the second criterion appears
   only in the run, and no claim of the paper turns on it.
5. **No bibliographic or attribution claim is checked:** the wording of Problem 1 and its
   numbering in either version of the source, and the paper's statement --- made in its own
   Scope section, of a search it describes there as "not exhaustive" --- that no earlier
   work answers the problem and that no tree with these properties was found, including its
   remark that `doi:10.7151/dmgt.2546` answers neighbouring problems of the same list but
   not this one.
6. **The reading of "uniquely radial" is not, and cannot be, validated by the program.**
   Section 1 says the phrase carries no numbered definition in the source and fixes it ---
   on the strength of the source's discussion of the tree $T_1$ of its Figure 3, whose
   radial subtrees are declared not uniquely radial --- as $\gamma_b=\operatorname{rad}$
   attained by exactly one broadcast; the program counts $\gamma_b$-broadcasts under that
   reading. The paper itself calls this the one clause of the hypothesis that is a reading
   rather than a quotation, and it is where a referee's judgement is required.
7. **Nothing about the wider family or about other trees.** The run states that the $k$-fold
   family of which $T$ and `T5` are the members $k=4$ and $k=5$ is not treated, that
   nothing is checked about trees other than $T$, `T5` and the three small controls, and
   that no repair of the equality is proposed or tested. The Scope section agrees.

## 4. How to check it

```sh
python3 verify.py
```

Python 3.9 or later, standard library only; no arguments and no input file --- the trees
and the broadcasts are hard-coded from the paper. The program prints one `PASS` line per
check, then the verdict and the `NOT RE-RUN HERE` paragraph, and exits 0 only if every
check passes.

The opening lines of `verify.output.txt` are a provenance header that the program itself
does not print: they name the program, give its SHA-256, and record `Python 3.9.25` for
that run. That digest is what pairs the transcript with the program, so run

```sh
shasum -a 256 verify.py
```

and compare the result,

    d0251bc5b97d1c224f2f18567bd7c2aa415e1e4046cb71cb5f610ea80e8aa724  verify.py

with the `sha256:` line of `verify.output.txt`; for the shipped files they agree.
