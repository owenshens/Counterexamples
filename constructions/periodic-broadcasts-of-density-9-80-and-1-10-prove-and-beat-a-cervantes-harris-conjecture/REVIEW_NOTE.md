# Referee note

Paper: *Periodic Broadcasts of Density 1/10, 1/8 and 5/36: Confirming and Strengthening Three
Conjectured Density Bounds of Cervantes and Harris.*

The folder holds exactly five files: `paper.tex` and `paper.pdf` (the paper), `verify.py` (the
checking program), `verify.output.txt` (the recorded run of that program), and this note.

## 1. What the paper claims

The subject is $(t,r)$ broadcast domination of the infinite truncated square tiling graph
$H_{\infty,\infty}$, the $3$-regular Archimedean $4.8.8$ net. The reception rule is equation (1)
of §1, with a strict distance cutoff, so at $t=4$ one tower contributes $4,3,2,1$ at distances
$0,1,2,3$ and nothing at distance $\ge 4$; the density is the limit along an exhausting sequence
of finite subgraphs and $\delta_{t,r}$ is the infimum over broadcasts.

Cervantes and Harris conjecture, in the form quoted in §1, that
$\delta_{4,2}(H_{\infty,\infty})\le 9/80$, $\delta_{4,3}\le 1/7$ and $\delta_{4,4}\le 1/6$; the
bibliography entry records this as Conjecture 4.10 of the published version.

**Theorem 1** settles all three by exhibiting the lattice-periodic broadcasts $T_5$, $T_6$, $T_7$
of §3 ("The three witnesses"), each written as a congruence in the coordinates $a@(x,y)$ of §2,
and gives
$$\delta_{4,2}\le\tfrac1{10}<\tfrac9{80},\qquad
  \delta_{4,3}\le\tfrac18<\tfrac17,\qquad
  \delta_{4,4}\le\tfrac5{36}<\tfrac16 .$$
So each conjectured inequality holds and each conjectured value is improved. Because $\delta_{t,r}$
is an infimum, each part is an existence statement, so one periodic broadcast suffices; equation (3)
of §2 makes each witness a finite check, an $L$-periodic set with $L$ of index $s$ having density
exactly $k/4s$ with no boundary correction and $f$ constant on each of the $4s$ classes. The paper
prints the whole table of $f$-values for each witness ($20$, $16$ and $72$ entries).

**Theorem 2**, proved by hand in §4 from the capping **Lemma 4**, gives $\delta_{4,r}\ge r/C(r)$
with $C(r)=\min(4,r)+3\min(3,r)+5\min(2,r)+8\min(1,r)$, hence $\delta_{4,2}\ge 2/26=1/13$,
$\delta_{4,3}\ge 1/10$ and $\delta_{4,4}\ge 4/31$. With Theorem 1 this gives
$1/13\le\delta_{4,2}\le 1/10$, the conjectured $9/80=0.1125$ lying strictly above the upper end.

**Proposition 3** states the baseline the paper insists on: every $(3,1)$ broadcast is a $(4,2)$
broadcast, so $\delta_{4,2}\le\delta_{3,1}\le 1/8$ already follows from the source's own theorem on
$\delta_{3,1}$ (Theorem 4.5 there, per the bibliography entry). At $r=2$ the paper therefore reads
its own contribution as the chain $1/8\to 9/80\to 1/10$.

No optimal density is claimed. §5 ("What is not settled") lists $\delta_{4,2}\in[1/13,1/10]$,
$\delta_{4,3}\in[1/10,1/8]$ and $\delta_{4,4}\in[4/31,5/36]$ as open, and notes that an optimal
broadcast need not be periodic at all.

## 2. What the program checks

`verify.output.txt` records one `PASS` line per check, no `FAIL` line, the closing
`VERDICT: ALL 99 CHECKS PASS`, and `program exited with status 0`. The program rebuilds the graph
from the six adjacency rules of §2, evaluates reception by a fresh breadth-first search at every
class representative, and uses `fractions.Fraction` (standard library only) throughout.

The statement the run bears on is **Theorem 1**, and 33 of the 99 checks carry it:

* the presentation of the graph and the shell sizes (2) — `graph.*` (5 checks: $3$-regularity on a
  $324$-vertex window, symmetry of adjacency, no loops, $\mathbb{Z}^2$-translation an automorphism,
  and the short-cycle profile at $36$ sampled vertices: no $3$-, $5$-, $6$- or $7$-cycle, exactly
  one $4$-cycle and exactly two $8$-cycles), `shells.*` (2 checks: shells $(1,3,5,8,11)$ from all
  four vertex types $a=0,1,2,3$, and the coordination prefix $3,5,8,11$) and `ball-sizes` (1 check:
  $|B_2|=1+3+5=9$, $|B_3|=1+3+5+8=17$);
* the three witnesses — `W5` (9 checks), `W6` (9), `W7` (7); the transcript identifies these with
  $T_5$, $T_6$ and $T_7$. For each: the class count ($4\times5=20$, $4\times4=16$, $4\times18=72$),
  invariance of membership under both stated lattice generators, the exact density ($1/10$, $1/8$,
  $5/36$), the broadcast property (min $f=2$, $3$, $4$ with $0$ violations over all classes), the
  $f$-range, reproduction by fresh BFS of *all* entries of the table printed in the paper, and the
  comparisons $1/10<9/80$, $1/8<1/7$, $5/36<1/6$. `W5` additionally records that $1/10<1/|B_2|=1/9$
  and that the class $0@(1,0)$ has no tower within distance $2$ and is served by two towers at
  distance exactly $3$; `W6` itemises its two tight classes $1@(0,0)$ and $2@(0,0)$, where $f=3=r$.

The other 66 checks are on objects the paper neither exhibits nor tabulates: further patterns
written out inside `verify.py` — `W1` (13 checks, density $9/80$), `W2` (6, $1/7$), `W3` (6, $1/6$),
`W4` (17, $1/9$) — plus `T31` (10) and `baseline` (2), the `lower` block (7) and the `record`
consistency block (5). The Verification section of the paper says as much: the verdict line covers
patterns and comparisons the paper does not use, so "that line is not by itself a certificate of the
statements above: only the groups named in the previous paragraph are", namely `graph`, `shells`,
`ball-sizes`, `W5`, `W6`, `W7`.

## 3. What the program does **not** check

* **Theorem 2 and Lemma 4 are hand proofs; the program is a control on their arithmetic only.**
  The `lower` block computes $C(2)=26$, $C(3)=30$, $C(4)=31$ and the fractions $2/26=1/13$,
  $3/30=1/10$, $4/31$; each of its three `...-arithmetic` lines carries the text
  `NOTE (capping/double-counting lemma ... not verified here)`, and the closing `NOT RE-RUN` block
  repeats that the capping/double-counting step behind $r/C(r)$ is not verified. The lower end of
  $1/13\le\delta_{4,2}\le 1/10$ therefore rests on §4 alone, and the paper's Verification section
  says exactly this.
* **Proposition 3 is not verified.** The `T31` block runs one $(3,1)$ broadcast of density $1/8$,
  written out in the program, through the $(3,1)$ predicate and then the $(4,2)$ predicate. That is
  one witness, not the general implication, which is proved by hand in §1. Nor does the program
  check that this pattern is the source's own; the transcript describes it only as "the $(3,1)$
  broadcast T31 written out in this file".
* **No lower bound and no optimality of any kind.** The `NOT RE-RUN` block states that the program
  does not search over period lattices, establishes no lower bound and no optimality, and that no
  claim that any density is least is checked. Nothing asserts $\delta_{4,2}=1/10$; the `record`
  lines repeat that each lower end rests on the capping lemma not verified there.
* **The ambient graph is not identified up to isomorphism.** The vertex-figure check is at $36$
  sampled vertices, and its own line records "no classification over Archimedean nets and no global
  isomorphism is checked here"; §2 likewise says "we do not verify that no other Archimedean net has
  that vertex figure". Everything else is checked on the presentation given by the six adjacency
  rules of §2 and should be read that way.
* **Which quantifiers are exhaustive and which are sampled.** The quantifier "every vertex of
  $H_{\infty,\infty}$" in the broadcast property is *not* sampled: equation (3) of §2 together with
  the lattice-invariance checks reduces it to the $20$, $16$ and $72$ classes, and those are
  evaluated in full. The sampled quantifiers in the run are the vertex figure ($36$ sampled
  vertices) and the $324$-vertex window carrying $3$-regularity, symmetry and absence of loops.
* **Numbers transcribed from the cited source.** The conjectured $9/80$, $1/7$, $1/6$; the source's
  $\delta_{3,1}\le 1/8$, $\delta_{3,2}\le 1/6$, $\delta_{3,3}\le 1/4$; the source line numbers and
  published labels (Conjecture 4.10, Problem 5, Theorems 4.5 and 4.6) — all are transcribed from the
  cited paper and none is verified by the program. The coordination sequence $3,5,8,11,\dots$ and
  the count $|P_v(4)|=1+3+5+8=17$ are also the source's integers (§2 attributes the first to the
  OEIS entry and the second to a figure caption there); the program recomputes them from its own
  presentation of the graph, which is a consistency check and not a reading of the source.
* **The source's figures are not read.** The `NOT RE-RUN` block states that the program reads no
  external file and no figure and makes no claim about the figures of any cited paper. §5 says those
  sets appear in the source only as raster images, were not transcribed, and that nothing whatever
  is asserted about them; so the literal form of the source's companion problem — a proof that *its
  own* three figures dominate $H_{\infty,\infty}$ — is untouched, as §5 states. No computation
  outside this folder is involved, and no claim in the paper rests on one.

## 4. How to check it

```sh
python3 verify.py            # prints 99 PASS lines and the verdict; exits 0 iff all pass
shasum -a 256 verify.py      # pair the program with the transcript header
```

Python 3 with the standard library only; the header of `verify.output.txt` records `Python 3.9.25`
and, beside the program name, the SHA-256 of the program that produced the output, so the two can be
paired. Computed here from the shipped `verify.py`:

    5a09a0dd5307f140bb1a3f9a7fecbd03a2729c56a49f339b6c9e5f29ce6183b9

which is the digest in that header. The tables of §3 are short enough to check by hand from the
congruences and the six adjacency rules of §2, independently of the program.
