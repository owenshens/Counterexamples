# Review note

**Paper.** *A 12-Vertex Tree Refuting the $GRM_{-2}$ Lower Bound for Trees of Maximum Degree Four*
(`paper.tex`, `paper.pdf`).

**Everything a referee has here.** `paper.tex`, `paper.pdf`, `verify.py`, `verify.output.txt`, and
this note. There is no data file: `verify.py` builds every graph it examines from edge lists written
inside itself, takes no arguments and reads no input.

## 1. What the paper claims

With $GRM_\lambda(G)=\sum_{uv\in E(G)}(\deg u+\lambda)(\deg v+\lambda)$, so that
$GRM_{-2}(T)=\sum_{uv\in E(T)}(\deg u-2)(\deg v-2)$, the paper refutes one sentence of one source.
Section 1, "The claim under test, and its locator", quotes that sentence character for character
from line 838 of the version-1 e-print of Bašić and Ilić, arXiv:2604.06044, and reads it against the
class fixed at line 797 as: every tree $T$ of order $n$ with $\Delta(T)=4$ and $n\equiv 0\pmod 4$
satisfies $GRM_{-2}(T)\ge-n$.

Theorem 1 (Section 2) exhibits a tree $W_1$ on $n=12$ vertices, its eleven edges listed in the
proof, with $\Delta(W_1)=4$ and $GRM_{-2}(W_1)=-13<-12$, so the quoted claim is false. Proposition 2
(Section 3) constructs, for every $k\ge3$, a tree $F(k)$ of order $n=4k+4$ with $\Delta=4$ and
$GRM_{-2}(F(k))=-(n+2)$, so the claim fails by $2$ at every $n\equiv0\pmod4$ with $n\ge16$; $F(3)$
is the 16-vertex tree $W_2$ whose edge list is printed after the proposition. Corollary 3 concludes
that the source's claimed extremal family $TT^4_{opt}(k)$ (line 845), whose two claimed parameter
sets both give $GRM_{-2}=-n$ through the source's own formula at line 795, is not the set of
extremal trees of order $4k+4$ for any $k\ge2$.

Both proofs are hand computations: an eleven-edge table for Theorem 1, and a degree count for
Proposition 2.

## 2. What the program checks

`verify.output.txt` is the recorded run of `verify.py`. Every check line reads `PASS`; the run ends
`VERDICT: ALL 70 CHECKS PASS` and `program exited with status 0`. The blocks, and the claim each
one serves:

* **4 forced positive controls** (`control-K13`, `control-K14`, `control-K15`, `control-P5`):
  $GRM_{-2}$ of $K_{1,3},K_{1,4},K_{1,5},P_5$ against $-3,-8,-15,0$, values written as literals in
  `verify.py`. These pin the definition and serve no claim of the paper.
* **12 checks on $W_1$ ($n=12$) — Theorem 1.** $W_1$ is a tree (12 vertices, 11 edges, no loop, no
  repeated edge, connected); degree sum $22=2\cdot11$ with sequence $[4,4,3,2,2,1,1,1,1,1,1,1]$;
  $\Delta=4$ attained at $0$ and $8$; the edge partition $m_{13}=1,m_{14}=6,m_{23}=2,m_{24}=2$ and
  counts $n_1=7,n_2=2,n_3=1,n_4=2$ printed after the proof; and $GRM_{-2}=-13$ obtained three ways
  forced to agree — the direct edge sum, the source's pre-substitution formula at line 795
  (equation (4) of the paper) and its equation at line 800 (equation (3)). Two further identities,
  $-(n+3)+S$ with $S=2$ and a counting identity, are checked but are not asserted in the paper. One
  check evaluates the source's own degree-count system at lines 785–790 in exact rational arithmetic
  and confirms it returns $W_1$'s actual $n_1,n_2,n_4,m_{14},m_{24},m_{33}$, so the counterexample
  sits inside the source's own parameterisation. The status line reports `l.838 VIOLATED by 1`.
* **12 checks on $W_2$ ($n=16$)** — the $k=3$ tree of Section 3, same battery: $m_{14}=9$,
  $GRM_{-2}=-18$, `l.838 VIOLATED by 2`.
* **12 checks on `W3` ($n=20$)**, same battery: $m_{14}=11$, $GRM_{-2}=-22$,
  `l.838 VIOLATED by 2`. Its edge list is in `verify.py` only; the paper does not print it, and
  Section 4 refers to it as "the tree $F(4)$ on $20$ vertices".
* **11 checks on `D8` ($n=8$)** — the double star discussed after the proof of Theorem 1, with
  $m_{14}=6$, $m_{44}=1$, $GRM_{-2}=-8$. This is an anti-control: the run reports `l.838 HOLDS`
  here, so the checker does not flag every order in the residue class.
* **5 forced negatives** — the source's earlier, weaker bounds must survive exactly where the paper
  says they do. $W_1$ respects line 822 ($\ge-(n+2)=-14$) and meets line 830 ($\ge-(n+1)=-13$) with
  equality, so $W_1$ refutes line 838 alone; $W_2$ violates line 830 at $n=16$ and meets line 822
  with equality; `W3` meets line 822 with equality. This is the content of the paragraph beginning
  "Two remarks fix the scope of Theorem 1" and of the sentence after the $W_2$ edge list. (The paper
  uses no numbered remark environment.)
* **12 checks on the family — Proposition 2, instantiated.** For each $k=3,\dots,12$ the constructed
  $F(k)$ is a tree of order $4k+4$ with $\Delta=4$, $n_3=1$ and $GRM_{-2}=-(n+2)$ (at $k=12$:
  $n=52$, $m_{14}=27$, $GRM_{-2}=-54$); plus two isomorphism checks, by AHU canonical form, that
  $F(3)\cong W_2$ and $F(4)\cong$ `W3`. These two are the only comparisons of graphs in the run.
* **2 checks on the source's claimed extremal parameters — Corollary 3.** The parameter sets at
  lines 843 and 844, put through the source's line-795 formula for $k=2,\dots,12$, give $-n$ in
  every case.

Arithmetic is exact throughout: integers, and `fractions.Fraction` for the lines 785–790 system. No
floating point, no randomness, no network.

## 3. What the program does not check

* **The main theorem is a hand proof and the program is a control.** Theorem 1 is a finite
  computation on eleven edges that a referee should redo from the table in its proof. `verify.py`
  re-derives the same sum from the same edge list and cross-checks it against the source's two
  formulas; it supplies no step the paper does not print. The same holds for Proposition 2 and
  Corollary 3.
* **No enumeration, hence no minimum and no minimiser.** The closing SCOPE block of
  `verify.output.txt` states that the run recomputes $GRM_{-2}$ and the degree and edge counts of
  the listed graphs from their edge lists, that it does *not* enumerate trees, and that it therefore
  "establishes no exact minimum at any order and identifies no minimisers". The paper claims no
  exact minimum either: Section 4, first bullet, gives upper bounds on
  $\min\{GRM_{-2}(T):|V(T)|=n,\ \Delta(T)=4\}$ only and claims no minimisers.
* **The quantifier of Proposition 2 is proved, only sampled by the program.** The proposition is
  stated for all $k\ge3$; the run instantiates $k=3,\dots,12$ and says so in its own words, "the
  general k >= 3 claim is NOT machine-checked here". The general statement rests on the printed
  proof.
* **$n=8$ is not decided.** The SCOPE block says the run "does not decide the case n = 8 beyond the
  single double star D8", and the paper likewise declines to decide whether $GRM_{-2}(T)\ge-8$ holds
  for every tree of order $8$ with $\Delta(T)=4$, adding that nothing here depends on it.
* **Nothing about the source is verified against the source.** The quoted line-838 sentence, the
  class at line 797, the formulas at lines 795 and 800, the system at lines 785–790, the bounds at
  lines 822 and 830, and the parameter sets at lines 843 and 844 are transcriptions from the cited
  e-print into `paper.tex` and `verify.py`. The program checks that those transcribed formulas agree
  with each other and with the exhibited trees; it cannot check that the transcription is faithful.
  Corollary 3 in particular is only as strong as the two transcribed parameter sets. A referee
  should read arXiv:2604.06044v1 against Section 1. No computation outside this folder is involved,
  and no claim of the paper rests on one.
* **Section 4 lists the paper's own further abstentions**, which the program does not address
  either: no corrected replacement for $TT^4_{opt}(k)$ and no repair of the structural
  characterisation at line 846; no claim about the source's $\Delta=3$ results, its $\lambda\ge-1$
  material, or its $\Delta=4$ conclusions at $n\equiv1,2,3\pmod4$.

## 4. How to check it

From this folder, with Python 3.9 or later:

```
shasum -a 256 verify.py
python3 verify.py
```

The first command must print

```
d9e8b947ed29f57d8326deda51cdcbf9780e2b8d95d2702cfbe4388feedc0a46
```

This is the digest on the `sha256:` line of the header at the top of `verify.output.txt`, which also
records the program's file name and the Python version used, so transcript and program can be
paired. The second command should reproduce `verify.output.txt` from the line beginning
`CHECKER for` onward — the four lines above `=== program output follows ===` are the header, not
program output — ending in `VERDICT: ALL 70 CHECKS PASS`, and exit 0.

To rebuild `paper.pdf`, run a LaTeX engine on `paper.tex`; it needs `amsart` with `amsmath`,
`amssymb`, `amsthm`, `array`, `hyperref`, `microtype`, `fontenc` and `lmodern`.
