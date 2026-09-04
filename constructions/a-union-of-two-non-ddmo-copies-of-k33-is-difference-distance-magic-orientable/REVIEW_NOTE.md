# Review note

Paper: *A Union of Two Non-DDMO Copies of $K_{3,3}$ Is Difference Distance Magic Orientable*

Files a referee has here, and nothing else:

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper (5 sections, 5 pages) |
| `verify.py` | the verification program named in Section 5 |
| `verify.output.txt` | the recorded run of that program |

## 1. What the paper claims

The source paper of Aceska, Arcila-Maya, Carlson, Marr, Parnes, Ryan, Schuerger and Vasquez
(reference [1], arXiv:2601.20492v1) closes with a paragraph of open questions; its third
question is unnumbered prose in the Conclusion (Section 13) of that version and asks, verbatim,
"is it possible for a union of non-DDMOGs to be DDMO?". The paper answers **yes**.

Section 1 observes that the question unions oriented objects but applies an undirected
predicate, so it has two inequivalent readings (A) and (B), and that an instance of (B) yields
one of (A). **Theorem 1** is the whole claim: reading (B), and hence (A), holds with $k=2$ —
two disjoint copies of $K_{3,3}$ form a DDMO graph on 12 vertices, although $K_{3,3}$ is not
DDMO.

Two statements carry it. **Proposition 2** (Section 2) exhibits the witness: $\vec D$ is the
printed orientation of $K_{3,3}$ on $A=\{0,1,2\}$, $B=\{3,4,5\}$ with the nine arcs
$0\to3$, $4\to0$, $5\to0$, $3\to1$, $1\to4$, $1\to5$, $3\to2$, $2\to4$, $2\to5$, and the
labelings $(6,1,5,10,2,8)$, $(12,3,9,11,4,7)$ of the two copies form a DDM labeling of
$\vec D\sqcup\vec D$. Equation (1) of Section 2 collapses the six weight equations of $\vec D$
to $f_0=f_1+f_2$ and $f_3=f_4+f_5$, so the proof is four additions: $1+5=6$, $2+8=10$,
$3+9=12$, $4+7=11$. **Lemma 3** (Section 3) says no orientation of $K_{3,3}$ has a DDM
labeling; the paper states plainly that this is not new — it is [2, Theorem 5] at $m=n=3$, and
the parity proof given is the necessity half of that argument, reproduced only for
self-containedness. Section 4 adds one observation: the witness lies outside the hypothesis of
Theorem 5.1 of [1], since here no component is a DDMOG at all.

The paper says explicitly that **no minimality of the order 12 is claimed**.

## 2. What the program checks

`verify.output.txt` closes with

    VERDICT: ALL 49 CHECKS PASS
    === program exited with status 0 ===

49 `PASS` lines, no failures, plus two informational `NOTE` lines. Mapped onto the paper:

* **Part 1 — 8 checks.** Five bear on Section 2: the nine printed arcs are an orientation (no
  repeated and no opposite arcs); the underlying graph is $K_{3,3}$ on the printed parts, all
  nine $A\times B$ pairs occurring exactly once; $\vec D$ is 3-regular; its in-degrees are
  $(2,1,1,1,2,2)$ and out-degrees $(1,2,2,2,1,1)$ as printed; and every row of the weight
  matrix is $\pm(f_0-f_1-f_2)$ or $\pm(f_3-f_4-f_5)$, which is equation (1). The other three
  checks concern an order-5 oriented graph $W$ that the paper does not print (see §3).
* **Part 2 — 10 checks: Proposition 2.** The two label sets $[1,2,5,6,8,10]$ and
  $[3,4,7,9,11,12]$ are disjoint with union $\{1,\dots,12\}$; all six weights of each copy
  vanish; the disjoint union is an oriented graph on 12 vertices with all twelve weights zero
  and $f$ a bijection onto $\{1,\dots,12\}$; the four printed sum triples are correct and
  partition $\{1,\dots,12\}$; and every part sum is even ($12,20,24,22$, copy sums 32 and 46),
  which is the closing sentence of Section 3.
* **Part 4 — 7 checks, of which 4 are Lemma 3.** $K_{3,3}$ has 512 orientations; none admits a
  DDM labeling, by brute force over all $512\times 6!=368640$ bijections; independently, an
  exact `Fraction`-RREF kernel enumeration agrees with the brute force on all 512; and every
  $B$-vertex has odd outdegree minus indegree, which is exactly the parity obstruction of the
  paper's proof against the odd sum $1+\dots+6=21$. The remaining three checks concern $W_4$
  (see §3).
* **Part 6 — 1 of 3 checks: Section 4.** Theorem 5.1 of [1] is not an if-and-only-if, witnessed
  by $\vec D\sqcup\vec D$: its hypothesis fails while its conclusion holds.

## 3. What the program does not check

**The paper's answer is a hand proof and the program is a control.** Section 5 says so in those
terms: the 12-vertex witness is four additions and Lemma 3 is a parity count, and nothing in
the paper depends on a machine. The program re-reads the printed objects and re-derives what is
asserted about them.

**Lemma 3 is transcribed, not established as new.** The general statement of [2, Theorem 5]
($K_{m,n}$ with $m,n\ge3$ is orientable to a DDMOG iff $(m+n)(m+n+1)=4s$) is quoted from that
paper, not reproved; only the instance $m=n=3$ is proved by hand and swept by Part 4.

**Much of the run corresponds to no claim in the paper, and the paper says so.** Part 3 (5
checks: an 11-vertex zero-weight union $W\sqcup\vec D$), Part 5 (16 checks and one `NOTE`:
sweeps at orders 4 and 5, the order-5 DDMOG census, and the $5+5$ splits of $\{1,\dots,10\}$),
the three $W$ checks of Part 1, the three $W_4$ checks of Part 4, and the last two checks of
Part 6 (the arithmetic equivalence "$3k(6k+1)$ even iff $k$ even" for $k=1,\dots,24$; that
$k=2$ is attained and $k=1$ impossible) all test material the paper never asserts: it prints no
11-vertex witness, no $W_4$, no census, no witness count. Section 5 states this and adds that
where the transcript names a lemma or a remark of the paper for a step it did not re-run, no
such numbered statement appears in the paper — item (4) of the run's closing `NOTE SCOPE`
speaks of "the k-even remark", and the paper contains no such remark; it numbers exactly
Theorem 1, Proposition 2 and Lemma 3. A referee should read those blocks as surplus, not as
support.

**Carried over from the run's own closing `NOTE SCOPE` block, which lists what it does not
re-run.** (1) The order-6 census of [2] — their 22 unique 6-vertex oriented DDMOGs, and the 19
isomorphism classes of 6-vertex graphs of minimum degree $\ge3$ — is not enumerated and no
agreement with it is asserted; the only order-6 graph used is $K_{3,3}$, whose 512 orientations
are swept exhaustively in Part 4. (2) Component orders 7 and above are never enumerated. (3) No
lane containing an isolated (order-1) component is enumerated, and no minimum order and no
witness count is established. (4) The sufficiency half of the $k$-even statement rests on a
classical triple-packing construction of Guy (1976) and is not verified; only its necessity
half and the case $k=2$. (5) The order-4 check bounds labels by 12 and establishes nothing for
unbounded labels — the transcript flags it in-line as a bounded statement. Two Part-5 checks
are likewise flagged in-line as weaker than they might read: 18 of the 1792 order-5
orientations "realise some label set inside $\{1,\dots,12\}$", which is not the same as being a
DDMOG, and the 29 five-subsets used by the $5+5$ join are a subtotal of 99.

**Still open, by the paper's own statement (Section 5).** Which orders admit witnesses; how
many witnesses there are; whether one exists below order 12; the 7- and 8-vertex counts of [2],
which play no role; and the other two open questions of the same paragraph of [1] — which
graphs are DDMO, and what distinguishes DDMOGs from non-DDMOGs.

## 4. How to check it

```sh
python3 verify.py
shasum -a 256 verify.py
```

Python 3.9 or later, standard library only, no external data file; exact integer and `Fraction`
arithmetic throughout, no floating point. One `PASS` line per check; exit status 0 only if every
check passed. The second command prints

    3843314e74fc5216dbd75b60c476982cce85b4e537f39f97770a7494605e5a3c

which is the `sha256` recorded in the header of `verify.output.txt` beside the program name, so
the transcript and the program can be paired. That header also records the interpreter as
Python 3.9.25. The program is deterministic, so a re-run should reproduce the transcript apart
from the version string it echoes.
