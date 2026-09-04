# Referee note

**The Maximum Euler--Sombor Index at Diameter Three, and Both Extrema at Nine Vertices**

Files here: `paper.tex`, `paper.pdf` (6 sections; Theorems 1 and 2; Remarks 3, 4, 5; equations
(1)--(4)); `census.c` (enumeration kernel, modes `classes` and `labelled`) and `census.py` (driver:
exact 80-significant-digit `Decimal` arithmetic, the controls, the isomorphism and automorphism tests,
the verdict), with the recorded run `census.output.txt` (**43 checks, 43 PASS, 0 FAIL**, exited 0); and
`verify.py`, which enumerates nothing, with its recorded run `verify.output.txt` (**ALL 92 CHECKS
PASS**, exited 0, 60-significant-digit `Decimal`). `verify.output.txt` calls the paper's Theorem 1
"Theorem B" and its Theorem 2 "Theorem A"; `census.output.txt` also says "Theorem A" for the latter.

## 1. What the paper claims

With `EU(G) = sum over edges xy of sqrt(d(x)^2 + d(y)^2 + d(x)d(y))` (equation (1)) and `G_{n,D}` the
connected graphs of order `n` and diameter **exactly** `D`:

* **Theorem 1.** For every `n >= 4`,
  `M(n,3) = sqrt3*(n-2)^2*(n-3)/2 + sqrt(n^2-3n+3) + (n-3)*sqrt(3n^2-15n+19)`, attained by exactly one
  graph up to isomorphism, namely `H(1,n-3)`. §4 proves this in closed form and, in the words of the
  abstract, "uses no computation". At `n = 9` it reads
  `147*sqrt3 + sqrt57 + 6*sqrt127 = 329.777869165403581...`, margin `6.5941623645...` over `H(2,5)`.
* **Theorem 2.** `m(9,3) = 8*sqrt3 + 6*sqrt19 = 40.009800121795059...`, attained by exactly one class,
  the **bicyclic** theta graph `Theta(3,3,4)` of §3; so the minimiser of that cell is neither a tree nor
  unicyclic, and the best tree `D(3,4) = sqrt61 + 3*sqrt21 + 4*sqrt31 = 43.8290342121...` is only third,
  behind the unicyclic runner-up `4*sqrt21 + 8*sqrt3 + 4*sqrt7 = 42.7697144846...`.

The source is the open problem of Khanra and Das, item (1) of their concluding list ("the class of all
connected graphs with given diameter"), reproduced verbatim as **Problem 1** by Sekar, Balachandran,
Elumalai and Liu and posed with the order fixed as **Problem 3.1.6** of the 2026 survey of Ali, Gutman,
Réti, Albalahi and Hamza, which is the reading the paper adopts. One column and one cell are settled,
not the problem: the "Status" paragraph of §2 and all of §6 say that Problem 1 is **not** resolved, that
the minimum half at `D = 3` is undecided already for `n = 10`, and that no cell with `D >= 4` and
`n > 9` is touched. §2 makes no priority claim: it argues by an exact first-Zagreb computation that the
maximiser of the corollary in the paywalled paper of Vetrík (Discrete Appl. Math. **333** (2023) 59--70)
must be `H(1,n-3)` as well, and §6 says that if that paper's hypothesis does admit a concave
1-homogeneous edge function then Theorem 1 is a rediscovery, leaving Theorem 2's minimum half and the
two numerical values.

## 2. What the programs check

**`census.output.txt`: 43 checks, 31 in pass 1 and 12 in pass 2.** This is what carries the
*exhaustiveness* and the *uniqueness* of both extrema of `G_{9,3}` — the part of Theorem 2 that no hand
argument in the paper reaches — with Theorem 1's conclusion at `n = 9` as a by-product. The two passes
share no enumeration machinery, and the driver reproduces its controls before printing any extremum.

* *Pass 1 controls (9) and arithmetic safety (2): §5's counting claims.* graph6 input read class by
  class — `261080 = A001349(9)` classes, none malformed or disconnected, each in exactly one bucket; the
  eight measured diameter strata `1, 91518, 148229, 19320, 1818, 180, 13, 1` reproduce §5's table, so
  the cell holds exactly `148229` classes, as do §5's 21 per-size counts for `m = 8..28`; the cell
  reduces losslessly to `30966` distinct exact values whose smallest gap is
  `1.038011241280298561547138384E-7`, so at 80 digits no comparison here is delicate.
* *Minimum (5) and maximum (7): Theorem 2, and Theorem 1 at `n = 9`.* The measured extrema equal
  `8*sqrt3 + 6*sqrt19` and `147*sqrt3 + sqrt57 + 6*sqrt127` to `1e-60`, truncate to the digits the paper
  prints, and each is attained by **exactly one** class, zero ties: the argmin isomorphic to
  `Theta(3,3,4)` (`m = 10`, bicyclic, diameter 3, `|Aut| = 4`, margin `2.759914362838`), the argmax to
  `H(1,6)` (`m = 28`, diameter 3, `|Aut| = 720`, runner-up `H(2,5)`, margin `6.5941623645`).
* *Four further values (6): Theorem 2's ordering and §5's two pins.* Best tree of the cell
  `43.8290342121` and unicyclic minimum `42.7697144846`, giving the strict order bicyclic < unicyclic <
  best tree; the unicyclic *maximum* of the cell `62.7871695845`, equal to the right-hand side of
  Theorem 3.2 of Sekar et al. at `n = 9`; and the minimum of `G_{9,4}`, `18*sqrt3 = 31.1769145362`
  (Lemma 2 of Kizilirmak), attained by the 9-cycle alone.
* *Two negative controls (2).* The minima of `G_{9,6}` and `G_{9,7}` are each reported as attained by
  two classes, so the counter that says "exactly one" at both ends of `G_{9,3}` can say more than one.
* *Pass 2 (12): uniqueness again with no isomorphism machinery.* No generator, no canonical form, no
  isomorphism test: all `2^36 = 68719476736` labelled graphs are swept, `66296291072` connected,
  matching the exponential formula in exact integers; **nothing** in the cell lies below the minimum or
  above the maximum; exactly `90720 = 9!/4` labelled graphs lie within `1e-6` of the minimum and
  `504 = 9!/720` within `1e-6` of the maximum, each window holding a single edge profile — one orbit at
  each end and nothing else. Cross-pass checks bound `classes <= labelled <= 9!*classes` per stratum.

**`verify.output.txt`: 92 checks.**

* *Steps 1--2 (12 + 12): the two objects of §3,* rebuilt from the vertex-named edge lists printed
  there — order, size, degree sequence, edge profile, `d(u,v) = 3` with `u`, `v` non-adjacent, diameter
  exactly 3, the extremal pairs `d(a1,b2) = d(a1,c2) = 3`, bicyclicity of `Theta(3,3,4)`, and each index
  against its surd form (to `1e-50`) and against the summands §3 derives by hand.
* *Steps 3--4 (9 + 10): Theorem 1.* The ladder `H(1,6) > H(2,5) > H(3,4)` at `n = 9` with each apex part
  in surd form; the closed form against a direct evaluation of `EU(H(1,n-3))` for `n = 4..14`, where
  `a = 1` is also the unique maximising split; the polynomial identity and positivity behind the strict
  convexity of `phi`; inequality (3) for `n = 4..199` with slope `0.1339746`; and `n = 4`, where the
  only member of `G_{4,3}` is `P4`.
* *Steps 5--7 (15 + 4 + 7): §5's hand arguments.* Three double stars and one hand-built unicyclic member
  of the cell, giving `40.0098... < 42.7697... < 43.8290...`; the cheapest 10-edge profile with a single
  `(3,3)`-edge, which would have beaten the minimum (`39.9522562726`), against the diameter 4 of
  `Theta(1,2,7)`, `Theta(1,3,6)`, `Theta(1,4,5)`, and the next degree sequence `(4,2^8)` at
  `41.9506201793`; Remark 3's four theta graphs, equal in profile and index and separated only by
  diameter; the floor (4) at `m = 11`, `n = 9`, giving `46.5729217146`; the Step 3 tier cap
  `189*sqrt3 = 327.3576026305`; and `K_{1,8}`, of diameter 2 with `Delta = 8`, which is why Step 3 needs
  `diam >= 3`.
* *Steps 8--10 (5 + 10 + 4) and `[ANTI]` (4).* §5's strata identities re-added; witnesses attaining the
  two pinning values, with `C_9`, `K_9 - e` and `P_9` placed outside the cell; §2's `M_1` computation
  and the failure of the prior art's `f`-condition,
  `sqrt7 - sqrt3 = 0.913701 > sqrt39 - sqrt31 = 0.677234`. The four `[ANTI]` items are objects the
  checker must reject: `H(1,6)` plus the edge `uv` (diameter 2), `Theta(2,4,4)` (diameter 4, though its
  index equals the minimum exactly), a disconnected graph (diameter reported `-1`), and a deliberately
  wrong candidate maximum of 320.

## 3. What the programs do **not** check

* **Theorem 1 is a hand proof and `verify.py` is a control.** §4 proves it for every `n >= 4`; the
  program only samples that quantifier — closed form and unique split for `n = 4..14`, inequality (3)
  for `n = 4..199`, the convexity identity and its positivity on the finite box `2 <= c <= 29`,
  `0 <= t <= 29`. Neither program proves the theorem.
* **`verify.py` cannot see the exhaustiveness or the uniqueness of the minimum,** and its closing block
  says so: it enumerates nothing over `G_{9,3}`, and the claim that all `148229` classes were evaluated
  with exactly one attaining each extremum rests on the census, which it does not re-run; its Step 8
  strata counts are therefore **transcriptions** of the census, not independent re-reads. §6 states the
  same division of labour.
* **The two objects are transcribed by hand from §3** into `verify.py`: its Step 1 header records that
  the transcription is not machine-verified and that maximality is not re-checked there, its Step 2
  header that minimality and uniqueness are not re-run there.
* **Pass 2 of the census decides in double precision inside a `1e-6` window.** `census.c` tabulates the
  edge values as C `double`s and reports double-precision extrema and the two window populations; every
  exact 80-digit comparison, including the surd identities and the tie counts, is made by `census.py` in
  pass 1. The margins `2.759914362838` and `6.594162364481` are far outside that window, but pass 2
  alone would not separate values closer than it.
* **The two pinning values are surd expressions taken from the cited statements,** not recomputed from
  their proofs; what is checked is that a graph the program builds attains them. Theorem 3.2 of Sekar et
  al. and Lemma 2 of Kizilirmak are not re-proved.
* **Neither program reads Vetrík (2023),** whose full text §6 says was not accessible. The `f`-class
  exclusion is an exact computation on the condition **as printed in a second paper**, not a reading of
  Vetrík's hypothesis; the transcript's closing `SCOPE` line and §6 both say so.
* **Parts of §5's `m = 10` argument are checked by no program:** the diameter claims about the
  `C_a`--`C_b` dumbbells, and the classification of a connected graph with degrees in `{2,3}` and
  exactly two vertices of degree 3 as a theta graph or a dumbbell. §5 says itself that the remaining
  `m = 10` degree sequences, and all of `m = 9`, are settled by the census rather than by hand.
* **Also not re-run** (transcript `SCOPE`): any cell with `n > 9`, and the minimum half at `D = 3` for
  `n >= 10`, which the paper does not claim.
* Pass 1 consumes graph6 input from `nauty`'s `geng`, the public generator named in §5; it is not part
  of this folder and must be obtained separately. Pass 2 needs only the compiled kernel.
* Three items are labelled as if they matched paper text the shipped paper does not contain: Step 8's
  two counting identities (`12*2^21 = 25165824`, and `C(7,2) = 21` free pairs); the `47 = A000055(9)`
  unlabelled-tree control, a control on the generator rather than a claim of §5; and the census's
  `G_{9,6}`/`G_{9,7}` ties, called "the paper's two genuine ties" although the paper discusses neither
  cell. No claim of Theorem 1 or Theorem 2 depends on any of them.

## 4. How to check it

```sh
python3 verify.py                             # standard library only, no data file, no network
gcc -O3 -march=native -fopenmp -o census census.c -lm
python3 census.py <path-to-geng> ./census     # optional third argument: --skip-labelled
```

`verify.py` exits 0 if and only if every check passes; `census.py` exits non-zero if any check fails,
and on a failed control stops with `VERDICT: CONTROL FAILED -- no extremum is asserted by this run`
before any extremum is printed. The header of `census.output.txt` records the build line and the
reproduction command, and the transcript records the two pass-level commands, one piping `geng -q -c 9`
into `./census classes`, the other
`./census labelled 40.00980012179506 329.77786916540356 1e-06`; `--skip-labelled` omits pass 2.

Each transcript opens with a header carrying the SHA-256 of the file or files that were run, so program
and transcript can be paired. Recomputed here from the shipped files, and equal to the header values:

```
9a9a0ef5dad23b0739205259c2882af6138ea62b92d421f79cf1a4962cceba38  verify.py
8dc2c88908dc62af38dd046f82626824f88eeb106d7b3c8bdf964e398d0f684e  census.c
a9ba702473f35f7c41d1f77169e8e1e39ad7582e5eda2c144ace5c8067b8e01d  census.py
```
