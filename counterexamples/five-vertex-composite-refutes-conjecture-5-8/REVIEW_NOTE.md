# A Five-Vertex Counterexample to Conjecture 5.8 of Curto, Geneson, and Morrison

`five-vertex-composite-refutes-conjecture-5-8`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **42 checks, all passing**:

    VERDICT: ALL 42 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    fc71bceab7074b048ae1e7c9c143d8e3d85e57420e4371d3630d6c279ff8a4cd

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: no external catalogue or search is needed for this paper, and every census it reports (1-, 2- and 3-vertex skeletons, all corners of Table 1, the 12 and 15 determinant lists) is re-enumerated above. The reduction lemma's general form (the paper's Lemma 2), for components of unbounded order, is checked only on the finite family reported above: the largest component order exercised anywhere in this program is 3, and the exhibited counterexample's own components have order at most 2. Two consequences of that, spelled out because they are not visible in the check lines. (i) The minimality claim (Theorem 3: every counterexample has at least five vertices) is NOT fully machine-checked here. Its exclusion of two and three components quantifies over components of arbitrary order; what this program establishes is the censuses of permitted 1-, 2- and 3-vertex skeletons, the multiaffine corner argument giving A(d)^{-1}1 > 0 on the whole open box (alpha,beta)^n for every permitted 2- and 3-vertex skeleton (the one-component case being tautological), and the reduction lemma itself only for components of order <= 3. So minimality rests, for components of larger order, on the paper's hand proof of that lemma. (ii) The total-activity bound beta^{-1} < 1^T(I-W_i)^{-1}1 < alpha^{-1}, which the paper quotes from Lemma 3.1 of the cited work in order to get alpha < d_i < beta, is likewise only re-verified here for components of order <= 3; it is not proved in general and this program does not check it beyond that finite family. The forward implication of the conjecture (the cited Theorem 5.2) is also not re-proved, being irrelevant to the refutation.
