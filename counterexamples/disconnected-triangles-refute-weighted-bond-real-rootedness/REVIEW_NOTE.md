# Disconnected Counterexamples to a Real-Rootedness Conjecture for Weighted Bond Posets

`disconnected-triangles-refute-weighted-bond-real-rootedness`

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
passes. The recorded run reports **32 checks, all passing**:

    VERDICT: ALL 32 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    3da8de4440e9a8a4aac4dc73eba495909121804444346982a2589d05e67f21d1

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: (i) the Moebius polynomials mu_{K_3} = 2 + 5t + 2t^2 and mu_{P_3} = 1 + 3t + t^2, and the multiplicativity mu_{rX} = mu_X^r, are INPUTS here, quoted by the paper from Equation (1.5) and Proposition 2.2 of the cited article; the definition of the weighted bond poset is in that article and is not reproduced in this folder, so no poset is enumerated and the identification of A and B with mu_{K_3} and mu_{P_3} is not checked. (ii) Also transcribed rather than checked: the numbering of Conjectures 4.13(2), 4.13(1) and 4.12(2), and the reading that 4.13(2) is stated there for arbitrary k with no connectivity hypothesis, which is the scope claim on which admitting disconnected graphs depends. (iii) The exact real-zero count with multiplicity is carried out for r = 1..12 only; for larger r the theorem rests on the paper's Q_r argument, whose ingredients are verified above in the parameter c -- the sum-of-squares identity in Z[c][t], Res_t(A - cB, B) = 1 for every c, the leading coefficient 5 - 4c, and no real zero at 47 rational c in (-1,1) -- while the factorization A^r - B^r = prod_j (A - zeta^j B) over C and the value c = cos(2 pi/r) are taken as read, except at r = 3, where c = -1/2 is rational and the factor is verified exactly. (iv) Two further truncations of quantifiers over r: the graph-family checks -- 3r vertices, exactly r components by union-find, rP_3 a proper spanning subgraph with one edge deleted per copy, and the even sign exponent 3r - r -- run for r = 1..8 only, and the gamma-expansion identity together with the positivity of its coefficients C(r,j)(2^{r-j} - 1) runs for r = 1..12 only, whereas the paper asserts the graph-family facts for every r >= 3 and gamma-positivity for every r >= 1; both constructions are uniform in r, but neither quantifier is exhausted here and no induction is machine-checked. (v) No trigonometric or complex arithmetic occurs anywhere in this program: that zeta = e^{2 pi i/r} is nonreal with |c| = |cos(2 pi/r)| < 1 for every r >= 3, that the leading coefficient |2 - zeta|^2 equals 5 - 4c, and that Q_r divides A^r - B^r are all taken as read -- the last verified only at r = 3 -- and the c-grid is 47 rationals k/24, which contains a genuine cos(2 pi/r) only in that same r = 3 case; the passage from the ingredients checked above to the paper's conclusions that Q_r(t) > 0 for every real t and hence that A^r - B^r fails to be real-rooted for every r >= 3 -- beyond the members r = 3..12 decided outright above -- is its own inference, printed here as a NOTE and not as a check. (vi) Nothing bibliographic or editorial is fetched or tested: the existence, authorship and contents of arXiv:2608.08692v1, the novelty or priority of the refutation, and the paper's negative scope claims -- that the connected restriction is not addressed, that no assertion is made about Conjecture 4.13(1), and that these examples do not refute the gamma-positivity Conjecture 4.12(2) -- are transcribed, as is the convention that 'real-rooted' means all 2r zeros real, which is the standard against which the counts above are compared. (vii) Checks named 'consistency' compare two printed quantities, as disclosed above.
