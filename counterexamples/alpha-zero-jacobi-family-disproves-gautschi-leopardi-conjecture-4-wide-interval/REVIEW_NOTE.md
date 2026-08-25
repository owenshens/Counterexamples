# Counterexamples to the Wide-Interval Case of Conjecture 4 of Gautschi and Leopardi

`alpha-zero-jacobi-family-disproves-gautschi-leopardi-conjecture-4-wide-interval`

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
passes. The recorded run reports **22 checks, all passing**:

    VERDICT: ALL 22 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    f438bcab7b69f431baf7b8b41f532cf3e2147c277d1a776ca48d7e5cd8fa6bec

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: the two hypergeometric representations quoted from the standard theory, P_n^{(0,beta)}(1-2z) = 2F1(-n, n+beta+1; 1; z) and its general-alpha analogue, are taken here as the definitions and are not re-derived from any other normalization of the Jacobi polynomials, so the triviality of the alpha = 0 normalization -- P_n^{(0,beta)}(1) = 1, which the proof uses to drop the tilde, and likewise Ptilde_n^{(alpha,beta)}(1) = 1 in the Remark -- is inherited from the z = 0 value of those series rather than tested; the elementary trigonometric identities of the Corollary (cos 3phi = 4c^3-3c, cos 2phi = 2c^2-1, and the monotonicity of cos that turns 7/8 > sqrt3/2 into 0 < phi < pi/6) are imported in the same way, corroborated numerically at beta = 2 to 1e-45 but not proved here; every bibliographic and literature-history claim is unchecked -- the passage quoted from Gautschi and Leopardi and its placement after their Conjecture 4, the fidelity to the source of the statements of Conjecture 4 and of Conjecture 3 themselves, including the identification of their (3.6) with the displayed comparison and of their (3.7) with the two intervals, the description of their verification grid (alpha = 0.9, 0.8, ..., -0.9 with beta = 20, 19, ..., 0 and then negative values down to -0.999) and hence the inference that (alpha,beta) = (0,2) is one of its nodes, the attribution of the theta-grid theta_nu = nu*pi/1001, nu = 1..1000, to their routine for 0 < theta < pi with N = 1000 (the arithmetic on that grid is reproduced here, its provenance is not), the claim that the curve described in that passage runs from (-1,0) to (1,-1) and so nowhere exceeds beta = 0 and that the node therefore lies well inside the region reported positive, the contents of the Matlab verification core of their Section 2.2 (the line setting th1 = acos(-1/(2a+1)) and the commented-out th1 = pi) together with the circumstantial explanation drawn from it and the statement that the script actually used for Conjecture 4 is not printed, Koumandos's theorem at (1/2,1/2), (1/2,-1/2) and (-1/2,1/2), the history of Conjecture 3 including its disproof on a lens-shaped region and the later revised and extended parameter domains, the assertion that none of those papers revisits the polynomial-value implication, the remark that recent work distinguishes this problem from the corresponding questions for Jacobi zeros, the priority claim that no counterexample to the implication has previously been recorded, and the reference list itself (author attributions, titles, journal volumes, pages, DOIs and the arXiv identifier); also unchecked: the continuity step of Theorem 1 and of the Remark, which turns the endpoint sign Delta_2(pi) < 0 into the existence of epsilon_beta > 0 and into failure near theta = pi (the endpoint identity and its sign are verified for all beta, but the punctured neighbourhood is exhibited only at beta = 2, by the Corollary's interior witness and the 1000-node grid), the elementary range statement 0 < u < 1 for 0 < theta < pi (verified only at the 1000 grid nodes), the well-definedness of Theta_1^{(alpha,beta)} in (0,pi) as the point where P_1^{(alpha,beta)}(cos Theta_1) = 0 and of cos Theta_n^{(alpha,beta)} as the largest zero of P_n^{(alpha,beta)}, the Remark's general-alpha endpoint value, which is computed only at the four exhibited rational pairs (-1/2,1), (1/4,3/2), (1/2,2), (1,3) and not as an identity in (alpha,beta), so that neither the inference that those pairs exhibit a two-dimensional region of failure near theta = pi nor the non-minimality of the alpha = 0 family is established here, the behaviour of Delta_2 between consecutive grid nodes -- the count of 85 violators, the contiguity of the block nu = 916, ..., 1000 and the singleness of the crossing near theta = 2.8725 are facts about the 1000 samples and about one bracketed root, not exhaustive statements on (0,pi) -- any assertion about n >= 4, about the 0 < theta < Theta_1 branch of Conjecture 4, about Conjecture 5, or about zeros of Jacobi polynomials, and the sweep itself at any beta other than 2 or on any grid other than theta_nu = nu*pi/1001, nu = 1..1000.
