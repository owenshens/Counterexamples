# Counterexamples to the Zhang–Wan Deep-Hole Completeness Conjecture

`counterexamples-to-the-zhang-wan-deep-hole-completeness-conjecture`

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
passes. The recorded run reports **24 checks, all passing**:

    VERDICT: ALL 24 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    21fa80814bd3f7b3ba0e6d58560594576dd4e2780e975ed98d17f8bee9ed27d8

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: the endpoint theorem in the generality claimed -- it is asserted for every odd prime power q, every elliptic curve E/F_q with #E(F_q) >= q+4 and every choice of O, whereas this program derives it in full only for the single curve y^2 = x^3 + 2x + 1 over F_3 with O the point at infinity, plus a bounded sweep of curves over the PRIME fields F_3, F_5, F_7, F_11 (capped at two curves for q=11) with O at infinity: no non-prime q (F_9, F_25, ...) is touched at all, which is exactly the case the paper says the source's side conditions miss, and no origin O other than the point at infinity is used; that the projective model has exactly one point at infinity, namely O = [0:1:0], is assumed rather than derived -- N is formed as (number of affine solutions) + 1 -- and the paper's own particular singularity computations (F_X = Z^2 over F_3, and F_Z(O) = 1) are not reproduced step by step, being replaced by an exhaustive Jacobian scan of PG(2,q); the codimension-three lemma as a general statement about arbitrary [n,n-3] codes over F_q with d >= 3, and its PG(2,q) line-counting proof, are only instantiated here, never proved; the general theorems the paper invokes -- Riemann-Roch giving l(mO) = m in genus one, injectivity of evaluation when deg(kO) < n, the designed-distance bound d >= n-k, Hasse's bound, and the cited textbook -- are used as the monomial-basis construction and confirmed only by rank and weight computations on the finitely many codes built here; every statement about the source paper is unverified because no external document was fetched or parsed, including that its Conjecture 1.4, Theorem 1.2(ii)-(iv), Remark 1.3 and Corollary 4.11 say what is quoted, that their count (#E(F_q)-n)(q-1)q^k specialises to the layer, that their side conditions are n >= q+k or q prime or k <= sqrt(q) and that at k = n-3 the first forces q <= 3 while the third fails, that Corollary 4.11 concerns residue codes and does not conflict, and the remark about preprint versus printed numbering; the claim that at k = n-2 = N-3 the covering radius drops to 1 so completeness fails trivially is not computed; at k = 2 the paper's hand arguments -- that every word of F_3^6 agrees with a codeword in at least three coordinates, and the step A_0 + A_1 = F_3 -- are not re-derived, only their conclusions rho = 3, d(w,C) = 3, w outside C_L(D,3O) and the 24/22 census are recomputed, and only for this one F_3 curve; nothing is checked about the interior of the range, where the paper itself decides nothing for q > 3 and any k < N-4; on every curve of the bounded sweep only the single value k = N-4 is examined, so no k < N-4 -- in particular not k = 2 -- is examined for any curve other than y^2 = x^3 + 2x + 1 over F_3, and the sweep enumerates only models in short Weierstrass form y^2 = x^3 + a x + b, one instance per (a,b) pair with no reduction to isomorphism classes, so it covers neither every elliptic curve over those four fields nor every model of the curves it does reach; the prior-art and novelty statements -- no resolution of the conjecture in the literature, no earlier appearance of the exhibited word, no novelty claimed for the lemma -- involve no computation and no literature search, of the citing literature or otherwise, was performed; and every provenance and record-keeping assertion the paper makes about material outside itself is uncorroborated here: that a first implementation of the census exists, that a second one was written independently from the definitions of the note, that its 81 coset-weight determinations and the counts 24 and 22 agreed throughout with the printed values, and that both programs and the transcript of that run are retained in an auxiliary archive available on request, are statements this program cannot test -- no such archive, program or transcript was read by it, none is distributed beside it, and neither earlier program is re-run here; this program is a further, independent implementation written from the printed data alone, and its agreement with the printed numbers is the only evidence for the census offered in this bundle.
