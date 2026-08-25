# Counterexamples to Gyoda's Injectivity Conjecture for Generalized Markov Trees

`counterexamples-to-gyoda-injectivity-conjecture-for-generalized-markov-trees`

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

    60865f7b4608afbac1bc1e6b2f077f91eb5570d65ebc2180122c246283ae7d18

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: (1) the infinitude of the family -- q_m integrality and monotonicity were re-derived only for 0 <= m <= 40, and the full mutation of MT(q_m,0,0,id) from its root only for 0 <= m <= 10 (the paper's own independently recomputed range), so the statement for all m >= 0, which rests on the induction and the mod-10 period argument, is not machine-checked; (2) the spine lemma as a theorem -- n_{1/d}=F_{2d+1} and the alternating tags were mutated out only for q in (0, 1, 2, 3, 6, 7, 12, 35, 100) and d <= 30, and n_{2/3}=29+10q only for 0 <= q <= 40, with the Fibonacci identity instance checked for d <= 60 and the mod-10 period for r <= 700; the general proofs by induction and Catalan's identity in general are not verified; (3) tree-wide facts for MT(22,0,5,id) -- exactness of mutations, tag permutations, the equation at every node, and the scans for repeated values or repeated (n_t,i_t) pairs cover only the full binary tree to depth 7, not the whole infinite tree, and therefore establish nothing about collisions deeper than that; and the Status section's algebraic remark that the equation is symmetric in x,y,z only when k1=k2=k3 is exercised here only at the single triple (22,0,5), through the scan of the six orderings of {2,15,889}, and is not verified for general k; (4) everything external to the paper's own arithmetic -- that the quoted conjecture is Conjecture 7.6 of the cited reference and is stated there as quoted, that the cited reference records that no counterexample was known, the definitions of the generalized Markov tree, of the root, of the mutation rules and of the synchronization with the Farey tree (all taken here from the paper's own recollection of them and not from the sources), the relation to Question 7.3 and to the separate Chen-Jia Conjecture 8.2, the posting dates of the cited preprints, the literature search, and the novelty/priority claim that no earlier refutation or earlier appearance of the collision at 889 exists, and the paper's account of its own provenance -- that the printed data were recomputed once independently, on an earlier occasion, with every check agreeing -- which this program can neither confirm nor deny, though the recomputation the paper describes there (the six checks of Section 2, and for each m with 0 <= m <= 10 the integrality of q_m together with the two labels obtained by mutating MT(q_m,0,0,id) from its root) is carried out among the checks above, so what is unverifiable is the history of that earlier run and not the facts it is said to have confirmed; (5) claims the paper explicitly does not make and this program does not supply -- global minimality of the counterexample, any statement about the classical case (0,0,0) or the symmetric locus k1=k2=k3, trees with sigma other than the identity, injectivity of t -> (n_t,i_t) as a theorem, and infinitely many collisions inside one fixed tree.
