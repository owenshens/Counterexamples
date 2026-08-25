# Order-Minimal Pendant-Free Counterexamples to Conjecture 2 of Nurdin et al.

`order-seven-counterexamples-nurdin-irregularity-conjecture`

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
passes. The recorded run reports **28 checks, all passing**:

    VERDICT: ALL 28 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    aa2f88eca7e965b1bd745d8e7ee0716d3a0044ae9b58ce7514fbd86afb12da8e

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: Theorem 2 and Corollary 3 are infinite families, so only a finite window of each is verified above, and the window is narrower than the statements in three ways. (i) Theorem 2 is checked for the admissible (t,d) pairs with t <= 4 and d <= 12 only. (ii) The theorem asserts its conclusion for EVERY (d-t-2)-regular graph H on d vertices; every such core up to isomorphism is tested here only for d <= 7, where the regenerated census supplies them all, and for 8 <= d <= 12 a single circulant core per pair is tested instead. (iii) Corollary 3 is checked for orders 8..16 only, with the paper's own core H = complement(C_d); order 7 is covered by the census instead. For all pairs and orders except (t,d) = (1,4), (2,4), (1,5), tvs(G) > 2 comes from the computed endpoint-weight obstruction, not from exhaustion; on those three smallest members, and on every census graph where the obstruction fires, it is cross-validated against complete exhaustion. The order <= 7 classification (Theorem 1) is verified in full: the census is regenerated from scratch, not read from any catalogue.
