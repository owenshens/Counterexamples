# The Harborth Constant of C_9 C_9

`the-harborth-constant-of-c9-plus-c9`

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
passes. The recorded run reports **61 checks, all passing**:

    VERDICT: ALL 61 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    088e81e998a1dae913ee32236c53b58bbe257efef66c53e4ece47172b470f734

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE NOT RE-RUN: 595185279859131 of the 811509569804714 completions of eq. (2) (26.657% of the census was covered here). The paper's upper bound is an exhaustive search over all 831 feasible canonical prefixes; at the rate measured above (2e+11 completions/s in pure Python) finishing it needs of the order of 2971 s more, single process, which exceeds this program's budget -- and the prefixes not yet touched are the less efficient ones, so that is a lower estimate. Run 'python3 verify.py --full' to reproduce the census in full; no --full run accompanies this transcript, so the time above is an estimate extrapolated from the covered part and not a measurement.
> NOTE what IS re-derived exactly at n = 9: the 16-element lower-bound set and its maximality, the stabilizer of order 108, the 992 canonical triples, the 831 feasible prefixes and the integer 811509569804714 of eq. (2); and the identical search engine settles n = 3, 5, 7 exhaustively.
> NOTE the paper's node count 879672298 is not reproduced: this program prunes with degree peeling plus a greedy-colouring clique bound rather than an exact clique search, so its node count (72315063 here) is a property of the implementation, not a mathematical claim.
> NOTE NOT CHECKED: the authors' own census code, and the complete independent re-run of it reported in the paper. This program is standard-library only and reads no external file; it re-derives the mathematics instead.
