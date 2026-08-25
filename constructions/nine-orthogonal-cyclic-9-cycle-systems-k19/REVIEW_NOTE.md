# Nine Mutually Orthogonal Cyclic 9-Cycle Systems of K_19

`nine-orthogonal-cyclic-9-cycle-systems-k19`

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
passes. The recorded run reports **19 checks, all passing**:

    VERDICT: ALL 19 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    d209d90fff30183e1f30e18a08614ae40131cfcbb1e3aa1ff87c2bd9cd6a687b

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: (i) the report of eight mutually orthogonal cyclic 9-cycle systems attributed to earlier work, which lives in an external table; (ii) the upper bound mu'(l,n) <= n-3, quoted from the literature -- only its evaluation 19-3 = 16 is checked; (iii) whether some entirely different family of ten mutually orthogonal cyclic 9-cycle systems of K_19 exists -- that is a maximum-clique question on the 128304-vertex orthogonality graph and is out of budget; what is settled here is that THIS family of nine admits no tenth member; (iv) the census is exhaustive only modulo the structural lemma printed above: its arithmetic steps are verified exhaustively on Z_19 and its conclusion is confirmed against brute force only in the small analogue (p,l) = (7,3). No assumption-free enumeration of all 9-cycle systems of K_19 is attempted, and none of (i)-(iv) touches Theorem 1: the family of nine and its 36 orthogonal pairs are established directly from the definition, with no lemma and no census.
