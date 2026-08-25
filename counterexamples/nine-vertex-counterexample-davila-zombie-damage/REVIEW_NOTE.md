# A counterexample to Davila's claw-free zombie-damage conjecture

`nine-vertex-counterexample-davila-zombie-damage`

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

    e4a97b2c6fb66a395fe1fa30826a69346922bdbf89a2ce21906210fbe481d5ac

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the order-nine census itself -- the 261080 connected and 4494 claw-free class counts at order nine, the count of exactly 151 isomorphism classes of minimum-order counterexamples, and the paper's statement that all 151 of them have (dmg,zdmg)=(2,9) -- is beyond a pure-Python budget, and the order-nine enumeration the paper credits to an external generator is not reproduced here. Established above instead: no connected claw-free graph of order at most eight violates the conjecture, G_2 (order nine) does violate it, and any order-nine violator with dmg >= 2 has values exactly (2,9). The order-nine violators with dmg <= 1, if any exist, are not swept here.
> NOT RE-RUN: the paper's 11117 connected classes at order eight. The order-eight sweep above enumerates only the claw-free connected graphs, so it checks the 881 of the claw-free column but not that entry; no claim verified here depends on it.
