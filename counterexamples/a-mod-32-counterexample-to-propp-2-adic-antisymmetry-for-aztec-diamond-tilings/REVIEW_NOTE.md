# A mod 32 counterexample to the congruence form of Propp's Conjecture 2 for Aztec diamond tilings by dominoes and square tetrominoes

`a-mod-32-counterexample-to-propp-2-adic-antisymmetry-for-aztec-diamond-tilings`

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
passes. The recorded run reports **63 checks, all passing**:

    VERDICT: ALL 63 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    3434223c3f9a6a04aef3fe9aeedcfd2cf0129532d07d5de40a43039b3f26de5c

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> scope
> NOT RE-RUN: M(14) and M(15) EXACTLY. Both engines and this program carry them modulo 2^8 only; the exact integers are near 10^41 and 10^47 and are not computed anywhere here. Proposition 2 needs only v_2 < 8, so the residues settle the refutation.
> NOT RE-RUN: M(n) for any n >= 16. Nothing beyond index 15 is computed, so every instance of Conjecture 2 with an index >= 16 -- including the two remaining k = 5 cells (13, 16) and (12, 17) -- is left open in both directions, and (14, 15, 5) is shown least only within the box n, n' <= 15.
> NOT RE-RUN: k >= 6. Any such test needs max(n, n') >= 31 (check F1), a 2^63-state computation, and none was attempted.
> NOT RE-RUN: Propp's Conjecture 1, the 2-adic continuity of M. It is neither used nor tested here, and the antisymmetry restatement Mhat(-3-n) = -Mhat(n) falls only conditionally on it.
> NOT RE-RUN: Propp's theorem that M(n) is odd. Check E9 verifies that the 16 residues in hand are odd and that every k = 1 cell holds, which is a control on the data, not a proof of the theorem.
> NOT RE-RUN: the two C engines that produced the paper's data. This program is a third, independent implementation; it agrees with them at n = 13, 14, 15 mod 256, but it does not inspect or execute them.
