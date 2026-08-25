# An Order-Three Counterexample to the Barát–Korondi–Varga Projection-Area Conjecture

`an-order-three-counterexample-to-the-barat-korondi-varga-projection-area-conjecture`

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

    2c6ec98080d7fdb5f0a449ba09974bd4079a533ddaaf4dec25359fe4dc3bbd85

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> INFO SCOPE: the order-three census above is COMPLETE -- all 4224 independent 9-subsets of [3]^3 were enumerated and each tested for solution-hood, so the theorem s_3(3) = 20 < 23 is fully re-derived here, independently of the Barat-Wanless catalogue. NOT RE-RUN: n >= 4. The conjecture is stated for all n and the paper claims nothing beyond n = 3; a census at n = 4 would have to test 15183071352 independent 16-subsets of [4]^3 (counted here by the transfer-matrix DP), about 3594477 times the order-three workload, which is many orders of magnitude past a 25-minute budget. ALSO NOT RE-RUN: the mutation suite described in this program's header. Those 31 corruptions were run while the program was being written, by editing copies of it; no mutation harness, list or log ships with it and nothing above re-executes them, so the note's mutation claim is NOT evidenced by this transcript. ALSO NOT REPRODUCED: the note's attribution claims, which are about the literature and not about mathematics -- that the inequality transcribed here is the one posed by Barat-Korondi-Varga (the note quotes it at second hand, as the formulation printed as Conjecture 1 of Barat-Wanless, and states that it cannot confirm what number it carries in the original; this program asserts no number for it), that the seven totals appear in Figure 6 of Barat-Wanless in the printed left-to-right order (only their MULTISET is checked above, which is all the theorem needs), and that B is isometric to the representative there of total projection area 20. Those require reading the cited papers.
