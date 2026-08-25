# A Counterexample to a Rainbow-Domination Implication of Žerovnik

`a-counterexample-to-a-rainbow-domination-implication-of-zerovnik`

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
passes. The recorded run reports **34 checks, all passing**:

    VERDICT: ALL 34 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    2b3eba795fa097dcf60590d3ae4a39fd9bf80eac840a64b27efa38f7cd625219

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: (1) the three double-counting premises behind the lower bounds -- tau >= k(30-c) at the empty vertices, tau <= 3w-2e, and the cut bound 3c-2e <= 3(30-c) -- are transcribed from the paper as premises, and only their integer consequences are scanned exhaustively; they are not re-derived symbolically here. (2) No solver and no exhaustive search over all labellings was run, so the paper's report that its 0/1 program was computed to optimality with optima 16, 21 and 25 is not reproduced by an independent optimisation: the upper bounds here come from the three exhibited labellings and the lower bounds from those transcribed premises plus finite enumeration. (3) The cited cubic bound gamma_rk(G) >= ceil(k|V(G)|/6) for all cubic G and all k <= 6, and the cited equality lemma in its published generality, are not verified beyond the order-30, k=3,4,5 instances used here. (4) Of the identification as the Tutte 8-cage, only girth 8 is verified; the quoted uniqueness of the cubic graph of girth 8 on 30 vertices is not. (5) The self-duality of the duad-syntheme configuration is not verified; its only use, the perfect-class statement on the syntheme side, is verified directly instead. (6) No bibliographic or attribution claim is checked: the wording of the question being refuted and the |V|=2n reading of it, the cited remark that this graph is not 3-rainbow domination regular, the claim that no earlier work determines gamma_r4 of this graph, and the paper's own acknowledgement that its literature search was not exhaustive. (7) The paper's account of its own process is not checked: that the 93 elementary conditions were done by hand and independently re-done, and that the machine cross-checks the paper reports were carried out there as described. This program is an independent re-derivation of those cross-checks, not a transcript of the paper's own run, and no such transcript accompanies it. (8) Nothing is checked about graphs other than this one, including the second implication of the question for other cubic graphs and any characterisation of the cubic graphs attaining the bound.
