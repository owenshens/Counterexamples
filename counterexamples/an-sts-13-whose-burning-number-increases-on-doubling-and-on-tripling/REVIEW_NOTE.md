# A Steiner Triple System of Order 13 Whose Burning Number Increases on Doubling and on Tripling

`an-sts-13-whose-burning-number-increases-on-doubling-and-on-tripling`

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
passes. The recorded run reports **41 checks, all passing**:

    VERDICT: ALL 41 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    587a6553fc2b98e115b114fc451e9f23d7eb8de7609544b0484b2b7fb6467480

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the target e-print is not fetched or parsed here. The verbatim question, its
> NOT RE-RUN: byte offset, the theorem labels and the nine printed values of h are quoted in
> NOT RE-RUN: the paper from the source; only the h values are re-derived above.
> NOT RE-RUN: Constructions 2.15 and 2.17 of the Colbourn-Dinitz Handbook are not consulted.
> NOT RE-RUN: The doubling and tripling rules are implemented from the standard statements
> NOT RE-RUN: reproduced in the paper, and are validated only by (i) the resulting designs
> NOT RE-RUN: being Steiner triple systems and (ii) the control triple(AG(2,3)) = AG(3,3).
> NOT RE-RUN: No isomorphism census is performed. The remark that there are exactly two
> NOT RE-RUN: STS(13) is neither used nor checked, and the non-projective non-degenerate
> NOT RE-RUN: STS(15) that would give a second witness order are not enumerated.
> NOT RE-RUN: b(T) is bounded below only; its exact value, 6 or 7, is not determined here.
> NOT RE-RUN: The separate observation about the target's Corollary with a missing rho >= 9
> NOT RE-RUN: hypothesis is not tested; nothing above depends on that Corollary.
