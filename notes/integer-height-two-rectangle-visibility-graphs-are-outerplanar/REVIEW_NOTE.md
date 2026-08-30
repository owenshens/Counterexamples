# Integer Height-Two Rectangle Visibility Graphs are Outerplanar

`integer-height-two-rectangle-visibility-graphs-are-outerplanar`

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
passes. The recorded run reports **45 checks, all passing**:

    VERDICT: ALL 45 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    200e8e2bbf37b84391b7590475594f9a63fd82dff8d533a92da9e43efb963dc7

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the enumeration is exhaustive for bounding boxes [0,w]x[0,2] with w <= 7 (711,202 layouts) and NOT beyond; the proof in the paper is what covers all w and all n, and this program checks that proof only on those widths.
> NOT RE-RUN: the brute-force cyclic-order route runs on w <= 4 only (1,638 layouts); for 5 <= w <= 7 outerplanarity is certified by the proof's own construction, not by an independent search.
> NOT RE-RUN: nothing here touches height(G) for h >= 3 beyond the two exhibited anti-control layouts, and nothing here verifies that height(K_{1,5}) equals 3 (only that K_{1,5} is realisable at height 3 and at no height 2).
> NOT RE-RUN: no literature claim is checked by this program. The prior-art status of the result -- including the reading of the 2014 GD poster abstract cited in the paper, which was obtained and searched by hand and not by this program -- is outside its scope.
