# The Filtered Integral Span at Weights 18 and 20:\\B(18)=13, B(20)=14, and a Suggested Closed Form

`the-bachmann-yu-filtered-span-closed-form-fails-at-weights-18-and-20`

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
passes. The recorded run reports **90 checks, all passing**:

    VERDICT: ALL 90 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    fe9f087d151babcaab2b995d698d4da238c5a4b96b324c2cc2c776d0c2022f65

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the following are OUTSIDE this program and are not asserted by it.
> NOT RE-RUN: (a) B(k) for k = 22 and beyond, and for odd k -- only k <= 20 is computed here;
> NOT RE-RUN: (b) any replacement closed form for B(k) -- nine data points are printed, no
> NOT RE-RUN:     formula is fitted or tested;
> NOT RE-RUN: (c) Bachmann-Yu mainconj:integralspan (R_g = M_Z) in general -- what is verified
> NOT RE-RUN:     is only its weight-18 and weight-20 instances, which follow from the attained
> NOT RE-RUN:     index 1 at n = 13 and n = 14;
> NOT RE-RUN: (d) any explicit element h of F_{<=18}M_Z \ Lambda_12 -- no witness element is
> NOT RE-RUN:     claimed or checked, because 31 or 259 printed coefficients cannot pin a
> NOT RE-RUN:     53-dimensional space; the deciding object is the index chain itself;
> NOT RE-RUN: (e) the literature: novelty, priority and the Craig near-miss are not machine
> NOT RE-RUN:     checkable and nothing here speaks to them;
> NOT RE-RUN: (f) the byte offsets quoted from the arXiv source -- this program reads no
> NOT RE-RUN:     external file and cannot confirm a line number in someone else's TeX.
