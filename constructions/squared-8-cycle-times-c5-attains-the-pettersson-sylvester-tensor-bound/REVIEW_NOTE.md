# An Ordered Pair of Graphs of Nonzero Twin-Width Attaining the Pettersson–Sylvester Tensor-Product Bound

`squared-8-cycle-times-c5-attains-the-pettersson-sylvester-tensor-bound`

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
passes. The recorded run reports **58 checks, all passing**:

    VERDICT: ALL 58 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    8b9574f9f180cfebe774455c313be8ddcd880119bbd9934864ef1f4b526e1332

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: no EXHAUSTIVE search for a contraction sequence of the 40-vertex product. STEP 3b builds one greedily and certifies it at width 8, which proves tww <= 8; it does not prove that no width-7 sequence exists. That direction is not needed and is not claimed here -- it is Ahn et al. Lemma 3.1, via b_min = 8, that rules width 7 out.
> NOT RE-RUN: no MINIMALITY. Nothing here bears on whether 40 is the least order of a tight product, nor on whether a tight pair with non-regular H exists. The only completeness statement made anywhere is STEP 6, the 121 ordered pairs of non-cographs on at most 5 vertices; the witness lives on 8 vertices, outside that band. Two larger censuses (all ordered pairs of non-cographs on at most 6 vertices, and a minimal-tight-pair sweep) were dispatched in the discovery run and never delivered; those channels are UNREAD, not negative, and nothing above leans on them.
> NOT RE-RUN: the family claim for general n > 12 is the hand argument printed in the paper. Machine coverage here is n = 5..12 only, for both bounds.
> NOT RE-RUN: no prior-art search. OpenAlex citation listing, Semantic Scholar keyword search and MathSciNet were all UNREAD in the accompanying prior-art pass; see the paper's bibliography for what was read.
