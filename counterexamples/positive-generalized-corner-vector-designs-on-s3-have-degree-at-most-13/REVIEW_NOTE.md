# Positive Generalized Corner-Vector Designs on ^3 Have Degree at Most 13

`positive-generalized-corner-vector-designs-on-s3-have-degree-at-most-13`

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
passes. The recorded run reports **48 checks, all passing**:

    VERDICT: ALL 48 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    5e1de55cadafb6d86c4ea4d3a1230c92392a107ff9e42d4b7d880d1e1e450758

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: (i) the SHARPNESS half of the corollary -- the two positive-weight 13-designs of this type reported in the cited paper, whose nodes and weights are given there numerically to six significant figures -- is not verified here, as the paper itself states; (ii) claims about the cited literature (its definition of a weighted design, and its record that signed 15-designs exist on S^3) are read, not computed; (iii) two elementary sign facts are used as definitions rather than computed: z >= 0 and 1-z >= 0 for z in [0,1], which is what the Bernstein interval means, and a^2 > 0 implies 0 < a^2/(a^2+s) < 1; (iv) this folder ships no mutation log, no corrupted-run transcript and no second implementation -- the falsifiability evidence here is the controls printed beside the checks, which are re-run on every invocation.
