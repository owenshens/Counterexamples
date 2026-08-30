# Counting Collinear Triples, Four Distinct Triangles Force at Most Seven Points, and the Regular Heptagon

`four-distinct-triangles-force-at-most-seven-points-and-the-regular-heptagon`

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

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    8235ca9efc5ae55e3f028461e873eeedb9e3188573cfc5fa60084b2ffa431e92

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: Schoenberg's criterion itself (a point set with squared-distance matrix M embeds in R^d iff -M is positive semidefinite of rank <= d on the sum-zero hyperplane) is quoted from the literature, not proved here.
> NOT RE-RUN: part (c) of Theorem 1.1, the step from seven points to eight.  It is a hand argument, Proposition 5.1 of the paper, and contains no computation.
> NOT RE-RUN: any claim about the NONCOLLINEAR convention beyond the single hexagon-plus-centre configuration checked in B.  In particular this program does NOT decide the maximality half under that convention, and the paper does not claim it.
> NOT RE-RUN: minimality or uniqueness of the certificates used above.  The one-signed test and the forced-zero test are SUFFICIENT conditions for non-realizability; nothing here claims they are the shortest such certificates, or that the 30 orbits could not be eliminated by a different route.
> NOT RE-RUN: the published classifications this result is adjacent to (Shinohara 2004 on planar three-distance sets, Erdos-Fishburn 1996 on g(k)).  The two-distance lane is re-derived above; the three-distance classification is NOT used anywhere in this proof and is NOT re-derived.
