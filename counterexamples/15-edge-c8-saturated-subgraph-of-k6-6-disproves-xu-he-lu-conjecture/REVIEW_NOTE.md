# Counterexamples to the Xu–He–Lu Conjecture on Partite Cycle Saturation

`15-edge-c8-saturated-subgraph-of-k6-6-disproves-xu-he-lu-conjecture`

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
passes. The recorded run reports **135 checks, all passing**:

    VERDICT: ALL 135 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    bbafe01c2d7ff733c105dc5bb26667e9996e15cb0ba5d3efe95b52cd1e7b914f

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE NOTE -- exactly what is and is not machine-decided here.
> ESTABLISHED IN FULL: Conjecture 6.2 asserts an EQUALITY, so a single
> saturated subgraph with fewer edges refutes it.  The l = 4 instance
> alone does that, and it is decided here exhaustively: 15 < 17 with
> C_8-freeness and all 21 non-edges witnessed twice over.  l = 5,6,7
> are three further independent refutations, not supporting evidence.
> NOT DECIDED HERE (1): the 'for every l >= 4' half of Theorem 1, i.e.
> the cases l >= 8.  What is checked at each tested l is every general
> ingredient of the proof separately -- Lemma 2 and the tightness of
> its two index ranges, the sigma symmetry, the exhaustiveness of the
> four-case partition of the non-edges, the summand-by-summand edge
> count, and the structural cycle prediction (R-cycles cap at 2r,
> Q-cycles start at 2r+4) -- so the general argument is corroborated
> at each instance, but induction on l is not a computation.
> NOT DECIDED HERE (2), and not computable: whether the printed
> Discrete Math. 349 (2026) 114802 text agrees with arXiv:2410.11194v2
> in the quoted wording and in the labels 'Theorem 1.1'/'Conjecture
> 6.2'.  Control (d) transcribes Definition 2.2 from the v2 source, so
> it inherits the same v2 scope the paper itself discloses.  What (d)
> DOES corroborate is that the formula n1+n2+l^2-3l+1 attributed to
> them is the edge count of their own construction, at three instances.
