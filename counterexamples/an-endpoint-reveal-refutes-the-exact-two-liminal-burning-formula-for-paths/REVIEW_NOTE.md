# An Endpoint Reveal Refutes the Exact Two-Liminal Burning Formula for Paths

`an-endpoint-reveal-refutes-the-exact-two-liminal-burning-formula-for-paths`

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
passes. The recorded run reports **53 checks, all passing**:

    VERDICT: ALL 53 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    4723944b34def8a71151614d7a43d4ebf608784666a21280d614b8e1a8d2d7d8

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the provenance of the target statement. The label `thm:2-liminal burning path graph`, the file name section_3/section_3.tex, its line 4, the numbering setup in theorem_setup.tex, the 24,463-byte e-print size, and the line numbers quoted from arXiv:2505.10727 were established by fetching and reading those sources; this program re-fetches nothing and checks mathematics only.
> NOT RE-RUN: the two published closed forms are used as EXTERNAL CHECKSUMS, not proved here. b(P_n) = ceil(sqrt(n)) is Bonato-Janssen-Roshanbin's. CL(P_n) = floor(n/2)+1 is quoted from arXiv:2606.12330; its primary source (Bonato-Marbach-Milne-Mishura, WAW 2024) was not accessible to us, so the b_1 column is our own computation agreeing with a secondary quotation, NOT an independent external check.
> NOT RE-RUN: n >= 15. The census is complete for 1 <= n <= 14 and says nothing beyond it. No closed form for b_2(P_n) is established, and none is claimed.
> NOT RE-RUN: the reveal-cardinality reading. The engine has the saboteur reveal exactly min(k, |pool|) vertices. Groups B2-B6 prove b_2(P_4) = 3 under both the exactly-k and the at-most-k readings, but the census values for n >= 5 are computed under exactly-k only.
> NOT RE-RUN: readings outside the four tabulated. Four combinations of the two ambiguous clauses were tested; a fifth reading nobody has proposed is not covered.
> NOT RE-RUN: the prior-art channels (arXiv, Semantic Scholar, OpenAlex, zbMATH, OpenCitations) that bound the novelty of this correction. OpenAlex citer enumeration returned HTTP 429 and was never read, and MathSciNet was never consulted.
