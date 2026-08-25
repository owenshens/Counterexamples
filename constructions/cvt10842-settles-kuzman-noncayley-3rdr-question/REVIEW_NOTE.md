# A Non-Cayley Vertex-Transitive\\ 3-Rainbow Domination Regular Graph

`cvt10842-settles-kuzman-noncayley-3rdr-question`

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
passes. The recorded run reports **26 checks, all passing**:

    VERDICT: ALL 26 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    986e591a9bd5909fbbce2731fc0abed2848218c9e43dbfec99f6c8392079e805

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: (i) the identification of X with the entry CVT[108,42] of the cubic vertex-transitive census -- asserted in the paper's abstract and in the first sentence of its proof -- and any minimality over that census, both of which would need the external census catalogue; (ii) the general lemma that a 3-rainbow dominating function of weight w yields a dominating set of size w in G box K_3, equivalently gamma_r3(G) = gamma(G box K_3), which is cited rather than machine-proved: its conclusion is checked above on every automorphic image and every one-colour extension of the exhibited function, and the equality is checked exhaustively on 5 small cubic graphs, but it is not proved for all cubic graphs.  This program verifies the exhibited graph itself.
