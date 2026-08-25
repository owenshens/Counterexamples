# A Five-Vertex Counterexample to Tang's Bunkbed Self-Avoiding-Walk Question

`house-graph-refutes-tang-bunkbed-self-avoiding-walk-question`

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
passes. The recorded run reports **20 checks, all passing**:

    VERDICT: ALL 20 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    93938a4ecfd5face28807004dfa7f5c70d623fab3902f2a637aaf8deddb742b8

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> note: not re-run here: no exhaustive census beyond order 6 was attempted, so the search for further counterexamples stops at six vertices; the paper makes no claim past order five, whose minimality is settled by the complete order <= 4 census above.
