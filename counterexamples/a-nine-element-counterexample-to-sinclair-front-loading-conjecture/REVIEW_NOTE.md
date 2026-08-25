# A Nine-Element Counterexample to Sinclair's Front-Loading Conjecture

`a-nine-element-counterexample-to-sinclair-front-loading-conjecture`

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
passes. The recorded run reports **72 checks, all passing**:

    VERDICT: ALL 72 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    25a9c62def7860a0b6eb0ad8bef5fa12bbbd9317b90f342749f5fe331524e34e

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> CERTIFICATE SCOPE: the load-bearing dimensions 31 (k=2) and 32 (k=3)
> CERTIFICATE SCOPE: do not rely on the ITERATIVE chain -- their
> CERTIFICATE SCOPE: generators are all words of length k applied to
> CERTIFICATE SCOPE: g_E -- but they DO use the same elimination as
> CERTIFICATE SCOPE: span_basis to SEARCH for their witnesses:
> CERTIFICATE SCOPE: independent_indices selects the rows of the minor,
> CERTIFICATE SCOPE: rref selects its columns, kernel_basis proposes
> CERTIFICATE SCOPE: the annihilators.  Only the WITNESSES are verified
> CERTIFICATE SCOPE: independently: a nonzero Bareiss determinant for
> CERTIFICATE SCOPE: the lower bound, and every annihilator dotted
> CERTIFICATE SCOPE: against every generator plus a Bareiss minor for
> CERTIFICATE SCOPE: their independence for the upper bound.  A bug in
> CERTIFICATE SCOPE: the elimination can lose a witness, not invent one.
