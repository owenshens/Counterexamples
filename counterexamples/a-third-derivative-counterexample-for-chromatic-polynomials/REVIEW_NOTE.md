# A Third-Derivative Counterexample for Chromatic Polynomials

`a-third-derivative-counterexample-for-chromatic-polynomials`

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
passes. The recorded run reports **80 checks, all passing**:

    VERDICT: ALL 80 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    b8cc8534d4fb3c646e9376ecc25709a87f5efcdf5fb5eefa9e96d47078b91040

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE of Corollary 3: the bare criterion 'L^(3)(-1/2) >= L_{K_n}^(3)(-1/2)' is ALSO met by much smaller graphs -- the least order among noncomplete generalized theta graphs is 4, attained by Theta(1, 3), Theta(2, 2) (Theta(2,2) and Theta(1,3) are both C_4). So the upper inequality of Conjecture 3 is not specific to H; what is specific to H is that its L^(3) is POSITIVE, which is what refutes Conjecture 2.
