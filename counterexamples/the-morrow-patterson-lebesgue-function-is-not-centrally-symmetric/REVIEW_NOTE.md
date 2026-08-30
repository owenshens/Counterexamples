# The Morrow–Patterson Lebesgue Function\ Not Centrally Symmetric at n=2

`the-morrow-patterson-lebesgue-function-is-not-centrally-symmetric`

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
passes. The recorded run reports **39 checks, all passing**:

    VERDICT: ALL 39 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    a5aec2e60725a47e1f7174e35c87fb9c59ea45429cb5d00a6d2ef80702bacbba

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: (1) the FAILURE of central symmetry at n >= 4. Only the integer parity mechanism of section 5.8 is re-derived for even n up to 40; the numerical minima of lambda_n over -MP_n reported for n = 4..30 need multiprecision transcendental arithmetic and are outside the standard library, so they are NOT recomputed here.
> NOT RE-RUN: (2) the uniform candidate family and the negativity of its one alternating sum, which the paper states is unproved and only computed; nothing here bears on it.
> NOT RE-RUN: (3) the 101x101 grid census at n = 30 and the diagnosis of the authors' plotting script; both are corroboration reported elsewhere and neither is a claim of this paper.
> NOT RE-RUN: (4) the corner values of section 5.9 are checked but are NOT claimed as a result of this paper -- they belong to the separate corner-attainment question and are carried only as an independent second witness to the asymmetry.
