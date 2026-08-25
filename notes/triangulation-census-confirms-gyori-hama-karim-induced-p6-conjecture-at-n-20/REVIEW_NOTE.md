# The n=20 Case of the Induced Planar P_6-Turán Conjecture

`triangulation-census-confirms-gyori-hama-karim-induced-p6-conjecture-at-n-20`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |
| `verify_planarity_and_p6_freeness.output.txt` | a recorded run of `verify_planarity_and_p6_freeness.py` |
| `verify_planarity_and_p6_freeness.py` | An independent re-derivation of the paper's two exhibited graphs from their graph6 strings, certifying planarity by exhibiting a sphere face set and confirming induced-P6-freeness by two mutually independent exhaustive searches. |

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **32 checks, all passing**:

    VERDICT: ALL 32 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    68e116cfc183e62700f30b90a57ce1f94811ff0877c007a1616890e47375f8ed

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the plantri census behind the UPPER bound -- no plantri run was made here, none of the 78,435,562 minimum-degree-four isomorphism classes over 6 <= n <= 20 was generated or tested, neither published count 11,284,042 nor 64,719,885 was reproduced, and the recursion giving (|L_n|)_{n=4..20}, in particular L_19 = {T} and L_20 = empty, was not re-run; the 34 facial insertions of T ARE checked above, but they close the degree-three branch of the upper bound only given L_19 = {T}, which is not verified here.
