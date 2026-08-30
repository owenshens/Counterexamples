# Lin's Hot–Cold Maximum-Determinant Correspondence at Order Ten

`lins-hot-cold-maximum-determinant-correspondence-proved-at-order-ten`

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
passes. The recorded run reports **43 checks, all passing**:

    VERDICT: ALL 43 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    92a2e214a19ae8a9424e337e08c5266b98b2427f40b992c51cdf2e92883aeb2a

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: only the cell n = 10 is decided here. Orders 18, 22 and 30 -- the remaining cases of Lin's question below 30 -- are not touched, and the bordering-plus-covering method does not reach them: n = 18 would need the isomorphism classes of 17-vertex tournaments.
> NOT RE-RUN: the four censuses reported in our own working notes (two C programs, one numpy batch pass, one min-code canonical form) are not re-executed; this program is a fourth, independently written, pure-Python census, and it agrees with them.
> NOT RE-RUN: nothing here checks a bibliographic claim. That f_k(10) = 33489 was already proved by Klanderman, Montee, Piotrowski, Rice and Shader, and that both values were already reported by Alvarez, Armario, Frau and Gudiel, are statements about the literature; the only literature datum this program touches is the published order-10 Pfaffian spectrum, which it reproduces.
> NOT RE-RUN: the identification of this witness with the order-10 cocyclic matrix of Alvarez et al. is deduced from uniqueness, not computed -- that paper prints no order-10 matrix to compare against.
