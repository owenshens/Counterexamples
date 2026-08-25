# The Saturation Number of K_3 2K_2

`the-saturation-number-of-k3-union-2k2`

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
passes. The recorded run reports **60 checks, all passing**:

    VERDICT: ALL 60 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    101f80647e43bad25182fab24646cf35de9d2ca2e45fc189f903135af4215cf3

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: this program is not the program that produced the paper's census; it is an independent and deliberately narrower re-verification written from the paper's text, and it reads no catalogue file.  The paper's census is exhaustive over McKay's catalogues of all 1,008,629 isolate-free graphs with 1..14 edges (orders up to 28) and over the 25,937,431 one-edge / one-leaf / one-K_2 extensions of the 740,226 fourteen-edge representatives.  This program re-runs that census only for orders <= 13 and <= 14 edges (self-generated), which suffices for sat(n) with n <= 13; it does NOT re-run orders 14..28 of the 14-edge layer nor any part of the 15-edge layer, so the values sat(n)=15 for n >= 14 and the uniqueness of K_6 u I_{n-6} for n >= 14 are NOT reverified here.  All three enumeration figures (1,008,629, 740,226 and 25,937,431) ARE recomputed here by Polya counting, but that fixes only the SIZE of those layers, not the F-saturation test over them.  Finally, the independent Polya cross-check of this program's own unpruned generator covers e=1..10 of the 14 levels; levels 11..14, including the top level of 90506 representatives on which the bound sat(13) > 14 rests, are supported by the generator's extension argument together with the certificate-pinning and small-order cross-validation checks, not by a second independent count.
