# A Counterexample to the Hibi–Seyed Fakhari Conjecture on Join–Meet Ideals

`a-counterexample-to-the-hibi-seyed-fakhari-conjecture-on-join-meet-ideals`

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
passes. The recorded run reports **64 checks, all passing**:

    VERDICT: ALL 64 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    e7141ba8053f3179c56fd7a3707f2b6e357ab5a3139691883b5cd76bc1117c81

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE NOT RE-RUN: this program does NOT re-derive the paper's full claim. The published theorem asserts beta_{1,j}=0 for EVERY j>=4, and only finitely many degrees can be recomputed. beta_{1,j}=0 for j>=7 is NOT recomputed here; the paper obtains it from the cited [HSF, Lemma 2.8] regularity bound reg(I_L)=4, a literature input, and no finite computation can cover infinitely many degrees. Degrees j=4,5,6 were recomputed in full; the paper itself computes only j=4,5. Also NOT re-run: full target rank of mu_j over a general field of characteristic p>=3 -- the ranks printed above are certified over F_2, over Q and hence over every field of characteristic zero, and over F_2147483647, and over no other characteristic.
> NOTE re-run with --deep to add degree 7 (mu_7 is 429144 x 121549 over F_2; needs about 1.3 GB of RAM and a few more seconds).
