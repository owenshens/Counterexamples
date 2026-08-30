# An Order-49 Partition Refutes the Heavy-Chopped-Rectangle Interval Conjecture of Gottlieb, Krnc and Muršič

`an-order-49-partition-refutes-the-heavy-chopped-rectangle-interval-conjecture`

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
passes. The recorded run reports **87 checks, all passing**:

    VERDICT: ALL 87 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    73088e9c7083b0932e6269e611e8a9e9eae8f3d1dbbb91196ae8a3f046c41449

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> Scope
> NOT RE-RUN: (1) HEAVINESS OF WHOLE CELLS OUTSIDE THE ORDER-49 WINDOW. Step 12
>   sweeps every admissible cell whose least member has order <= 49, but only over
>   that cell's members OF ORDER <= 49. Cells with b >= 8, and the parts of cells
>   with b <= 7 lying above order 49, are NOT claimed clean -- they are only shown
>   unable to contain a counterexample of order <= 49. Only cells (10,4), (11,4)
>   and (8,7) are swept in full (Step 11), so the inventory of 43 is a lower bound
>   on the number of non-heavy interval members, not a total.
> NOT RE-RUN: (2) THE INFINITE FAMILY. The last check of Step 11 evaluates
>   [(a+1)^3, a^2] at a = 11, 19, 27 only. That every a == 3 (mod 8) gives deficit
>   exactly 7 is an OBSERVATION suggested by three data points and is NOT proved
>   here; the paper states it as a conjecture of ours and nothing else in this
>   program depends on it.
> NOT RE-RUN: (3) THE SOURCE TEXT AND THE PRIOR-ART SEARCH. The two appendix tables
>   and the closed forms of thm:rect, prop:bound, prop:hook, prop:padded_stair,
>   prop:2rows and thm:partition_interval are TRANSCRIBED from arXiv1.tex of
>   arXiv:2506.04991v1 into the constants at the head of this file; this program
>   cannot fetch that file and does not verify the transcription. Nor does it
>   check any literature claim: that no published theorem covers a 5-row order-49
>   partition is an editorial finding, not a computation.
> NOT RE-RUN: (4) THE FINITE RANGES ARE FINITE. thm:rect is checked for r,c <= 9
>   (and its closed form used, not reproved, elsewhere), prop:hook for r,c <= 12,
>   prop:padded_stair for c,r <= 10, prop:2rows for c1 <= 14, obs:lucas as a
>   binomial identity for a,b <= 64, thm:partition_interval as an iff to order 26,
>   and the two-engine agreement over the reachable set of lambda plus every
>   partition of order <= 9. Nothing outside those ranges is asserted.
