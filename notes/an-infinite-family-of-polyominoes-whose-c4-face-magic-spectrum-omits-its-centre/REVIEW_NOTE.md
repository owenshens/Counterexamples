# An Infinite Family of Polyominoes Whose\\ C_4-Face-Magic Spectrum Omits Its Centre

`an-infinite-family-of-polyominoes-whose-c4-face-magic-spectrum-omits-its-centre`

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
passes. The recorded run reports **99 checks, all passing**:

    VERDICT: ALL 99 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    11012681b63de78e60d4b901fabcde4a724dfa71488b114a92cd3ae83b2f1395

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: this program checks the paper's claims and nothing else.
> NOT RE-RUN: the SPECTRA of the census of Section 6 -- which found 37
>   non-interval spectra among the 163 shapes, every gap exactly at the
>   centre -- are CORROBORATION and are not recomputed here. Only the
>   SHAPE population of that census is re-derived above (163 shapes, and
>   exactly one of them with a bounded C_4 face that is not a cell). The
>   census's k = 7 stratum rests on a single implementation that was never
>   independently re-exhausted, and its 37 is the count under the CELL
>   reading. No claim of the paper's theorems depends on any of it.
> NOT RE-RUN: attainment for L_m with m >= 8 is NOT verified. Lemma 3
>   excludes the centre for EVERY even m and Lemma 4 confines Spec(L_m) to
>   six values, both with no search; that Spec(L_m) equals all six values
>   is established here only for m = 2 and m = 6, by exhibited labelings.
> NOT RE-RUN: minimality. Nothing here asserts that L_2 is the smallest
>   polyomino with a non-interval spectrum, nor that the L_m are the only
>   such shapes.
> NOT RE-RUN: for the holed shape of Remark 9 this program checks only
>   that the hole is a bounded C_4 face and that one labeling is magic
>   over cells AND over the hole. Whether the two readings have DIFFERENT
>   attainable sets on that shape was decided by the single unreplicated
>   census above and is deliberately not claimed by the paper.
> NOT RE-RUN: the odd-m half of the parity dichotomy is not touched. It
>   would need Shiu-Low-Liu 2024, Theorem 2.3, whose 2025 corrigendum we
>   could not obtain; the paper therefore states the family only for m
>   even, which is forced by 4 | n and not chosen.
