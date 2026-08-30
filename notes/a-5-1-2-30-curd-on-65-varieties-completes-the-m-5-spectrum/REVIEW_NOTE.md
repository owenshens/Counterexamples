# A 5^1 2^30-CURD on 65 Varieties, Completing the 5^1 2^b-CURD Spectrum

`a-5-1-2-30-curd-on-65-varieties-completes-the-m-5-spectrum`

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
passes. The recorded run reports **53 checks, all passing**:

    VERDICT: ALL 53 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    58aa5e7e893828542d79e0002c4bf6a628592a5293979f9217c1c4a661dc5986

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the SEARCH. The object above was found by Algorithm X on a 360-element exact cover
>   (240 (class, variety) elements plus 120 leave orbits, 5,280 candidate rows), at seed 1 in 130
>   nodes. This program does not repeat that search and could not reproduce its choices; nothing here
>   depends on it, because a positive certificate of this kind is self-checking.
> NOT RE-VERIFIED: the 'only if' half of the source's Theorem 2.9, that a 5^1 2^((n-5)/2)-CURD forces
>   n in {5, 25, 65}. It is quoted, not re-proved. What is re-derived here is only its arithmetic at
>   n = 65 and the two bounds of the source's Theorem 2.5 at t = 4.
> NOT REPRODUCED: the source's own objects at n = 5 (its Example 1.2) and n = 25 (its Theorem 4.7 /
>   Corollary 4.11 at q = 5). Only their parameter arithmetic is re-derived below; no object at either
>   order is built, parsed or checked here, so the completed spectrum {5, 25, 65} rests on the source
>   for two of its three members.
> NOT COVERED: uniqueness, enumeration, or the number of solutions at n = 65. This is a construction,
>   not a classification; the exhausted fraction of the unrestricted space is 0, and no non-existence
>   claim is made or checked anywhere in this file.
> NOT COVERED: the source's OTHER open cell, a 4^1 2^12-CURD on n = 28 (its line 361). This program
>   re-derives only that its parameters differ from ours (m = 4, w = 21, t = 3) and asserts nothing
>   about its existence in either direction.
> NOT TRANSCRIBED: Inequality (6) of Danziger-Stevens (2004), their general two-block-size extremal
>   bound. Its text extraction is layout-mangled -- the fraction bar and the scope of the leading
>   lambda are unrecoverable -- so it was not transcribed with certainty and is not evaluated here.
>   Only the published bound p_2 > 1 => v <= p_2^2 is checked.
> NOT SEARCHED: prior art beyond the channels listed in the 'Bibliographic channels were incomplete'
>   paragraph of the accompanying note
>   which full texts were unreachable. OpenAlex returned HTTP 429, MathSciNet and two Wiley full texts
>   were inaccessible, and the Handbook CURD section is print only. Nothing in this program bears on
>   those gaps.
