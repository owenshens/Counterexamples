# A Proof of the Bhagat–Kulkarni–Larsson–Murali\\ Smallest-Discrepancy-Heap Formula

`a-proof-of-the-bhagat-kulkarni-larsson-murali-smallest-discrepancy-heap-formula`

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
passes. The recorded run reports **32 checks, all passing**:

    VERDICT: ALL 32 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    893c15fea8e299b10e65fcae5d51b7d95a4998467948dbf8f192884a050b3412

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN -- what this program does NOT cover:
>   * Heaps h > H.  Above H the outcome genuinely is convention-dependent and neither
>     the paper nor this program claims anything there.
>   * The excluded ratios s_1/s_2 = (k+1)/k (i.e. d | a): the target's COMPANION
>     conjecture.  Only one control touches them, and only to show that the paper's
>     Lemma 4(F2) genuinely fails there; no claim of any kind is checked at d | a.
>   * Subtraction sets with |S| >= 3.
>   * The induction itself.  Theorem 2 is proved by hand for every n >= 1; this program
>     audits it on the finite box s_2 <= 100 and cannot certify the general case.  It is
>     not a proof assistant and the proof has not been formalised in one.
>   * Pairs with s_2 > 100.  The box is finite and stated; a wider census (s_2 <= 200,
>     19002 pairs) was run elsewhere for the CONJECTURE only and is not re-run here.
>   * The bibliographic exposure named in section 5 of the paper: the in-preparation
>     manuscript Kulkarni-Larsson cannot be fetched, and no program can settle that.
