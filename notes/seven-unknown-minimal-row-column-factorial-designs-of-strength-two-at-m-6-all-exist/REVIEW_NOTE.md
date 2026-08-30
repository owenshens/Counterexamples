# Seven Row-Column Factorial Designs of Strength Two at m=6 Listed as Unknown by Rahim and Cavenagh

`seven-unknown-minimal-row-column-factorial-designs-of-strength-two-at-m-6-all-exist`

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

    174f359ef2ec03c35080bee8b4c09262f0d9d862860faa4ad9914889916811c7

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the cells 12 <= k <= 20 at definition level. Their arrays have up to 2^20 columns, which is not referee-sized; what is re-checked here for every K in 8..23 is the hypotheses of Theorem 1, from which those cells follow. They are in any case NOT new -- they are Rahim and Cavenagh's own published theorem for 12 <= k <= 20.
> NOT RE-RUN: nothing here bears on the whole m = 6 case of the conjecture. This program checks the seven MINIMAL cases the source lists as unknown; the reduction of a general admissible (k,n) at m = 6 to those minimal cases is not verified here, and no such reduction is proved in the source.
> NOT RE-RUN: the search over CHOICES OF COLUMN SET. The generating column set is the first eight columns of M and is fixed in the source of this program; only the 21,252 admissible generator triples {h1,h2,h3} on THAT one column set -- the 759 subspaces D they span -- are enumerated. The exhaustive C(23,8) = 490,314 census over column sets was never completed (the original run sampled 20,000 subsets, 4.079%), and existence needs none of it. Nothing here claims this witness is unique, canonical or least.
> NOT RE-RUN: the source's own printed examples at the level of their DIGITS. Those arrays are the source's, they are not transcribed into this program, and no cell of them is counted here. What section 5b checks about the reading of the definition is their printed DIMENSIONS, transcribed with their labels, and the parameter arithmetic of the two readings. The passages of the source quoted or paraphrased in Section 1 of the paper were read by hand and are not machine checkable at all.
