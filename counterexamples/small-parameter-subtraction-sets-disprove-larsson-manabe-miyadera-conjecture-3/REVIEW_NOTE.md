# Small-Parameter Counterexamples to Conjecture 3 of Larsson, Manabe, and Miyadera

`small-parameter-subtraction-sets-disprove-larsson-manabe-miyadera-conjecture-3`

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
passes. The recorded run reports **24 checks, all passing**:

    VERDICT: ALL 24 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    f7f9867eb7576863ceea2ab5675437b91a65d8b764785f5ac180fef395bf35ce

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: the conjecture quantifies over all a and n and all x above the threshold.
> Refutation is computed exactly for a = 1..20 (both pass conventions); the reduction
> to a = 1 for all larger a rests on the scaling lemma, which is verified numerically
> for a <= 8 and k <= 60 but not proved here.  The census over n is finite (n <= 12) and
> is also bounded in x, where the paper's census sentence is not: that sentence says that
> for every other n the two identities hold above the stated threshold, i.e. at every
> x >= threshold, whereas each non-exceptional n is certified here only for threshold <= x
> <= 400, so a break at larger x would not be seen.  The three exceptional n the theorem
> exhibits are settled outright at their witnesses and do not rest on the census.  The
> eventual-periodicity conjuncts are not tested, and the paper asserts nothing about them.
> No external catalogue or table is read.
