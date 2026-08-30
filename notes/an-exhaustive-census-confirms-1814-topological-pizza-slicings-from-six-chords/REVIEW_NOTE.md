# An Exhaustive Census Confirms a(6) = 1814 for Topologically Distinct Pizza Slicings

`an-exhaustive-census-confirms-1814-topological-pizza-slicings-from-six-chords`

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

    bdca61d83544b22776d7eed1c9e90778a8cb967b86a14d22e927d345059dc932

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN -- what this program does NOT establish:
>   * NOT RE-RUN: the SUPERSET LEMMA itself.  That two straight chords meet at most once and that a straight chord meets each convex face in at most one segment are proved in Section 3 of the paper by hand; this program re-runs the enumeration those facts license, not the facts.
>   * NOT RE-RUN: any n >= 7.  A272906(7) and beyond are untouched, as are A241600 (known for n <= 7) and A090338 (known for n <= 9).
>   * NOT RE-RUN: the sister sequence A273280 (chords of a SQUARE; 1, 1, 2, 5, 19, 129, 1806, Giovanni Resta, May 2016), which carries the identical "unconfirmed" caveat and is not attempted here even though the method applies.
>   * NOT RE-RUN: the identification of the 43 classes of the X = 15 shard with the 43 chirotopes of Christ's database of simple arrangements of 6 lines.  Only the two COUNTS are compared, 43 = A090338(6).
>   * NOT RE-RUN: Jon Hart's May 2016 guided random trials, which produced the value 1814 in the first place.  The value is his; what is re-derived here is the exhaustive upper bound that makes it a theorem.
>   * NOT RE-RUN: any claim that Table 1 is minimal.  47 point sets suffice; nothing here says fewer do not.
> checks: 53 run, 53 passed, 0 failed, 151.8s
