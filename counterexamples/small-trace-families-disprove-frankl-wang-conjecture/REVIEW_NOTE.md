# Small Counterexamples to a Trace Conjecture of Frankl and Wang

`small-trace-families-disprove-frankl-wang-conjecture`

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
passes. The recorded run reports **56 checks, all passing**:

    VERDICT: ALL 56 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    332a38f5e2cf986cfbebbec1d479d36bbc40731662020daf69e4b3afae4ae2b3

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE not re-run: the squashing lemma (down-set reduction) is quoted from the cited literature, not proved; it is exercised on 300 random 40-member families at n=6 and confirmed exhaustively at n=4. Two results above rest on it. First, the upper bound m(6,5,25)<=40, which uses that one quoted lemma plus the exhaustive down-set census run above. Second, the auxiliary equalities m(5,3,7)=13 and m(5,4,13)=19, whose searches range over the 7581 down-sets of 2^[5] only, so without the lemma they are only the lower bounds m(5,3,7)>=13 and m(5,4,13)>=19; of the auxiliary values, just m(3,3,7)=7 and m(4,3,7)=10 come from unrestricted brute force over every family. The three lower bounds 40, 58, 80, which are what refute the formula, use no external input. No census is attempted at n=7.
