# A Length-14 Negative Answer to Alderson's q=4 Optimality Question

`a-length-14-negative-answer-to-alderson-q-4-optimality-question`

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
passes. The recorded run reports **71 checks, all passing**:

    VERDICT: ALL 71 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    381a65cf99c1c9cabd3ab031fefa06645b4ccd0637a035f0e64176f68532d1d0

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE NOT RE-RUN: five things above are inputs taken from the literature and are recomputed by no check in this transcript. (1) Alderson's Problem 9.2 itself, together with his Theorem 6.6 and its parameters (112,2,104)_{16/4}: this program reads no external paper, so the statement being answered and the comparison length 112 are transcriptions, and a misreading of either would be caught by nothing above. (2) The definitions of 'extendable' and of 'additively maximal': the 7 properties counted in problem_9_2_hypotheses_all_verified_for_this_code are this program's rendering of those definitions, not a derivation of them. (3) The Remark's projective-equivalence claim. The 112 lines external to B counted above are external to THIS B; that Alderson's scattered F_2-linear set for t=2 is a Baer subgeometry projectively equivalent to it, and hence that these fourteen lines are fourteen of his own 112 external lines, is not checked here. (4) The spread classifications of van Dam and of Mellinger. The fourteen explicit lines are verified above to partition Pi together with B, which is a self-contained substitute, but no classification is confirmed and it is not tested that this spread lies in one of their three orbits. (5) Faithfulness, tested only in the paper's own sense that every coordinate projection is onto F_16; a stronger notion is not tested. Also NOT attempted: any lower bound, and any search for the global minimum length -- what is established above is 14 < 112, not that 14 is optimal, and the nonsquare, nonprime clause of Problem 9.2 is untouched. Apart from the five items above, the only paper values this program reads are the three determinant tuples and the weight enumerator W_C, and those are quoted solely as comparison targets printed beside independently recomputed values.
