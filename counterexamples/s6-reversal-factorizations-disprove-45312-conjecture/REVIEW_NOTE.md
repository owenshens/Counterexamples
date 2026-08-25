# Counterexamples to Conjecture 5.4 of Parisi–Spahiu–Skandera–Wang

`s6-reversal-factorizations-disprove-45312-conjecture`

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
passes. The recorded run reports **35 checks, all passing**:

    VERDICT: ALL 35 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    09cb79d93fa2900c3dd28c08a4bac83d40eb2c10149812e1e0f955a8a3b633f0

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN (each item below is asserted by the paper and is NOT covered by any check above):
> 1. The 3412-gap values quoted for N1 (2 for twelve elements, 3 for 564312) are not recomputed: the gap statistic is defined in the cited reference and not in this paper. The gap>1 hypothesis of the cited Theorem 5.2 is therefore NOT verified here; only 4231-avoidance and 3412-containment are.
> 2. The 'only if' half of the classification -- that the nineteen remaining elements have NO reversal factorization -- is not re-derived by an exhaustive search over all interval sequences. This program verifies only the hypotheses on which the two cited obstruction theorems act, and neither obstruction theorem itself is reproved here.
> 3. The corollary is stated for every n >= 6, but only n = 7 is computed (pattern containment is checked for n = 7, 8, 9). The step from n = 7 to general n rests on the cited parabolic compatibility of the Kazhdan-Lusztig basis, not on an induction carried out here.
> 4. The degree bound deg P_{v,w} <= (l(w)-l(v)-1)/2 is enforced by the construction of the basis rather than tested against it; the independent evidence that the computed basis is the right one is positivity, P(0)=1, the two known values, the smooth count, the Bruhat supports, and the inverse/w0 symmetries.
