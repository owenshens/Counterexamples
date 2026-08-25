# Nine additional counterexamples to Conjecture 5 of Sondow, Nicholson, and Noe

`ramanujan-exceptions-refute-sondow-nicholson-noe-conjecture-five`

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
passes. The recorded run reports **28 checks, all passing**:

    VERDICT: ALL 28 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    5fc1042f36e2e49f542a864e2244bdbd2c50714f998b316479189061b42daf25

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: (i) pairs with mn > 21317 are not tested one by one; that range is covered by the two analytic ingredients above, whose quoted numerical inequalities are verified but whose reduction to W_m(n) and whose estimate rho(k) < 2k(1+gamma(k)) for k >= 5225 are quoted from the cited literature and are not reproved here; nor are the closed forms of Q', Q'', H, u_1', u_2' or the identity Q(t) = t^3 q(e^t), which this program evaluates in certified rational interval arithmetic but does not re-derive. (ii) Three bounds the argument needs on unbounded ranges -- Q'' > 2.46 for t >= t_0, gamma(x) < 101/1244 for x >= 5225, and the decrease of u_1 and u_2 for t >= t_0 -- are established here for those full ranges, by monotone-envelope reductions to finitely many rational inequalities that are all checked; the grids that remain (gamma's strict decrease at 25 geometric points up to about 1.4e18, Q's increase at 10 points) corroborate only and carry no load. (iii) Nothing here verifies the paper's attribution that the pair (38,9) was asserted in the earlier literature to be the unique exception: that (38,9) is an exception and that nine further pairs are exceptions is verified by exact integer arithmetic, but the content of the cited claim is an input taken from the paper and is not checkable from the files shipped here.
