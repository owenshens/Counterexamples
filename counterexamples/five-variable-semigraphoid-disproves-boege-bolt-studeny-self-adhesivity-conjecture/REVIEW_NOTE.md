# A Counterexample to a Self-Adhesivity Conjecture of Boege, Bolt, and Studený

`five-variable-semigraphoid-disproves-boege-bolt-studeny-self-adhesivity-conjecture`

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
passes. The recorded run reports **20 checks, all passing**:

    VERDICT: ALL 20 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    aeebbb009953639a47a42cc822a3d837f05d91f6655fd00fa8d3cac5d3f202a8

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> ASSUMED (a) one canonical copy suffices [BBS Def 4.3 Rmk 1, Lem 4.5]
> ASSUMED (b) sg^diamond(M|L) = sg(A_L)|N  (definitional in [BBS])
> ASSUMED (c) only |L| in {2,3} need computation [BBS Cor 4.12, Lem 4.15-4.16]  -- NOT relied upon: this verifier computes all 32 overlaps
> ASSUMED (d) marginal equality plus the adhesion condition force delta_s = 0 for every s in F and every supermodular g inducing the adhesion
> ASSUMED (e) 'f is not a coatom' / not among the five-variable coatoms of [BBS Sec 7.6] -- the rank-19/face-dim-7 arithmetic IS checked here, but the step from face dimension 7 to 'not a coatom' is not
> ASSUMED (f) minimality of five variables [BBS Cor 4.17, Sec 7.4] -- nothing about |N| = 4 is computed in this paper or in this verifier
