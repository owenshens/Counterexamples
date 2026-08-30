# A Double-Suspension Counterexample to the Composition Converse for Strong CW-Regular Subdivisions

`a-double-suspension-counterexample-to-the-composition-converse-for-cw-subdivisions`

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
passes. The recorded run reports **52 checks, all passing**:

    VERDICT: ALL 52 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    ff1a0718dc0c37dc3251476944f579d9860b6086fd261e3f22602124e4314d76

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the double suspension theorem, |P| * S^1 = |P| * S^0 * S^0 homeomorphic to S^5 (J. W. Cannon, Ann. of Math. 110 (1979) 83-112). This is the ONE input taken from the literature and the only support for the sphere clause of hypothesis (iii); nothing here reproves it, and no homeomorphism type is computed anywhere in this program.
> NOT RE-RUN: ex:joinballs (that Delta^n * bd Delta^{n'} is a ball with boundary bd Delta^n * bd Delta^{n'}) is quoted from the target paper; the program checks its combinatorial consequences -- dimensions, boundary face sets, pseudomanifold conditions and Euler characteristics at all 393 y -- not the ball structure itself.
> NOT RE-RUN: the classical identification of pi_1(|P|) with the binary icosahedral group of order 120. What is proved here is strictly weaker and strictly sufficient: pi_1(|P|) surjects onto A_5, hence is not trivial.
> NOT RE-RUN: the case rank(sigma) = 1 of the Question. No computation in this program bears on it; it is OPEN, not empty.
> NOT RE-RUN: any minimality. No census of smaller complexes, fewer vertices or other homology spheres was performed, and none is claimed; Bjoerner and Lutz already give non-PL 5-spheres on 18 vertices against the 19 of K_X.
> NOT RE-RUN: the integral homology of K_X, and the exact-iff statement on the Z = B_0 slice. K_X is checked only to be a closed 5-pseudomanifold with chi = 0.
> NOT RE-RUN: the literature. This program performs no search and no citation check.
