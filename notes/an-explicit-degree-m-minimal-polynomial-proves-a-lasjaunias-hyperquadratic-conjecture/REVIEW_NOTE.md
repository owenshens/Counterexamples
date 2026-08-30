# An Explicit Degree-m Minimal Polynomial:\\ the Factorisation Half of Conjecture 11 of Lasjaunias

`an-explicit-degree-m-minimal-polynomial-proves-a-lasjaunias-hyperquadratic-conjecture`

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
passes. The recorded run reports **81 checks, all passing**:

    VERDICT: ALL 81 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    bc2b807c53a2e058ff9df68ba50e69b66b0b747e400629df1ab04a6096fe0eca

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the continued-fraction CONSTRUCTION of alpha in Proposition 5.1 of the source.
> NOT RE-RUN:   PART E does not follow that construction; it solves H(V)(alpha) = 0 in F_p((1/T))
> NOT RE-RUN:   by Newton from alpha_0 = w_{k+1}T and then verifies the two properties the paper
> NOT RE-RUN:   actually uses -- H(V)(alpha) = 0 and first partial quotient w_{k+1}T, hence a pole
> NOT RE-RUN:   at T = infinity.  Uniqueness of the root with that polynomial part is not proved.
> NOT RE-RUN: exactness of alpha. PART E certifies its identities to order u^20 only; the step
> NOT RE-RUN:   from that to an exact root is Hensel's lemma in the complete field F_p((1/T)),
> NOT RE-RUN:   whose hypothesis v(H) > 2 v(H') is re-derived (E02) and whose conclusion is cited.
> NOT RE-RUN: PART E for p > 60: the alpha census stops there, PARTS A-D run to p <= 300.
> NOT RE-RUN: irreducibility of P over F_p(T). Only the hypothesis of Capelli's criterion is
> NOT RE-RUN:   re-derived here (A20, B26: v_{T-1}(a) = +-1, coprime to m); the criterion is cited.
> NOT RE-RUN: nu(alpha) = m, the third clause of the conjecture. That is Corollary 4.2 of the
> NOT RE-RUN:   source (published Corollary 5, from Ayadi-Lasjaunias 2016) and is not our result.
> NOT RE-RUN: the full factorisation multisets of H(V) over F_p(T) -- no factorisation algorithm
> NOT RE-RUN:   is implemented; only exact division by a given monic P.
> NOT RE-RUN: primes p > 300. The theorem is proved for the whole infinite family; the census
> NOT RE-RUN:   above exhausts 0% of it and is corroboration, not evidence for the theorem.
