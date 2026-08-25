# A Stabilization Counterexample to Conjecture 4.2 of Boix and Lima–Pereira

`tor-stabilization-refutes-boix-lima-pereira-conjecture-4-2`

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
passes. The recorded run reports **26 checks, all passing**:

    VERDICT: ALL 26 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    a950df1bbd027c8746abdf9aaabcf0c3414199a838c3410e2fc374858a984278

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the two cited results that force k >= 3 (the k = 1 and
> k = 2 cases of the identity) are literature statements; they are
> supported here on a finite census, not proved.  The stabilization
> statement is a homological argument valid for every gap; it is
> instantiated here only at gaps 1 and 2 for this example, and no
> exhaustive census over complete intersections was attempted.
> NOT RE-RUN, SCOPE OF THE CENSUS: the census runs in three
> variables, so every instance has n = 3.  Its k = 3 instances have
> k = n, which is the unrestricted assertion already disproved in
> the literature, NOT the k <= n-1 range of the conjecture under
> test; the k = 3 failures counted above are therefore not that many
> independent counterexamples to that conjecture, and they
> enter its range only through the stabilization proposition, which
> is applied here only to the paper's own instance (gaps 1 and 2).
> The in-range refutation rests on identity_fails_refutation with
> stabilization_gap_1 and stabilization_gap_2 alone.  Only the k = 1
> and k = 2 halves of the census are in range, and they are
> confirmations of the cited results, not tests of the conjecture.
