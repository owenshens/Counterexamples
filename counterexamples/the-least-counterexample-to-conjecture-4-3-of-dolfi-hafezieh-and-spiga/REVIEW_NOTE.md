# The Least Counterexample to Conjecture 4.3 of Dolfi, Hafezieh, and Spiga

`the-least-counterexample-to-conjecture-4-3-of-dolfi-hafezieh-and-spiga`

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
passes. The recorded run reports **36 checks, all passing**:

    VERDICT: ALL 36 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    a23614b99c13aaf54da969bf59accb6221bdad7e728dbbcf9035f7e1d87a4515

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE: the census in E and F is the paper's COMPLETE minimality
> claim: all 169 integers below 170 are enumerated, the 23 with
> omega >= 3 are tested against the gcd hypothesis, and each of
> the 4 survivors has 2^n-1 factored completely.  Nothing in the
> paper's claim is left unverified or narrowed.
> NOT re-run: (i) the supplementary sweep G1 stops at n = 100,
> because factoring 2^n-1 for prime exponents near 170 (e.g.
> n = 167) is far outside the budget -- G1 is an extra, not part
> of the paper's claim, whose Lemma 2 is proved, not computed;
> (ii) the bibliographic attributions in Remark 3 ([LLRS]
> Example 3.8, OEIS A046800) are not checkable offline.
