# A Counterexample at n=6 to Theorem 1 and Conjecture 36(3) of Kundu and Velmurugan

`s6-covering-number-refutes-kundu-velmurugan-conjecture-36-3`

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
passes. The recorded run reports **46 checks, all passing**:

    VERDICT: ALL 46 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    1632b8cf0cc14fde08be5c8c965f69c0931c9c885b9a52fedbcab59987fc9fec

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE: the 6 assertions of the paper NOT covered by the checks above
> 1. The abstract's priority claim, that the value 5 is already
> implicit in the power formulas of Sun-Zhang-Zhu, is a
> statement about the literature and is not tested here.
> 2. The remark's claim that the LOWER inequality of Conjecture
> 36(3) remains valid at n = 6 is not tested: the paper does
> not reproduce the lower bound, so there is no expression to
> evaluate.  Nothing above depends on it.
> 3. The remark's disclaimer that the counterexample does not
> address even n >= 8 is a non-claim; the census above stops
> at n = 11 and asserts nothing beyond it.
> 4. That Kundu and Velmurugan's SageMath run in fact produced
> the bound for n <= 11 cannot be checked from outside their
> paper.  The census shows only that the bound itself is
> false at n = 6, hence that any verification affirming it
> over n <= 11 must be in error.
> 5. Miller's theorem is used here only at n = 6, where it is
> re-derived by direct computation; the general n > 4
> statement is not tested and no check relies on it.
> 6. The WORDING of the two statements acted on -- Conjecture
> 36(3) and Theorem 1 of Kundu and Velmurugan, including the
> latter's exclusion set -- is transcribed from their paper.
> This program reads no literature, so it instantiates those
> quotations at n = 6 and compares them with the derived
> value 5; it cannot confirm that the quotations are faithful.
> Everything else in the corollary is derived here.
> NOT RE-RUN: the paper's full text.  What is re-derived from scratch, with no paper value taken on trust, is Theorem 1 -- both covering numbers, from an S_6 character table built twice by independent routes (Murnaghan-Nakayama and Jacobi-Trudi) and checked against row and column orthogonality -- every arithmetic step its proof displays, and, for Corollary 2, the two right-hand sides at n = 6, the exclusion set, and the strict inequality 5 > 4.  The 6 assertions listed under SCOPE above lie outside that: (1) the abstract's priority claim to Sun-Zhang-Zhu, a statement about the literature; (2) the remark's claim that the LOWER inequality of Conjecture 36(3) still holds at n = 6, which the paper never reproduces; (3) anything beyond n = 11, the census range, so the remark's even n >= 8 disclaimer is untested; (4) whether Kundu and Velmurugan's SageMath verification was in fact carried out, which is unobservable from outside their paper; (5) Miller's general n > 4 theorem, re-derived here only at n = 6, the sole case any check uses; (6) the faithfulness of the quotations of Conjecture 36(3) and of Kundu and Velmurugan's Theorem 1, which are transcribed here, not read from the source.
