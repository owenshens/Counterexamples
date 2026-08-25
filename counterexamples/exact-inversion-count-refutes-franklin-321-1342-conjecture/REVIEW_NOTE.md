# A Counterexample to a Conjecture of Frankl\'in on (321,1342)-Avoiders

`exact-inversion-count-refutes-franklin-321-1342-conjecture`

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
passes. The recorded run reports **101 checks, all passing**:

    VERDICT: ALL 101 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    8a618a677b3fcf918542ebf777134a89bb5699c16b81c84f20cc61a2e3881137

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE / WHAT IS NOT COVERED BY ANY CHECK
> 1. Every claim is verified on the finite ranges printed above.
> The theorem's 'for every k>=0, n>=1' is not established by
> exhaustion; the finite evidence is the enumeration, and the
> algebraic chain (eq. (J) = the 4-fold sum -> telescoping ->
> coefficient extraction) is verified as series identities to
> the printed orders, not symbolically for all k.
> 2. Lemma P11 (indecomposable => inv >= len-1) is imported from
> the literature.  It is verified here on ALL permutations of
> length <= 8 and on every indecomposable class member found,
> plus the boundary length 12, but not proved.
> 3. The BGU structural cover (P8) is verified as a FACT for
> n <= 8, not derived from BGU's proof.
> 4. Asymptotic statements in the paper ('grows like n^3/6',
> '1+C(k,2) and 1+C(k+1,2) are both asymptotic to k^2/2') are
> not finitely checkable and are NOT checked here.
> 5. Bibliographic claims -- the quoted sentence, the page and
> proposition locators, and the attribution of the conjecture
> to numerical evidence -- are outside this program.
> 6. In particular, the calibration of STAGE 7 does NOT close the
> gap left by item 5.  It shows that THIS program's counter
> indexes by inversions with I_0 = {1}, in Av(132), Av(231),
> Av(12) and Av(321) at once, against four models built here
> from unrelated definitions.  That these same four
> enumerations, at these same indices, are the ones proved in
> the paper's reference [2] is an attribution of the same
> bibliographic kind as item 5 and is NOT checked.  Hence the
> paper's sentence 'in [2] as here, k is the number of
> inversions of the permutation and I_0={1}' -- on which the
> whole refutation depends, since a different convention in [2]
> would make the conjectured count correct -- rests on the
> quotation and locator of item 5, and on no check below.
> CHECK RESULTS
