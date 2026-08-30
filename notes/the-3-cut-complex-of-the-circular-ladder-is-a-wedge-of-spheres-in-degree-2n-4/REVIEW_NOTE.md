# The 3-cut complex of the circular ladder is a wedge of spheres in degree 2n-4, for n 5

`the-3-cut-complex-of-the-circular-ladder-is-a-wedge-of-spheres-in-degree-2n-4`

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
passes. The recorded run reports **45 checks, all passing**:

    VERDICT: ALL 45 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    9d8f711929de37c4e7603f103c604a1660e99e920fab1a587e12c630ddb15178

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN, and therefore not evidence of anything here:
>   * n >= 25.  The exact homology stops at n = 24 and the combinatorial claims at n = 40.
>     The statement for ALL n >= 5 rests on the PROOF in the paper, every quantitative step of
>     which this program confirms at 20 separate values of n.  No finite run closes it, and
>     this one is not an induction.
>   * The DIRECT (absolute, non-relative) homology of Delta_3 beyond n = 7.  At n = 7 the
>     complex already has 16228 nonempty faces and at n = 9 it has 261909; n = 4, 5, 6, 7 are
>     computed BOTH ways and those four are what license the relative route at larger n.
>   * PART (i) OF THE SOURCE CONJECTURE.  No m >= 3 graph is touched anywhere in this program.
>     The n = 4 numbers above are the m = 2 case of part (i) and nothing more; the first open
>     cell of part (i), m = 6 and n = 4, is not attempted.
>   * The GENERAL LEMMA for cubic triangle-free graphs is CHECKED at four control graphs, not
>     proved and not censused over cubic triangle-free graphs of any given order.
>   * THE HOMOTOPY TYPE.  This program computes homology and the 2-skeleton fact; the passage
>     to a wedge of spheres is Hurewicz plus Whitehead, a hand argument, not a computation.
>   * PRIOR ART.  Nothing here searches the literature.  The paper names the near miss and the
>     calibration in section 9 above is against a number that paper prints, but attribution is
>     a reading task and this program does not perform it.
