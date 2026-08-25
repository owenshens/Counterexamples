# A Negative Answer to Mohammadi's Question on Powers of Cover Ideals

`eight-vertex-graph-disproves-mohammadi-cover-ideal-powers-question`

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
passes. The recorded run reports **123 checks, all passing**:

    VERDICT: ALL 123 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    9bf539e18fc331bfd8a090233da71fdd2fe550894c79bee3a20cca48ea148dfc

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> --- what this program does NOT establish (stated, not hidden)
> 1. The wording of Mohammadi's Question 4.1 (Collect. Math. 65 (2014), paywalled, no preprint) and of the Ha-Van Tuyl Question 6.2 restatement quoted in the paper are TEXTUAL claims; no program can check them.  What is checked here is that the mathematical object satisfies the hypothesis (CM, 5-linear cover ideal) and violates the conclusion (square not linear).
> 2. \cite[Theorem~1.34]{MillerSturmfels} is used as the definition of the upper-Koszul engine in sections 4 and 7.  It is NOT assumed for the headline: beta_{1,a}(J(G)^2)=1 is decided again in section 6 by syzygy linear algebra and once more field-free, and section 7b shows the two engines agree on beta_1(J(G)).  Higher Betti numbers of J(G) (beta_{2,7}=4) still rest on Theorem 1.34 as implemented.
> 3. 'Cohen-Macaulay over every field' is carried by SHELLABILITY, which is field-free and is checked in section 2; the step shellable => Cohen-Macaulay, and Eagon-Reiner (CM <=> linear resolution of the dual), are cited theorems, not verified here.  The six-field Betti table of section 7 corroborates but cannot by itself reach 'every field'; section 7b closes the beta_1 column field-free.
> 4. The Remark's comparison with Ficarra-Moradi (their Theorem D, the hitting set {a,c,f,y_1,...,y_{n-6}}, non-quadratic Alexander dual) is about THEIR construction and is not reproduced here.  Only the two claims about G itself are checked: the induced 4-cycle 1-6-5-7-1 and gcd of the generators = 1.
> 5. Only beta_0 and beta_1 are computed for J(G)^2.  The paper needs no more: generation in degree 10 plus one first syzygy in degree 12 already contradicts linearity.  Nothing here computes the full resolution of J(G)^2, and nothing needs to.
