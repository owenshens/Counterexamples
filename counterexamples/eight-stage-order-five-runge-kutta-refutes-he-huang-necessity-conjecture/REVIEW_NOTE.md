# A Counterexample to He and Huang's Runge–Kutta\\ Necessity Conjecture

`eight-stage-order-five-runge-kutta-refutes-he-huang-necessity-conjecture`

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
passes. The recorded run reports **47 checks, all passing**:

    VERDICT: ALL 47 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    222c9e5d5688dcd31da1a91f2dbe096070c27ee758c93e92a533b0fb3442e6a9

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> GAPS -- steps NOT decided by this program (read before trusting it)
> G1 SOURCE FIDELITY.  q_r = A c^{or} - c^{o(r+1)}/(r+1), the generator list
> Q_2 = span{q_1, A q_0, q_0 o c, q_0}, the chain Q_1 <= Q_2 <= ..., the reading
> of QO(m) as the Hadamard b o Q_m = {0}, and the admissibility constraints
> m >= n-1, m+n+1 >= p are all TRANSCRIBED from the paper's quotation of
> He-Huang.  No line of this program can detect a misquotation.  If any of
> them is misquoted, every check above still passes and the paper is wrong.
> G2 Q_2 <= Q_m FOR m >= 2.  The witness q_1 is shown to lie in the printed
> spanning set of Q_2 and to satisfy b o q_1 != 0.  That it also lies in
> Q_m for every m >= 2 is taken from the quoted chain; this program never
> constructs Q_m for m >= 3, so non-monotonicity of the He-Huang spaces
> would defeat the 'for every admissible (m,n)' step while leaving the
> QO(2) computation intact.
> G3 INTEGRALITY OF (m,n).  min m = 2 needs m to be an integer; section 12
> shows the half-integer relaxation gives 3/2.  The paper attributes
> 'nonnegative integers m,n' to Theorem 1.2, while the source is reported
> to write only 'for some m, n'.  This program verifies the claim UNDER the
> paper's stated hypothesis; it cannot verify the hypothesis.
> G4 HADAMARD VERSUS SCALAR.  The refutation is specific to the componentwise
> reading.  Section 11 shows b_hat^T q_1 = 0 EXACTLY, so under a scalar
> reading of QO(m) nothing here refutes anything.  The paper says so; a
> reader who thinks He-Huang meant the scalar product should reject the
> paper, not this program.
> G5 ATTRIBUTION TO VERNER.  That this tableau is Table 4 of Verner (2014),
> RK(8-6:5)a, is NOT checked and cannot be checked here.  It is also not
> needed: the order-five property is derived from the tableau itself, so a
> misattribution would be a citation defect, not a mathematical one.
> G6 THE OTHER FOUR CONDITIONS.  DO(n), QD_weak(m,n), PR(n) and QR(m)
> are never evaluated.  This is sound -- refuting one member of a
> conjunction refutes the conjunction -- but it means the program says
> nothing about whether those four members of the Theorem 1.2 package
> hold.
> B(p) is NOT among them: it IS evaluated, at the derived p = 5, both
> directly (b_hat^T c^{o(k-1)} = 1/k for k = 1..5, all 5 components
> verified above) and again by the bush trees of the order-1..5 sweep,
> and it HOLDS on this tableau.  A method may of course satisfy B(p)
> and still fail QO(m), which is exactly what happens here.
> G7 NO PRIOR-ART CHECK.  Novelty of the counterexample is a literature
> question and is outside the scope of any arithmetic.
