# Four Mutually Unbiased Bases in Dimension Eight Give a Non-Decomposable Witness

`four-mubs-detect-bound-entanglement-in-dimension-eight-where-five-were-conjectured`

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
About a fifth of a second, single threaded. Every comparison it makes is between exact
integers or exact `fractions.Fraction` objects; floating point appears only inside
informational strings. The program prints one line per check and a closing verdict, and exits
0 only if every check passes. The recorded run reports **23 checks, all passing**:

    VERDICT: ALL 23 CHECKS PASS

It reads the object exhibited in the paper as input and nothing else: the index conventions
(2.1) and (2.2), the three symmetric matrices over `F_2`, the shift `s = 1`, the witness
(2.3), and the 288 integers of Table 1, pasted into the program verbatim. From these it
rebuilds the four bases, checks orthonormality and the 384 cross overlaps, assembles
`64 W^Gamma(M_4,1)` in `Z[i]`, verifies the projector identities and the kernel dimension 34
behind `lambda_max = 11/8`, parses Table 1, factors both `rho` and `rho^Gamma` by exact
rational `LDL^T`, and recomputes `tr[W^Gamma(M_4,1) rho] = -14588/524801` by two independent
routes. It carries a control of each polarity: `tr[W^Gamma . I_64/64] = 7/8 > 0` on a PPT
state, and `z^T W^Gamma z = -1/4 < 0` on the vector `z = e_1 + e_19`.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    e6462c241ab4d4f0a59e160f92ca6760138a50a17aa7e1d90d399442d37bd986

## Scope

**What kind of check this is.** The program re-derives every computational claim the
paper makes *about the exhibited object* — it rebuilds the bases from the printed
conventions rather than reading them, and recomputes the MUB property, the witness, the
projector and kernel facts, both `LDL^T` factorisations and the objective — and it then
takes the one deductive step from those to the paper's conclusion: `rho` is PPT and
`tr[W^Gamma(M_4,1) rho] < 0`, hence the witness is non-decomposable, hence at `d = 8` the
minimal number of MUBs is at most `m = 4 < 5` and the conjecture is false at `r = 3`.
What it does **not** do is re-derive the *existence* statement from scratch: it neither
searches for the object nor re-runs the numerical optimisation that produced `rho`. It is
therefore a checker of the exhibited object plus the deduction from it, not an
independent rediscovery — which is what the paper claims for it (§6) and no more. Its
inputs are exactly the objects printed in the paper, so a referee who distrusts the
program can redo the same work by hand from §2 and §3 alone.

Adversarially: corrupting the exhibited object makes the program fail. Changing one
integer of Table 1 by one, flipping the sign of one off-diagonal entry, substituting a
different stabilizer basis, or altering one printed `S` matrix each produce `FAIL` lines
and exit status 1.

The program's own statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the search that located the exhibited object: the reduction of the stabilizer
> MUB sets of C^8 to 1045 configurations and the 8360 non-decomposability decisions over all
> 8 shifts are NOT re-run here; only the printed object is.
>
> NOT RE-RUN: the numerical optimisation that produced rho is NOT re-run; rho is read from
> Table 1 of the paper and only re-verified.
>
> NOT RE-RUN: the exact minimum of tr[W^Gamma(M_4,1) sigma] over PPT sigma is NOT determined;
> the numerical upper bound -0.0350501 quoted in the paper is NOT certified here, and only the
> SIGN of -14588/524801 is used.
>
> NOT RE-RUN: lambda_min(W^Gamma) = -13/8 is NOT certified; only lambda_max = 11/8 is (via the
> rank of A + Q_0 + Q_13 + Q_19), and neither is needed for the result.
>
> NOT RE-RUN: m = 3 at d = 8, which would widen the margin to two, and the m = 5 = d/2+1
> sufficiency question at d = 8, are NOT addressed here.
>
> NOT RE-RUN: dimensions d = 2^r with r >= 4 are NOT addressed; nothing here is run at d = 16.
>
> NOT RE-RUN: MUB families of C^8 outside the stabilizer class (2.1) are NOT considered, and no
> claim is made that every 4-MUB set or every shift gives a non-decomposable witness.
>
> NOT RE-RUN: that W^Gamma(M_m,s) is nonnegative on separable states is quoted from Spengler
> et al., Phys. Rev. A 86, 022311 (2012); it is NOT reproved here.
>
> NOT RE-RUN: the d = 4 (r = 2) refutation mentioned in the paper's Status paragraph is NOT
> reproduced here.
