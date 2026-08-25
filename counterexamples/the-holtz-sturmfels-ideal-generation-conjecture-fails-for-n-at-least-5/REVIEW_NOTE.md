# The Holtz–Sturmfels Ideal-Generation Conjecture Fails for \(n 5\)

`the-holtz-sturmfels-ideal-generation-conjecture-fails-for-n-at-least-5`

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

    cebc860b3d14a2dbdba62a944bb217871f880f63945ae21163a02e9119c207ec

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the paper asserts J_n != P_n for EVERY n >= 5, an infinite family.
> Re-run here: n = 5 in full, exactly and symbolically (the exhibited object);
> the relations F_i^S for n = 6..12, over every one of the 1715 5-subsets S and
> every i, at exact integer points; and the HD(n) module facts (dimension and
> vanishing of H_S-invariants) for n = 3..14 over every (3-subset, 5-subset) pair.
> NOT re-run: every n > 14, and the relations for n = 13..14. Those rest on the
> paper's uniform argument (restriction of phi_n to a 5-subset), not on computation.
> NOT re-run: the descent of the Status paragraph to the REAL form of Conjecture 14
> (ring R[A_.], group SL_2(R)^n).  Supplied here: every coefficient of the pentad, the
> B_i, the F_i and the 283 constructed and sampled generators of HD(5)^{weight 0} is an
> integer (96648 coefficients scanned, 0 non-integral), and every rank and kernel above is
> computed over Q by exact elimination, so it is the rank of the same integer matrix
> over R.  NOT supplied: that HD(n) itself has a basis of quartics with integer
> coefficients.  That is taken from the paper (orbit span of an integral quartic under
> a group defined over Q) and not re-derived, so the real form is not re-verified here.
> --- GROUP 5: negative controls (each corruption must be detected) ---
> control pentad_term_sign_flipped               detected=True
> control pentad_term_deleted                    detected=True
> control pullback_constant_12_to_13             detected=True
> control pullback_t_exponent_10_to_9            detected=True
> control epsilon_made_symmetric                 detected=True
> control one_B_2_coefficient_perturbed          detected=True
> control lemma3_matrix_entry_perturbed          detected=True
> control Sym4_invariant_dim_forced_to_1         detected=True
> control hamiltonian_cycle_signs_dropped        detected=True
> control one_pentad_variable_substituted        detected=True
> control rank_test_admits_HD5_elements_and_rejects_F_2 detected=True
> control saturation_rank1_test_rejects_a_non_member detected=True
> control perturbed_cayley_slice_no_longer_pulls_back_to_zero detected=True
