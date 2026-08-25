# A Counterexample to the Etzion–Vardy–Yaakobi Conference-Matrix Conjecture

`a-counterexample-to-the-etzion-vardy-yaakobi-conference-matrix-conjecture`

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
passes. The recorded run reports **43 checks, all passing**:

    VERDICT: ALL 43 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    115bf26af674b922c35f128783603804025a3b45b11120ac93ed33cc4d14f5dc

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE NOT RE-RUN: the exact minimum distance of C. The paper claims only d(C) <= 13 and that claim is fully verified above. The census settles the minimum weight over 52531811444 of the 8629188747598184440948 nonzero codewords (information weight <= 5); codewords of information weight >= 6 are NOT enumerated, because depth 6 costs about 46 times depth 5.
> NOTE no exhaustive census over conference matrices of order 30 is attempted: the paper exhibits one matrix and claims nothing about the others.
> NOTE NOT VERIFIED: the provenance sentence (representative 3 of Spence's order-30 conference two-graph file). This program is standard library only and fetches nothing; the paper states that the provenance is not used in the proof, and the checks above are self-contained.
> NOTE MACHINERY SELF-TESTS (4 of the 43 checks; they test this program, not the paper, and cannot fail for a well-formed input): number_of_off_diagonal_pairs_is_435; census_enumeration_is_complete_for_its_range; census_witness_weight_equals_the_census_minimum; census_witness_lies_in_the_F_29_row_span_of_W. Each is a count checked against a closed form, or a second route that must agree with a first. A loop-bound or elimination error here reports FAIL, but a wrong matrix is caught by the Gram, minor, rank and codeword checks, not by these.
> NOTE DEPENDENT RESTATEMENTS (9 of the 43 checks): each of the following is a logical consequence of an earlier PASS line, so, given that line, it cannot fail on its own -- conference_constant_n_minus_1_equals_p (follows from p_is_prime_and_order_equals_p_plus_1, which already forces n = p+1); minor_is_nonzero_mod_p_so_rank_at_least_15 (follows from minor_residue_mod_p_equals_21, and 21 is not 0); dimension_is_half_the_length (follows from rank_over_F_29_equals_15 with order_of_the_matrix_is_30); code_is_euclidean_self_dual (follows from every_C_perp_basis_vector_lies_in_C, whose condition it repeats); xW_is_a_nonzero_codeword (follows from xW_has_hamming_weight_13); xW_information_part_is_the_first_unit_vector (follows from codeword_re_derived_from_rref_without_using_x with pivot_columns_are_the_first_15); derived_dimension_agrees_with_the_conjecture (follows from conjectured_parameters_p1_p1half_p3half_are_30_15_16); code_is_not_MDS (follows from minimum_distance_upper_bound_beats_the_conjectured_distance, the Singleton value having been checked equal to the conjectured distance); census_minimum_is_at_most_the_papers_bound_13 (follows from census_minimum_equals_12_derived_here, and 12 <= 13). They are printed because each names the step of the argument it licenses, but the number of lines carrying force independent of every other line is 30, not 43: a referee should read the VERDICT count below as 30 independent checks.
> NOTE KNOWN INSENSITIVITY: the checks are invariant under the switchings D W D, D = diag(+-1), whose -1 entries are confined to the coordinates where x and xW both vanish -- computed from the two vectors printed above, those are the 5 coordinates 5, 9, 18, 21 and 24, an orbit of 2^5 = 32 displays. Each orbit member is itself a conference matrix whose F_29 row span is a self-dual [30,15] code containing the same weight-13 word, so the theorem holds verbatim for every one of them; only the provenance sentence, which the paper excludes from the proof, distinguishes them. The decoded matrix is printed back above so the display itself can be compared character by character.
