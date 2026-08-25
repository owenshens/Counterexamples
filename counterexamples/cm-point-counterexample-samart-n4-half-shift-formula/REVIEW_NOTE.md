# A Counterexample to Samart's Half-Shift Formula for n_4

`cm-point-counterexample-samart-n4-half-shift-formula`

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
passes. The recorded run reports **51 checks, all passing**:

    VERDICT: ALL 51 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    a400331eee18a94a2f183746228e07fed48b769a99646a56d88552cef13e5df4

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: the paper's own run used a multiprecision interval library; this program reproduces the certificate in binary64 and confirms on a sample of 190 of the 262144 cells at 60 decimal digits that the floors do not depend on that. NO ARM OF THIS PROGRAM IS AN ENCLOSURE: the trusted computing base the paper's Reproducibility paragraph pins down -- an interval context at 80 decimal digits, with each binary64 result enlarged outward by nextafter -- is not exercised here, and neither is directed interval arithmetic for the per-cell lower bounds on |a_{c*} +- 2|; the 512^2 census is evaluated at binary64 point values and the 60-digit arm at Decimal point values. What stands in for those semantics is the pair census_binary64_not_deciding and floor_margin_beats_binary64_budget, the latter showing that the census's own minimal distance from 4096*arcosh(X) to an integer exceeds a 64-ulp binary64 error budget (2.98e-08) by a factor of 50.8. That is a budget argument about the floors, not a rigorous enclosure, so a referee wanting the interval semantics themselves must still consult the paper's own run. The values of n_4 used for the consistency and scan checks are midpoint quadratures, not enclosures; the certified bound 151680651/67108864 is the exact rational one and is what the refutation rests on.
> NOT RE-RUN HERE (covering hypothesis, sampling): the two checks that touch sup|a_c - a_{c*}| < delta empirically are SUBSAMPLES, not cell-by-cell verifications. cell_covering_within_delta visits 1600 of the 262144 cells (every 13th index in each direction) and evaluates a 9x9 sub-lattice inside each; cell_lower_bound_end_to_end visits 5476 of the 262144 cells (every 7th index) and takes a 7x7 minimum inside each. A sampled maximum is not a supremum and a sampled minimum is not an infimum, so neither check would detect a violation confined to an unvisited cell or to a point between sample nodes; they corroborate the hypothesis, they do not establish it. What establishes it is the exact-rational pair derivative_bound_arithmetic (|grad a_c| <= 2 + |c|/2 = 263/100 on the whole torus, from the triangle inequality) and delta_covers_cell_radius (263*(355/113)/25600 + (99/70)/200 < 123/3125 with pi and sqrt2 bracketed rationally), which together bound the deviation over every cell without sampling; derivative_bound_sampled is likewise only a 57600-point corroboration of that analytic bound. The census lower bound therefore rests on the analytic covering argument, not on these samples.
> OBSERVED HERE, AND NARROWER THAN THE PAPER'S CONCESSION TO EARLIER WORK: the paper concedes that the blanket statement it quotes from an earlier status report, that below Im tau = 1/sqrt2 the identity fails everywhere, already covers tau0. This program's scan does not support that statement on the half-shift line Re tau = 1/2: at tau = 1/2+0.7000i, with y^2 = 0.490000 < 1/2, the two sides agree to 3.52e-12, while the largest y on that line at which the scan exhibits a failure is 0.6800 (diff 3.9e-03); the earlier evidence was gathered on Re tau = 0 and Re tau = 1/4. See scan_holds_below_inv_sqrt2_on_half_line. This makes the paper's contribution larger rather than smaller and does not touch Theorem 1, whose refutation rests on the exact rational certificate at tau0; the scan values are midpoint quadratures, not enclosures.
