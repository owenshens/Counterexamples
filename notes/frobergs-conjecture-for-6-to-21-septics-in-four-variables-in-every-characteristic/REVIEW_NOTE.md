# A Characteristic-Free Integral Witness for Fröberg's Series for Six to Twenty-One Septics in Four Variables

`frobergs-conjecture-for-6-to-21-septics-in-four-variables-in-every-characteristic`

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
passes. The recorded run reports **222 checks, all passing**:

    VERDICT: ALL 222 CHECKS PASS

It holds the objects the paper exhibits — the 120 hexadecimal masks of the family, of which the
paper's Appendix A prints the 21 the theorem uses, the ten endpoint cells of Table 1, and the
twenty index sets and twenty integer determinants, whose bit lengths and gcds are the last two
columns of the paper's Table 1 — and re-derives every quantity
the paper claims: it decodes the family and
recomputes its `sha256`; re-derives `q(r)` for `r = 6..21`, the five blocks, and the claim that
the ten printed cells are exactly the endpoint cells; rebuilds each of the ten integer
matrices from the decoded forms, extracts the twenty named minors, and recomputes each of the
twenty determinants from scratch by Chinese remaindering against the integer Hadamard bound;
takes the ten gcds; and, independently of the gcd argument, recomputes the rank of all ten
matrices over `F_2`, `F_3` and `F_5`. All arithmetic is exact integer arithmetic — no floating
point enters any decision. The run took 168 s on 30 cores.

The recorded run is a run on a 32-vCPU Linux host rather than on a laptop, because the twenty
determinants cost about 2,500 core-seconds; `multiprocessing` is used only to spread
independent modular determinants over cores, and the values do not depend on the number of
workers.

Both polarities are controlled, each control named to the object it tests: two closed-form
determinants (`det(J_n - I_n) = (-1)^(n-1)(n-1)`, `n = 12, 30`) and a 40×40 block of a real
minor computed twice, by Chinese remaindering and by fraction-free Bareiss elimination, as
checks on the determinant engines; a forced-positive object (the monomial complete
intersection `r = 4`, `F_i = x_i^7`, degree 13, whose 336×336 minor must have determinant
`±1`); a provably deficient object (`r = 6` with the tuple `(x1^7, x2^7, x3^7, x4^7, x1^7,
x2^7)`, degree 13, whose rank is 336 and not the maximal 504, so the decider is able to say
no); and two tamper controls — flipping one bit of one printed mask changes the determinant it
feeds, and a determinant off by one fails the residue test — so the comparisons are not
vacuous.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    7aed8da9d5c7e949dd08fea339271119ce7aef99ccf9261a9301653f2fb26385

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN, and the paper says so in its Section 7 and in REVIEW_NOTE.md's Scope: (1) the
> identity of the printed family with the ancillary file published with arXiv:2608.24797 --
> the sha256 above is recomputed from the printed masks, but the external file is not fetched,
> since this program reads no file and no network; (2) r <= 5 and r >= 22, and every (n,d)
> other than (4,7) -- no claim is made and none is checked; (3) the 168-prime sweep reported
> in the source note: only p = 2, 3, 5 are re-run here as an independent confirmation, and the
> gcd argument, not that sweep, is what covers every prime; (4) the minor-selection heuristic
> -- which s x s minors to take is not re-derived, and it does not need to be, because the
> printed index sets are re-extracted and their determinants recomputed from scratch; (5) any
> claim about GENERIC forms beyond the specialisation argument of Section 2, which is a proof
> and not a computation.

Two further limits belong to the paper rather than to the program, and Section 5 states both.
The interval is `6 <= r <= 21` and nothing wider: nothing is claimed for `r >= 22`, and nothing
for `r = 5`. Note also that at `r = 21` the predicted series has support through degree 9 and
vanishes from degree 10 on; what the conjecture requires there is surjectivity in degree 10.
The program's Step 2 check for that cell is named
`r_21_still_has_positive_series_support_in_degree_10`; the name is inaccurate — the values it
prints, `a_9 = 10` and `a_10 = -134`, are the correct ones, and the paper states the correct
reading.
And what is added is narrower than the theorem: inside `6 <= r <= 21` the characteristic-2
slice already follows from the rank table over `F_2` published in arXiv:2608.24797 together
with the propagation lemma of that same paper, which carries no field hypothesis. What is
added here is the coprime-pair certificate, which covers the remaining primes, and the fact
that the witness is a single explicit tuple over the prime field `F_p`. The family of 120
forms, the block table and the propagation lemma are that paper's, not ours; the lemma is
proved in Section 4 rather than cited so that the argument is self-contained.
