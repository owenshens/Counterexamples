# A Degree-18 Counterexample to Conjecture 5.23 of Acevedo–Blekherman–Debus–Riener

`degree-18-vector-disproves-acevedo-blekherman-debus-riener-conjecture`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |
| `paper_appendix_listing.output.txt` | a recorded run of `paper_appendix_listing.py` |
| `paper_appendix_listing.py` | A self-contained standard-library program, identical to the listing printed in the paper's appendix, that regenerates the even partitions of 18 together with the partial-symmetry blocks and inequalities, prints their counts and the minimum first-order and partial-symmetry slacks, and exhibits the unique second-order inequality violated by the vector in question. |
| `search_trop_vs_T2_small_degrees.output.txt` | a recorded run of `search_trop_vs_T2_small_degrees.py` |
| `search_trop_vs_T2_small_degrees.py` | An exact rational-arithmetic search that, for degrees 2d = 6, 8, 10 and 12, tests every generator of the second-order superdominance system for membership in the cone of partial-symmetry inequalities describing the tropicalized sums-of-squares cone, returning a Farkas certificate for any generator that fails. |

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **19 checks, all passing**:

    VERDICT: ALL 19 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    0542cd17ce2613ee2971068ab588b1f00c362506a488a3314c800cd2806645ff

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT VERIFIED HERE: (i) the membership Lemma (see the TRUST BOUNDARY at
> C12) -- no Python program can check that reduction; (ii) the witness
> table, the degree and the definitions, which are transcribed INPUTS
> (P1-P4) -- a transcription error is invisible to this program; (iii)
> the ordered-alpha reading of note N1, which was CHOSEN because it
> reproduces the paper's 131, so the counts cannot discriminate against
> a paper that read alpha the same way.  Every other computational
> assertion of the paper's 'Exact verification' paragraph -- the five
> counts, the 417 first-order and 1056 partial-symmetry slacks, the
> second-order violation, the 7025-inequality T^(2) system with exactly
> one violated inequality, and the two further facts above -- IS derived
> and compared here.
> --- C0: coordinate set Lambda^{ev}_18 ---
> derived |EvenPartitions(18)| = 30
> derived |EvenPartitions(0)| = 1  ->  ()
> derived |EvenPartitions(2)| = 1  ->  (2)
> derived |EvenPartitions(4)| = 2  ->  (4) (2,2)
> derived |EvenPartitions(6)| = 3  ->  (6) (4,2) (2,2,2)
> derived |EvenPartitions(8)| = 5  ->  (8) (6,2) (4,4) (4,2,2) (2,2,2,2)
> INDEPENDENT DP: #even partitions of 18 = p(9) = 30
