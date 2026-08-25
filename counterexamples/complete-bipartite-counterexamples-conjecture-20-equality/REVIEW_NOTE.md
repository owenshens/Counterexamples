# Counterexamples to the Equality Clause in Conjecture 20 of Chen–Guo–Li–Wang

`complete-bipartite-counterexamples-conjecture-20-equality`

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
passes. The recorded run reports **20 checks, all passing**:

    VERDICT: ALL 20 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    19638bfb24ded3d13d48a20d53a3ebe22446d1b55764e1fb5bd18eb0c34271cc

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: the exhaustive census over n = 7, 8, 9 (no isomorph-free generation was performed above n = 6 and nauty/geng was not used), hence the published tally 'exactly five cases, all at n = 8, k = 6' is NOT confirmed as exhaustive, and the five named graphs are not checked to be pairwise non-isomorphic either (only that GCXnf_ is self-complementary, that G?zvf_ and GQhTUg are complementary, and that G?zvf_ and GCXnf_ differ), so the count five is not confirmed from below here either; the minimality claim that K_{4,4} has smallest order among counterexamples is confirmed only against orders n <= 6; the product inequality mu_k(G) mu_k(complement G) <= n(n-k) itself is nowhere proved here, only evaluated on the graphs enumerated or named above (and for 1 <= k < n/2 it is untouched by the family, which lives at k = 3n/4); the family K_{2t,2t} is verified for t = 2..5 only, the general t being closed-form; the transcription of Conjecture 20, of Theorem 18 and of Remark 1 from the published article was not checked against the published text (no external source was consulted, so the equality clause is taken exactly as quoted in the note -- the refutation is however also checked to survive the weakest reading, at least 2 components instead of at least k+1); and the remark that Chen, Guo, Li and Wang had verified the conjecture for all graphs on at most nine vertices is reported, not examined, here.
