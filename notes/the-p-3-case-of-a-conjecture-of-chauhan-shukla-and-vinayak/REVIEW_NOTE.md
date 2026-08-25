# The p=3 Case of a Conjecture of Chauhan, Shukla, and Vinayak

`the-p-3-case-of-a-conjecture-of-chauhan-shukla-and-vinayak`

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
passes. The recorded run reports **138 checks, all passing**:

    VERDICT: ALL 138 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    aebe853ec86c85852da1923fb2af01bb2e2782d456042fd852c191acfdcd86d5

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT REPRODUCED
> 1. The paper's OWN two orders, L-hat_13 and L-hat_14, are not
> rebuilt here, and cannot be rebuilt from the manuscript: the
> base order L_n is specified only as 'the order of Chauhan and
> Shukla, Section 3 / Definition 3.10, specialized to p=3', the
> 169 and 238 triples are never printed, and the paper's 'Exact
> verification' section states that no data file and no copy of
> the program that produced the certificates accompanies it --
> that program and its transcript are held in an archive
> available on request, and are not inputs to this file.
> Consequently the audit-table entries
> 'exch. fail. 0' and 'purity fail. 0' FOR THE PAPER'S ORDERS, and
> the claim that S_13, S_14 are the spanning facets OF THOSE
> ORDERS, are not re-run.  What is re-run instead: the whole of
> the repair mechanism from the printed blocks (Section D), and
> the paper's CONCLUSION -- shellability of Delta_3(C_13^3) and
> Delta_3(C_14^3) with 1 and 8 spanning facets -- via orders this
> program builds and verifies itself (Sections C, E, G).
> 1a. Three further identities are NOT pinned, each verified by
> deliberate substitution to survive every check in this program:
> (i) the ORDER of the printed block is pinned only up to
> permutations with the same in-block discharge structure -- e.g.
> exchanging the two obstruction targets {3,5,12} and {4,5,12} in
> equation (3) still gives 15 of 15 and still discharges all four
> listed obstructions, so no check on the printed data can see it;
> (ii) the IDENTITY of the earlier witness triples T_i in the
> Section 2 table -- replacing {2,5,9} by {2,5,11} at n=13, j=112
> yields another facet with the same dead vertex 2 and a live
> vertex that the same reordering discharges, so D2-D5b cannot
> distinguish them; the table's claim that these are the ONLY
> obstruction pairs of L_n is likewise unverifiable without L_n;
> (iii) the members of S_13 and S_14 -- see 1b.
> All three are consequences of one omission: L_n is not printed.
> 1b. In particular the IDENTITY of the members of S_13 and S_14 is
> not pinned by this program.  Their NUMBER is forced (Sections C
> and E: h_{d+1} = |chi-tilde| = 1 and 8), every member is a
> genuine facet passing the order-independent spanning test (C10),
> and both sets are realised exactly by a shelling verified here
> (C11) -- but a different realisable set of the same size would
> also pass C10 and C11, because which facets are spanning depends
> on the order and the paper's order is not printed.
> 2. The census stops at n = 32.  Theorem 1 claims every n >= 13;
> n > 32 is not machine-checked here (the paper obtains it from
> [CS, Theorem 1.4], not from computation).  The Euler-
> characteristic check stops at n = 22 and the GF(2) homology at
> n = 20, both for cost: they enumerate all 2^n vertex subsets.
> 3. Integral ODD torsion in degrees below the top is not excluded by
> the GF(2) computation of Section E; it is excluded instead by
> shellability plus [Kozlov, Thm 12.3], which Sections C and G
> establish computationally for the two boundary cases.
