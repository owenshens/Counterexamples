# Unimodal unicyclic graphs with several breaks of log-concavity

`an-affirmative-answer-to-the-levit-kadrawi-multi-break-unicyclic-closure-problem`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |
| `REVIEW_NOTE.md` | this file |

That is the whole folder: five files, no subdirectory and no data file.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
About three seconds on a laptop. The program prints one line per check and a closing verdict,
and exits 0 only if every check passes. The recorded run reports **56 checks, all passing**:

    VERDICT: ALL 56 CHECKS PASS

It does not read an exhibited object and confirm it: it RE-DERIVES the paper's claim.
Every graph is rebuilt from the pattern-graph parameters and the labelling rule printed in
section 2 --- nothing about a graph is taken from the paper except its parameters --- and every
integer the paper prints --- the orders, the independence numbers, the modes, the full 32-term coefficient
list of `H*`, the printed coefficient tails of the other five witnesses, all 17 printed
products `i_k^2` and `i_{k-1}i_{k+1}`, and the break sets in full --- is recomputed and
compared. The independence polynomial of each of the six unicyclic witnesses is computed twice,
by the two independent identities of Lemma 3 of the paper (delete a cycle edge; split at a
cycle vertex), and both are required to agree. The program also runs, beyond what the paper
claims, complete censuses of the two cells `T_{G,2,6}` and `T*` (3,240 and 1,596 closures);
the paper states no census. The `graph6` string printed in section 3 is
decoded and matched edge for edge against the constructed `H*`.

The program also carries its own controls, and they run in both polarities. Negative: the
forest recursion and both unicyclic identities are checked against brute-force subset
enumeration on all 1,441 trees and all 3,898 unicyclic graphs of order at most 6; and the
break detector must stay silent on 76 log-concave graphs, the paths and cycles of order up to
40, whose coefficients are simultaneously matched against the closed forms
`i_k(P_n) = C(n-k+1,k)` and `i_k(C_n) = (n/(n-k))C(n-k,k)`. Positive: published integers of
other authors are reproduced --- `GT_{2,5}` with degree 37 and breaks {34,36}, the laws
`alpha = 10m+4` with breaks {10m+2,10m+3} and `alpha = 12m+8` with breaks
{12m+5,12m+6,12m+7} for the two named families, and the order 26 at which trees first fail
log-concavity, attained by `T_{3,4,4}` with nothing smaller breaking --- this last inside the
family `T_{3,m,n}`, `1 <= m,n <= 8`, which is the search the program performs, not over all
trees.

The controls also run in the failing direction. Corrupting the exhibited data in a scratch copy
of the program --- one digit of `H*`'s printed coefficient list, one character of the `graph6`
string, one of the 17 printed products (specifically, restoring the pre-review value
`48867289141230114` for `i_39 i_41` of `H_1`), and one index added to `H*`'s break set ---
makes it print `FAIL` and exit 1 in each case, one failed check out of 56 each time. The
program is therefore checking the claims and not merely printing them.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    657b5af756aafcc731a13b5c9e0c2945a0dfdb45c7a8c4cfdcb5a77fa8fffa50

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: this program checks the paper, and the paper only. It does NOT enumerate the
> nonedges of the four cells other than T_{G,2,6} and T* -- in particular no census of
> (T_{1,3}:S_{2,4})_2^(16) or of (T_{1,7}:S_{2,5})_2^(13) is performed, and none is claimed;
> H_A and H_B are checked as the two individual graphs they are. It does NOT search for
> smaller witnesses, so no minimality claim is supported. It does NOT verify any statement
> about arbitrarily many breaks: the largest number of breaks it exhibits is 7, and the
> m = 1..7 table is a finite computation, not a proof of a law. It does NOT re-read the source
> e-print arXiv:2603.17114, so the quotations, line numbers and byte counts in section 1 of the
> paper are outside its reach and must be checked against the e-print.

Three further limits belong here rather than in the program, because they are about the reading
of the problem and of the literature rather than about arithmetic.

First, on what is and is not new. The paper is not a claim of unqualified novelty and should not
be read as one. The unicyclic graph `H_1` is one of the 3,240 enlargements of `T_{G,2,6}` that
the authors of the source problem computed themselves; what the paper records for it is the
number and position of the breaks. On the tree side the multiplicity is already a published
theorem (Bautista-Ramos, arXiv:2511.00334, Figure 1 caption), and the two- and
three-consecutive-break tree families are those of arXiv:2603.14204 (Graphs and Combinatorics,
doi:10.1007/s00373-026-03054-4). Section 1 of the paper disclaims any tree-side result.

Second, the trees underlying `H*` and `H**` are
pattern graphs in the sense of [2] but are not rows of that paper's table of named families;
under the strict reading of "the newer pattern-graph families" the consecutive-break clause is
carried by `H_A` and `H_B`, which do lie in those families. The paper's Remarks record this.
Third, the closure preserves the tree's break set in five of
the six witnesses; amplification is exhibited only by `H_7`, and the paper's abstract and
Remarks say so.

Fourth, the source statement itself. The program does not read arXiv:2603.17114, so the
paper's quotation of Problem 6.6 is not covered by the recorded run; the paper says as much,
and its theorem does not depend on the reading of that problem.
