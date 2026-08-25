# An Elementary Proof of Erlbacher's Counterexample to TxGraffiti Conjecture 9

`fourteen-vertex-graph-refutes-davila-conjecture-9`

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
passes. The recorded run reports **35 checks, all passing**:

    VERDICT: ALL 35 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    3721bb78c1301ac5cb966b4d2852cb131081386c085db95c85f798c07d0d110a

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE DISCLOSURE: the chain family as described in the paper is under-determined.  "two edges of each middle block" does not say WHICH two, and up to Aut(K_{3,3}) there are two readings: an independent pair and a pair sharing an endpoint.  Both give connected cubic triangle-free graphs on 8k-2 vertices, and the two readings agree at k=2, which has no middle block, so the 14-vertex graph and its gamma=4, Z=7 are unaffected.  They can differ once a middle block exists (k >= 3), and at k=4 they do: the excesses 4 at n=30 and 5 at n=38 checked above hold for the independent-pair reading, while the sharing-pair reading gives gamma=10 and Z<=13 at n=30, i.e. excess at most 3.  A reader who reconstructs the family the other way will not reproduce the paper's 4 at n=30.
> NOT RE-RUN: (i) none of the attributions the paper's framing rests on can be checked by an offline program, and none of them is checked here -- that the bound Z <= gamma + 2 for connected cubic diamond-free graphs is Conjecture 9, under that number, of arXiv:2406.19231v2; that Zenodo DOI 10.5281/zenodo.21269439 resolves at version 3 of July 8, 2026; and that the sub-path cited in the bibliography exists inside that record: a reader must check all three at the sources, and this program's PASS lines say nothing about them; (ii) the two independent implementations, the mutation tests and the independent re-verification reported in the cited artifact are not re-executed -- every value above is instead recomputed here from the edge rule; (iii) the chain family is described in the paper only up to an ambiguity: it does not say which two edges of a middle block are subdivided, so the excesses 4 at n=30 and 5 at n=38 are verified here for the independent-pair reading alone, and under the other reading (a pair sharing an endpoint) this program computes gamma=10 and Z<=13 at n=30, an excess of at most 3, so the paper's reported 4 at n=30 holds under one of the two readings and fails under the other; (iv) Z is determined here for 2 <= k <= 5 of the chain family and only gamma for 6 <= k <= 8, which is the range the paper itself reports.
