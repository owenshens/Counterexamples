# A Counterexample to the Aggregate-Condition Form of Karp–Sitnik Conjecture 1

`integer-parameter-3f2-disproves-karp-sitnik-aggregate-condition-conjecture`

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
passes. The recorded run reports **30 checks, all passing**:

    VERDICT: ALL 30 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    a2807bd5c9149cf306490caadfa1fc4b1182a04b2d8dc067c864f34b4f21293c

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the READING of Conjecture 1 that makes the object above a counterexample. This program decides inequalities about hypergeometric series; it cannot decide what the authors of [KS] meant. What is refuted above is the form of Conjecture 1 whose parameter hypothesis is sum(b_i-a_i)>0 together with positivity, i.e. the form in which the aggregate condition REPLACES b_i>a_i. It is NOT a refutation of the rival reading 'b_i>=a_i>0 together with sum(b_i-a_i)>0', under which the aggregate clause is not idle (it excludes only A=B) and under which the exhibited object is void, since b_1=1>=a_1=2 fails; the two checks exhibited-object-fails-b_i>=a_i-too... and bound-holds-on-non-strict-componentwise-parameters... record that shortfall as computed facts, and no counterexample to that rival reading is claimed here or in the paper.
> NOT RE-RUN: the quoted source text. No external table, catalogue, preprint or journal page is consulted by this program, and none is needed for the inequality; but the paper's contribution is a reading of other authors' text, so none of the following is machine-checked here: that Conjecture 1 of [KS] reads verbatim as printed at the head of this transcript; that Conjecture 2 of [KS] lists exactly the hypotheses quoted in the paper, with no componentwise clause; that those numbers are the ones carried by the arXiv version cited; that Theorem 3 of [KS], the bound credited there to Luke, and Theorem 5 of [KP] have the hypotheses assumed above; and that [Derbazi] is the preprint cited, of the version cited, containing the quoted sentence about the broader conjectural condition at the place cited. A referee must compare every quotation, arXiv identifier and statement number in the paper against the sources by eye.
