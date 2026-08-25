# A Counterexample to the Unrestricted Dargad–Larsson Domination Problem

`residue-b-domination-refutes-unrestricted-dargad-larsson-problem`

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
passes. The recorded run reports **31 checks, all passing**:

    VERDICT: ALL 31 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    fbcb362333b58849ca0ece0c80d459da6f39acbec5f0caa0904c0765af8fc6b1

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> scope note: the sweeps above are exhaustive over the stated finite ranges of (a,b) (b <= 20 for the order,
> b <= 10 for periodicity and least period, b < 20 for the boundary dichotomy, b <= 16 for the derivability closure of CHECK 7).
> The paper's claim for all 0<a<b rests on its Theorem 1 hand proof; this program certifies the computational
> content, not the induction beyond those ranges.
> gaps this program does NOT close, stated rather than hidden:
> (1) QUOTATION FIDELITY.  That Problem 7.1, Proposition 7.2 and Observation 5.2 are numbered and worded in
> arXiv:2607.27989v1 as the paper quotes them -- in particular that Prop 7.2's printed comparison range
> really ends at b-2, that Q really is {(m+1)q-1 : m>=1}, and that Problem 7.1 as printed carries no
> restriction on the heap size n -- is input data here, not a verified fact; the identifier itself is not
> resolved by this program either.  Every verdict of CHECK 1d, CHECK 5b and CHECK 7 is conditional on
> those three transcriptions being right; a single typo in any of them would invert the negative verdict.
> (2) UNRESTRICTED (a,b).  Nothing above is an induction; the finite sweeps cannot certify all 0<a<b.
> (3) DECORATION, declared.  These checks state true but definitional facts and could not fail for any input:
> CHECK 3a-iii (q|b <=> q|(a-1) is forced by b=q+a-1); CHECK 3a-iv (implied by CHECK 3a-ii on the same
> parameters); CHECK 3b-i's empty dominated sets (a=1 gives at most one option per player, so emptiness is
> structural); CHECK 3b-ii (q divides 0 always); CHECK 6b (B_0 ends at q-1 and tau=q-1 by definition);
> the 'spread = a-1' clause of CHECK 4a.  They are kept for the record, but the load-bearing checks are
> 1a-1c, 2a, 2b, 3a-i/ii, 3c-ii/iii, 4a's residue-collision clause, 4b, 5a, 5b, 6a and 7a-7c.
> (4) CHECK 4b is a consequence of CHECK 2a, not independent evidence: once the derived order is known to be
> 'same block and r<=s', the blockwise-extremum domination rule follows algebraically.
