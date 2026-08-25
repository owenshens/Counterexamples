# A seven-vertex counterexample to Conjecture 5.3 of Göbel and Misra

`seven-vertex-colored-paths-disprove-gobel-misra-conjecture`

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
passes. The recorded run reports **25 checks, all passing**:

    VERDICT: ALL 25 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    262d41a5b11c73087fc40dc5ca0f829ec5770124c2b273e410a2b09eaa288b25

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN: the full seven-vertex census over all canonical vertex colorings (11,481 x 1,602 = 18,392,562 paths) is not run here; the seven-vertex work is restricted to a uniform vertex coloring over all 1,602 canonical edge words, and to all 11,481 canonical vertex colorings against the 42 canonical edge words that use at most two colors.
> NOT RE-RUN: the two cited theorems of the source conjecture are not reproduced from their own statements; only the edge condition of the odd-order construction, as quoted in the paper, is implemented.  That condition is necessary for the construction, so refuting it refutes the alternative, but the 'constructions' column of the census table is checked only as the count of pairs left over after the reflection test, and the left-over pairs are then tested against a NECESSARY condition only -- the parity-monochromatic-exchange structure, on the edge words at odd order and on the vertex words at even order.  The table's zero residual is therefore corroborated in every row, never proved.
> NOT RE-RUN: the EVEN rows of the minimality table are not tested against the theorem the paper cites for them.  For even m the paper cites Theorem 3.6 of the source and does not reproduce its statement, so this program has no even-order condition with any textual authority to implement.  The condition applied to the even rows here -- parity-monochromatic exchange on the VERTEX words -- is THE PROGRAM'S OWN INVENTED ANALOGUE of the odd-order edge condition, inferred from the fact that at even order it is the vertex positions that split into two parity classes of equal size.  Nothing in the paper asserts that analogue, and this program does not verify that Theorem 3.6 implies it.  If it is not in fact implied by Theorem 3.6, then the 16-of-16 and 480-of-480 figures corroborate NOTHING about the even rows, and for m = 4 and m = 6 the 'constructions' and 'residual' columns are reproduced here only as their sum.  What those figures do establish, with no theorem involved, is a measured structural fact about the census output alone: every even-order pair the reflection test leaves over has both vertex words alternating in two colors, which is a vanishingly small fraction of the available vertex-word pairs -- the fraction is derived and printed by the selectivity check.  In consequence the EVEN-order half of the paper's minimality Proposition -- that no counterexample has four or six vertices -- is corroborated by this program and is NOT independently checkable from the material shipped beside it: the paper's proof tests each remaining even-order pair against Theorem 3.6 of the source, and nothing here reproduces that theorem's statement.  The ODD rows m = 3, 5 are not subject to this caveat, because there the paper does quote the condition and this program implements what it quotes.
> NOT RE-RUN: at odd order the quoted edge condition forces the second edge word to be the reverse of the first, so under a palindromic vertex word every pair meeting it is already a reflection.  The seven-vertex uniform sub-census leaves 12 pairs over and none of them can meet the condition for that structural reason; that zero is computed and reported as structural, and carries no independent weight.
> NOT RE-RUN: the cited source is never consulted by this program, at any point and by any means.  Everything the paper attributes to it -- the matching expansion, the two construction theorems, and the conjecture being refuted, together with the numbers they are cited under -- reaches this program only as the paper's own quotation of them, and the author names and the preprint identifier are likewise taken from the paper on trust.  This program asserts nothing about any of that: not that the attribution is right, not that it is wrong.  Two things therefore have to be checked at the source and cannot be checked here.  First, that the quoted numbering is the source's numbering; note in particular that the cited work's title concerns colored Gaussian CYCLES while the conjecture refuted here is stated for colored PATHS.  Second, that the reading the paper adopts of the odd-order theorem -- taking the edge-index range in its condition (3) to be {1,...,m-1} -- is the intended one, since it is that reading, exactly as the paper states it, which this program implements and finds violated by the exhibited pair.  The algebra of the counterexample does not depend on any of this and is verified here in full; the identification of WHAT it refutes does depend on it and is not verified here at all.
