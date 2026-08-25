# The Minimum Dimension of a Maximal Face of the 6 6 Completely Positive Cone

`cp6-maximal-face-minimum-dimension-is-seven`

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
passes. The recorded run reports **47 checks, all passing**:

    VERDICT: ALL 47 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    8fdea5ac3d5738e32c3aa356bdb7c5b972c545bf567d1f1c62278a4cd671d5fb

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> --------------------------------------------------------------------------------------
> GAP LEDGER -- premises used by the paper that this program does NOT verify
> --------------------------------------------------------------------------------------
> 1. Lemma 2, the extreme exposer (every maximal proper face of a proper cone K is K cap y^perp for some y generating an extreme ray of K^*): hand proof, no finite computation.
> 2. AHD's classification of the extreme rays of COP^6 into nonexceptional / zero-padded exceptional of order 5 / Cases 1-19 [AHD, Thm 5.1]: cited, not reproved here (38 published pages of case analysis).  If a further class existed, every check above would still pass.
> 3. Hildebrand's criterion that the minimal zeros of an exceptional extreme matrix span R^6 [Hildebrand, Thm 4.5]: cited.  It is the SPANNING hypothesis of Lemma 3; only its necessary consequence (a system of distinct representatives for the supports) is checked above.
> 4. The analytic step that the tabulated component of each piece is the generic one, that {I_alpha} and p are constant on a piece, and that J_alpha can only GROW off the generic locus (HA Sec. 2 / AHD Sec. 5).  This is what extends one tabulated certificate to every matrix of the piece; it is not checkable from support data.
> 5. Equation (2) itself: that Au >= 0 for a zero u of a copositive A, hence no cancellation in u_beta^T A u_alpha.  Two lines by hand; the check above is an identity on stand-ins and no matrix A is ever constructed.
> 6. The O5 branch as the paper argues it (from A e_6 = 0 alone, using no table row). The O5 checks above instead test the tabulated O5 row and its repair, which is a consistency check on different data, neither necessary nor sufficient for that paragraph.  That row is a single hand transcription which has not been collated with the published print, so its failure of the symmetry gate is a fact about the transcription and is NOT offered as an erratum in [HA].
> 7. The ENTIRE upper bound low(6) <= 7: the integer 7 is taken from Holmgren-Zhang and no 7-dimensional maximal face is constructed here (the paper exhibits none either).  If HZ were wrong, the checks above would still pass and only the lower bound low(6) >= 7 would survive.
> 8. Universality of the two finite sweeps: Lemma 3(1) is checked on c over {-1,0,1,2}^6 plus 8 rationals, and dim F_{aa^T} = 15 on all 62 sign patterns of {-1,1}^6 plus 2 further a's -- samples of infinite families, not proofs of them.
> 9. Fidelity of the two tables to the printed source.  PAPER_TABLE and HA_TABLE were both entered BY HAND from the same print of [HA, Table 1] (which the paper's Table 1 reproduces), in different encodings and read by different parsers.  Their agreement, checked above, excludes divergent typing and parser-specific slips; it does NOT exclude a misreading common to both, and no automated extraction from the published table was performed.
> 6 of the 47 checks are labelled STRUCTURAL / IDENTITY / BOOKKEEPING: they hold for any input
> and carry no falsification power; the remaining 41 are the tests of the paper.
