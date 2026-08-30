# A Completely Uniform Nested Steiner Quadruple System of Order 56

`a-completely-uniform-nested-steiner-quadruple-system-of-order-56`

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

Python 3.9 or later, standard library only: no third-party package, no external data file, no
solver and no network. Every input it consumes is printed in the paper -- the thirteen base
blocks of display (4), the 1-factorization formula (5), the doubling of Lemma 2, and the 29
base-64 lines of Section 4 -- and it reads nothing from disk. All of its arithmetic is exact
integer arithmetic; it contains no floating-point number and therefore makes no float decision.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **31 checks, all passing**, in about one second on one core:

    VERDICT: ALL 31 CHECKS PASS

It rebuilds the SQS(56) from the paper's recipe, decodes the witness, and re-derives every
quantity the paper asserts: the block and pair counts, the exact coverage of all 364, 3276 and
27720 triples at the three stages, the pair degree 27, both SHA-256 digests, the exact
termination of the base-3 decoding, the agreement of the base-64 witness with the redundant
decimal prefix printed beside it over all 200 digits, the twenty nested blocks displayed in the
paper, and the tabulation showing all 1540 pairs present at multiplicity 9. It also runs
controls in both polarities: two forced positives, re-deriving the settled orders v = 8 and
v = 14 of Lu's table by its own exhaustive search rather than from Lu's printed pairings; and
two anti-controls, confirming that its decision procedure rejects the witness with a single
ternary digit altered and rejects the trivial all-zero nesting.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    60a155f81c59b76ea33a8d7f5b9fd45a3fa7dbea4befb26d074ffa12dcd5345c

## Scope

**Which of the two kinds of program this is.** It RE-DERIVES the paper's claim; it does not merely
confirm an object handed to it. The SQS(56) is rebuilt from scratch out of the thirteen base blocks
of display (4), the 1-factorization formula (5) and the doubling of Lemma 2 -- the blocks are not
stored anywhere in the program or read from any file -- and only the ternary nesting vector is taken as
given data, decoded from the base-64 string printed in the paper. Every number Theorem 1 asserts is then
recomputed, including the exact triple coverage at all three orders. This was tested adversarially:
altering one base-64 character in a scratch copy makes the program exit 1 with 5 failed checks, and
altering that character *and* patching both printed SHA-256 constants to the corrupted values still
leaves 3 failed checks, among them `nesting-is-uniform-of-multiplicity-9` (multiplicity set becomes
3..17 instead of the constant 9). The program is not checking digests only.

The program's own statement of what it does not cover, quoted from its output:

> NOT RE-RUN: nonexistence of anything. The two remaining open cells of the v=56 row -- 924
> nested pairs at multiplicity 15, and 1260 at multiplicity 11 -- are only checked for
> arithmetic consistency above; no search for or against them was run, and they stay open.
>
> NOT RE-RUN: whether this SQS(56) is rotational, i.e. admits an automorphism of order 55
> fixing a point. No automorphism computation was performed on the 6930 blocks; the paper
> claims only that rotationality is not USED.
>
> NOT RE-RUN: uniqueness, isomorphism, or any enumeration of the nestings of this or any other
> SQS(56). One witness is exhibited; 0% of the space is exhausted, and no claim of minimality
> or canonicity is made.
>
> NOT RE-RUN: Table 4 of Lu at v = 20, 26, 32, 38, 44, 50. Only v = 8 and v = 14 are
> re-derived here as forced positives; the other six rows are taken on the published authority
> and bear on novelty, never on correctness.
>
> NOT RE-RUN: the search that produced the witness. This program re-checks the witness only;
> the randomised descent that found it is not reproduced, and none of the above depends on it.

Two further limits are not the program's to check and are stated in the paper instead. First,
v = 56 is a gap and not the frontier: a completely uniform nested SQS(62) is already published,
so this note settles a hole below a settled larger order. Second, the published journal text of
Lu, and the subscription reviewing databases, could not be read; a prior observation of this
cell recorded only there would be invisible to us.
