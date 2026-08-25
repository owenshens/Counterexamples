# An Infinite Family of Counterexamples to an Atomic-Weight Conjecture of Dargad and Larsson

`truncated-support-family-disproves-dargad-larsson-atomic-weight-conjecture`

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
passes. The recorded run reports **33 checks, all passing**:

    VERDICT: ALL 33 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    82674291f6c25f8cdbcdc21b52a4b7b82e00973c120337819d546f78c7936060

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> --- GAPS between the checks above and the paper's claims ---
> G1  SCOPE, and the one that matters most.  Every check is conditioned on P3, this program's transcription of Conjecture 7.1 and of its standing hypotheses -- above all that a = 1 is admissible.  Dargad-Larsson's preprint also carries a '3 <= a < b' hypothesis elsewhere, in an unrelated section on domination.  If the conjecture's real scope excludes a = 1, or if its hypotheses are not 0<a<=b, a<b, tau<b-a, n>=R_betahat, then EVERY CHECK HERE STILL PASSES and the paper is nevertheless wrong.  No offline computation can close this; it needs a fidelity audit against the preprint itself, which this program does not attempt because it would need a network fetch of a third-party document.
> G1a Two concrete ways G1 could bite, both invisible to every check here: (i) if DL's equations (5)/(13) give an R_betahat larger than 2t+4, then n = 2t+4 is below the smallest governed heap and the instances are out of scope; (ii) if DL's Definition 5.1 gives Right the subtraction set {tau,...,b} or {tau+1,...,b-1} rather than {tau+1,...,b}, the games computed here are the wrong games.  The rule anchors below test the built trees against the prose of THIS paper's own proof, which establishes internal consistency, not fidelity to Dargad-Larsson's Definition 5.1.
> G2  FINITE vs UNIVERSAL.  Theorem 1 is stated for every integer t >= 1.  It is verified here for t = 1..12 only, by direct computation; no induction on t is verified.  What IS fully settled is the Corollary, because falsifying a conjecture needs only ONE instance: any single t above suffices, so the refutation of Conjecture 7.1 does not depend on the universal quantifier.
> G3  ATOMIC WEIGHT IS NEVER COMPUTED.  This program implements the order relation, not the atomic-weight calculus, so it can REFUTE a predicted atomic weight but can never independently CONFIRM one.  The imported facts are aw(up)=1, aw(*)=0, additivity, and the characterisation P6(C): for an all-small G, aw(G)=0 exactly when G is infinitesimal with respect to up, i.e. when for EVERY integer j >= 1 the game j*G is less than or confused with up and greater than or confused with down.  The refutation used here is the contrapositive: if some j >= 1 gives up <= j*(G - P*up) or j*(G - P*up) <= down, then aw(G) != P.  Note that the STRICT form 'aw(G)=0 => down < G < up' is NOT the criterion and is false -- aw(*)=0 while * is confused with up and with down, as the engine self-test above records -- so it is not used or relied on anywhere in this program.  aw(down) = -1 and aw(down+down+*) = -2 are then derived, each from an identity the order engine decides here.
> G4  BOTH ENGINES SHARE ts().  The order test and the outcome induction are independent decision procedures, but they read the same game trees, so their agreement cannot detect a mis-transcribed subtraction set.  The rule-anchor checks above address that separately, by testing the built trees against move-by-move facts the paper's own proof states in prose.
> G5  THE REMARK IS THE WEAKEST CHECK HERE.  The Remark's appeal to DL Table 4 is answered by a direct computation of TS_2(n;3,7) for 12 <= n <= 34, but that computation can only fail to refute: it does not read Table 4, it does not confirm the tabulated weights, and at the largest heap sizes the j-sweep was cut to j = 1 by the node cap (the check's own detail line prints how many).  A refutation there would have contradicted the paper's Remark and its abstract, so the check is not empty -- but its PASS establishes consistency with the Remark, not the Remark.  The abstract's clause that tau = (b-a)/2 is not itself the obstruction therefore still rests on Dargad-Larsson's Table 4, which nothing in this program has read.
> --- NOT COVERED (declared, not checked) ---
> Fidelity to arXiv:2607.27989v1 itself -- that Conjecture 7.1, its standing hypotheses, Definition 5.1 and Table 4 are as transcribed in P1/P3 above -- is NOT checked: it would need a network fetch of a third-party preprint, and this program is offline and stdlib-only.  It is NOT counted as a check.
> The two atomic-weight facts P6(A)/(B) and the characterisation P6(C) are imported from the cited literature, not proved here; every check that uses them says so on its own line.
> The values of DL Table 4 are never read.  The paper's Remark, and with it the abstract's clause that tau = (b-a)/2 is not itself the obstruction, is therefore only tested to the extent that the order test found no refutation at the boundary instance (3,7,2) for 12 <= n <= 34 -- and at the largest of those heap sizes the sweep was cut to j = 1 by the node cap, as the check's detail line records.
> Theorem 1 is verified for t = 1..12 only, by direct computation over finite game trees; no induction on t is verified anywhere, and no second implementation of these values is shipped alongside this one.
> Total distinct game positions constructed: 83029 ; order-test memo entries: 66245
> Checks whose label says "strengthening (beyond the paper's claim)" test additional instances, not statements of the paper; the paper's own claims are the Theorem 1 and refutation checks.
