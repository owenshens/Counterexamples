# A Counterexample to a Local Tor-Length Conjecture of Boix and Lima–Pereira

`complete-intersection-counterexample-boix-lima-pereira-tor-conjecture`

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
passes. The recorded run reports **22 checks, all passing**:

    VERDICT: ALL 22 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    1ff0a0679d737a1d7de53871b0b4fb45f15898b3bc5d2971386f75d4b9a0ac8a

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT VERIFIED HERE, and asserted by the paper: (1) the localisation step S -> S_m -- every length printed above is an exact Q-dimension over the graded ring S, never an ell_R over R = S_m, and the paper's argument that finite length makes localisation preserve lengths and commute with Tor is taken on trust, so the theorem's own statement (which is written entirely over R) is reached only through that untested step; (2) the paper's own proofs of the two facts this program instead re-derives numerically -- the prime-avoidance argument for sqrt(J_S) = m and the appeal to Koszul duality for the symmetry of the Tor vector; (3) the second, intersection route to Tor_1 enumerates graded pieces only in degrees 0..7, so no degree above 7 is inspected by it (its total does equal the Koszul-homology total, which is what closes the gap, but that is an argument, not a computation over the higher degrees); (4) and, the largest gap of all, the CONTENT and the NUMBERING of the statements being refuted. The shape of the conjectured equality, the shape of the total bound, the labels Conjecture 4.2, Question 5.1, Question 5.4 and Proposition 1.9, and the bibliographic identifiers of the cited article (its DOI and arXiv number) are all transcribed from the paper, which transcribes them from a source this program cannot open: nothing here reads the cited article, and this run has no network access by policy. Every arithmetic fact above can therefore be correct while the corollary is attached to the wrong numbered statement; that attribution needs a human with the cited article in hand.
