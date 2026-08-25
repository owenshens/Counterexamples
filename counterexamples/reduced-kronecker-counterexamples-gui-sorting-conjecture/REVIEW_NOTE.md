# Counterexamples to Gui's Sorting Conjecture for Reduced Kronecker Coefficients

`reduced-kronecker-counterexamples-gui-sorting-conjecture`

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
passes. The recorded run reports **18 checks, all passing**:

    VERDICT: ALL 18 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    02d5d37e5cfd5c8442af0c4f09868e0fa33e5ad3cdd534a765a8f63dca2faa9a

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOT RE-RUN HERE: (i) the reduced Kronecker coefficients are limits, and stability is OBSERVED over finitely many N rather than proved -- the two load-bearing values are accepted for k=2..12 once 3 consecutive admissible N agree, and the longer 5-term constant window is exhibited only for k=2..7, so no window longer than 3 terms is shown for k=8..12; (ii) the character-based confirmation of the family covers k=2..12, and k=13..400 rests on the closed formula for nu=(1), which is itself validated against characters on all pairs of partitions of size <= 5; (iii) the exhaustive census over pairs runs to total size 6, beyond the total size 4 needed for the minimality claim, but its nu range is TRUNCATED at |nu| <= |lambda|+|mu|: that truncation is supported only by the evidence that every coefficient with |nu| in [total+1,total+2] vanishes, and |nu| > total+2 is untested; (iv) no external catalogue or table of Kronecker coefficients is consulted; (v) the multipartition corollary is NOT verified here beyond its two-factor content: it rests on identifying the two partitions written (lambda cup mu)^[1,2] and (lambda cup mu)^[2,2] in the cited multipartition conjecture with (sort1,sort2), which is a reading of an external paper's notation and is TAKEN ON TRUST -- no check above tests it, because coding the r-block index-residue split and comparing its r=2 output with sort1,sort2 would only compare one transcription of a single definition with another and could not fail; (vi) nothing here consults the cited paper itself (no network access and no external file is read), so the numbering and wording of the external statements quoted as Conjecture 5.5, Theorem 5.6 and Conjecture 5.7 of its published version are TAKEN ON TRUST, as is the attribution of the pair (1,)/(1, 1, 1) to that Theorem 5.6 -- the inequality for that pair is nevertheless computed here, only its attribution is not.
