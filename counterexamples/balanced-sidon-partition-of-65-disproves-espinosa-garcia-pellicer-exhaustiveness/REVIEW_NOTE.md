# A Counterexample for [65] to an Exhaustiveness Conjecture of Espinosa-García and Pellicer

`balanced-sidon-partition-of-65-disproves-espinosa-garcia-pellicer-exhaustiveness`

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
passes. The recorded run reports **93 checks, all passing**:

    VERDICT: ALL 93 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    4b3e3b2cfbf3aba126bf5b4220a2405793a31a8c408c3b35f29715a334f9ceee

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> gaps: steps between the checked facts and the paper's claim
> The checks above establish, exactly and without tolerance, that W
> and rho(W) are two distinct balanced Sidon-Ramsey partitions of
> [65] with profile (10,10,9,9,9,9,9) lying outside the ten-member
> family F.  Three steps from there to the paper's sentence are NOT
> covered by any check above and are not counted as passing:
> (i)  READING OF THE QUOTED SENTENCE.  'We conjecture that these
> and their reflected partitions ... are all balanced Sidon-
> Ramsey partitions with those parameters' is read as an
> exhaustiveness claim (these are all of them).  Under the
> weak reading (each of these is one) the sentence is true and
> the refutation does not apply.  No program can settle this.
> Two textual facts favour the strong reading and a reader
> should weigh them: the preceding sentence already gives the
> parameters, so the weak reading is redundant; and the
> paragraph opens 'So far we have found $5$ ...', i.e. it is
> about how many such partitions exist.
> (ii) TRANSCRIPTION OF THE FIVE DISPLAYED PARTITIONS.  F was
> built from partitions transcribed by hand from
> arXiv:2309.08553v1, Further.tex lines 34-77.  With no
> source file supplied, NOTHING here rules out a sixth
> displayed partition, or a dropped part, which would make
> F too small and could admit a genuine family member as a
> counterexample.  Run
> python3 verify.py /path/to/Further.tex
> to close this gap; section (8) then checks it in both
> directions, not merely that each transcribed part occurs.
> (iii) VERSION OF RECORD.  Whether the quoted sentence survives
> into Discrete Appl. Math. 378 (2026), 120-124,
> doi:10.1016/j.dam.2025.07.002 cannot be decided here: the
> article is paywalled.  The paper scopes its claim to
> arXiv:2309.08553v1, which is what is checked.  A human with
> journal access should read around p.123.
