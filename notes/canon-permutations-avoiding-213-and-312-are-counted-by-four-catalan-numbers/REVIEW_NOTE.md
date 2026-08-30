# Canon permutations avoiding 213 and 312 are counted by four Catalan numbers

`canon-permutations-avoiding-213-and-312-are-counted-by-four-catalan-numbers`

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
passes. The recorded run reports **196 checks, all passing**:

    VERDICT: ALL 196 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header
and followed by an exit status, both written by the run harness. The header records the
SHA-256 of the program that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    fb9289378fc4cf0abb3e08d91834f0f8227d480a8d90ceac2d3a6c2e87fe8608

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> Scope
> NOT RE-RUN, stated here rather than left to be discovered.
> NOT RE-RUN (a) RANGE: the exhaustive enumeration reaches k = 6 and the transfer-matrix ladder k = 24 (k = 40 for the witness sigma = 123 alone). Every larger k rests on the hand proof of Theorem 1, which is complete and needs no machine; nothing here is evidence for k > 40 beyond that proof.
> NOT RE-RUN (b) ALPHABET: the alphabet is fixed at three letters except in the two checks that deliberately exhibit the failure of literalisation and of Theorem 2 at n = 4. The paper claims nothing about c_n^k(213,312) for n >= 4 and this program computes nothing about it.
> NOT RE-RUN (c) THE PUBLISHED VALUES are TRANSCRIBED from Elizalde and Luo and are not recomputed from their paper, which this program cannot see. What is checked is that our own definitional census agrees with the transcribed integers -- and one transcribed row (their Fibonacci row) did NOT agree as transcribed; see the disclosed check above.
> NOT RE-RUN (d) Question 6.5 of the source is NOT settled: c_3^k(Lambda) is computed for every pair Lambda only at k = 4, and a closed form valid for all k is proved only on the 3-element orbit of {213,312}, together with the trivial cases Lambda empty and Lambda = S_3.
> NOT RE-RUN (e) NO BIBLIOGRAPHIC CLAIM IS CHECKED HERE. The novelty boundary "new for k >= 3" is a literature judgement, not a computation; the source is a seven-day-old preprint with no citers, so an independent same-week proof is excluded by nothing in this program.
> NOT RE-RUN (f) NO NETWORK, no OEIS lookup, no PDF or e-print parsing: the locator "tex line 946" and any digests quoted in the paper are recorded provenance, not bytes re-verified here.
