# The Ivan–Wang Conjecture on Induced N-Saturation for 1 n 5

`the-ivan-wang-conjecture-on-induced-n-saturation-for-1-n-5`

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
passes. The recorded run reports **105 checks, all passing**:

    VERDICT: ALL 105 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    728948996e3eacbcb3b9bf2886c2ad393fda62be91d55bb17d93eadf389f62a7

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> SCOPE: the census of Table 1 is reproduced in full for n=3,4,5 -- every one of the 57 / 6,476 / 8,656,937 candidates was tested, nothing was sampled or truncated.
> SCOPE: for n=1,2,3,4 sat*(n,N) is also obtained with no lemma at all, by scanning every subfamily of B_n. For n=5 that scan (2^32 subfamilies) is out of reach, so the n=5 minimum rests on the reduction of Lemma 3, whose computational content (no copy of N meets {} or [n]) is checked at n=5 and whose full set equality is checked at n=3,4.
> SCOPE: at n=6 only the 4-subset facts are checked -- the Lemma 4 edge count, the Lemma 2 P4 criterion, role uniqueness and the no-extremes fact. No census, no minimality and no orbit claim is verified at n=6, and nothing whatever for n >= 7. The paper claims nothing for n >= 6 either.
> NOT RE-RUN: the paper's own first implementation, which this bundle does not ship; this file is an independent reimplementation written from the paper's specification alone. Also not re-run: the n=5 brute-force scan of all 2^32 subfamilies of B_5, so sat*(5,N) here rests on the reduction of Lemma 3 as set out in the SCOPE line above; and no census, minimality or orbit computation at n >= 6.
