# Two non-bold wheel entries of Parida-Moura's Table 4 are equalities: t(W_11) = t(W_12) = 8

`two-non-bold-wheel-entries-of-parida-moura-table-4-are-equalities-t-w11-t-w12-8`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

Two entries of Table 4 of Parida and Moura, *Cover-free families on graphs*, arXiv:2605.12634v1
(12 May 2026), are printed non-bold as `8` under a caption reading "Upper bounds ... where the
bold entries indicate exact values". The paper shows both are equalities:

 t(W_11) = 8 and t(W_12) = 8,

where `W_n` is the wheel on `n` vertices (a hub joined to every vertex of a rim cycle `C_{n-1}`)
and `t(G)` is the least `t` admitting a `G`-cover-free family on `[1,t]` in the sense of that
paper (Definitions 3.1 and 4.1 there). The source poses no question about these two cells: the
absence of boldface records what was not proved there, and the two matching lower bounds are what
the paper supplies. This **confirms** the source's bounds; it refutes nothing.

## What was checked, and how

The two halves of the theorem have different evidential status, and the paper says so in its
Section 4.

* **The upper bound is hand-checkable and needs no code.** The two families that attain it,
 `R_10` and `R_11`, are printed in full on the ground set `{1,...,7}`. Checking them against the
 definition is 221 elementary set differences (20 Sperner conditions and 80 cover-free
 conditions for `R_10`; 22 and 99 for `R_11`), and the hub extension that turns them into wheel
 families at `t = 8` is proved by hand in half a page.
* **The lower bound is not hand-checkable, and no claim here pretends otherwise.** It is two
 exhaustive censuses: no `W_11`-CFF(7,11) and no `W_12`-CFF(7,12) exists, over 14,002,112 and
 14,015,792 search-tree nodes. A referee cannot redo that by eye. What `verify.py` provides is a
 re-run of that search from the objects printed in the paper.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **34 checks, all passing**:

 VERDICT: ALL 34 CHECKS PASS

It reads as input only what the paper prints — the families `R_10` and `R_11`, and the integers
transcribed off Table 4 — and derives everything else: the two graphs, the upper bound via the
hub construction, the two censuses in full (exact node counts and shard counts included), and
the controls the paper reports.

**Runtime.** About 280 s wall clock on 13 cores; the two decisive censuses dominate it and are
sharded over 1024 disjoint prefixes through a `multiprocessing.Pool`. On a single core expect
roughly 35 minutes. Nothing in the program is randomised and no floating-point value is ever
compared or branched on.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 a2fcc202bfcc73aef1eb60ea36ea13c3741bd40d80496208495af80ad33c0028

## Scope

* Only the two equalities `t(W_11) = t(W_12) = 8` are established. The five other non-bold cells
 of Table 4 — `t(C_10)`, `t(P_11)`, `t(C_11)`, `t(P_12)`, `t(C_12)` — are neither proved nor
 used, and the statement "Table 4 is exact throughout for n ≤ 12" is not established here.
* No `t(W_n)` with `n ≥ 13` is claimed.
* **No minimality.** Nothing says `W_11` is the smallest wheel whose table entry was not bold.
* **The transcription is the one unverifiable step.** The verbatim table rows and the two
 integers `8` and `8` are read by hand off the e-print of arXiv:2605.12634v1; no program can
 recompute them.

Note that the closing `NOTE SCOPE` paragraph of `verify.output.txt` records what the program does
not re-run: the five remaining non-bold path/cycle cells of Table 4, any `t(W_n)` with `n >= 13`,
any minimality claim, and the transcription of Table 4 itself. None of these bears on the
theorem, and none is part of the paper's claims.
