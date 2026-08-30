# Counterexamples

Output of an automated mathematics solver framework we developed. The framework mines the
literature for open conjectures and open existence questions, attacks them, and verifies
what survives. This repository holds the results: 182 short papers, each with the code that
checks it.

**[OPEN_RESULTS.md](OPEN_RESULTS.md)** — the list: **112 papers refuting a published
conjecture or question**, **23 exhibiting an object whose existence was open**, and
**47 proving a conjecture or recording a note**. Two statements each carry two independent
papers, which exhibit different objects and are listed separately.

| directory | n | contents |
|---|---|---|
| [`counterexamples/`](counterexamples/) | 112 | a published conjecture or question is refuted |
| [`constructions/`](constructions/) | 23 | an object whose existence was open is exhibited |
| [`notes/`](notes/) | 47 | conjectures proved, expository notes, and re-proofs of results due to others |

## One folder per paper

    <directory>/<paper>/paper.tex            the source
                        paper.pdf            the compiled paper
                        verify.py            a program checking its computational claims
                        verify.output.txt    a recorded run of that program
                        REVIEW_NOTE.md       what the program establishes, and what it does not

```sh
tectonic -X compile counterexamples/<paper>/paper.tex
cd counterexamples/<paper> && python3 verify.py
```

Every paper is checkable from the document alone: the object is printed in full, and any
computation is specified in the text. Each was checked against the cited source's own
definitions, and an independent reviewer then tried to break the conclusion.

The papers added after 25 August 2026 additionally went through two review gates, each run by
two different model families and each followed by a revision pass, and a literature search for
prior art before publication. Two further checks were applied from wave 22 onward: definitional
drift was judged against the cited source's own retrieved text rather than the paper's quotation
of it, and each paper's transcript was checked for *direction* — whether it establishes more than
the paper claims, which is harmless, or less, which is not. A paper was withheld if the statement
it settles turned out not to be the statement its source posed, if its program did not establish
its headline, or if the result was already in the literature. Of the 50 candidates in wave 22, 24
were published: 16 were withheld because the source had never posed the statement, and one because
the result was already in print.

## The verification programs

**147 of the 182 papers ship a program and a transcript.** The remaining 35 need none: their
decisive check is carried out in the paper itself, on the object it prints, and each of those
folders says so in its `REVIEW_NOTE.md`.

Each program is **dependency-free** — Python 3.9 or later, standard library only, no
third-party package and no external data file — so it runs with a bare `python3`. It prints
one line per check and a final verdict, and exits 0 if and only if every check passed. Across
the 151 programs — four folders ship a second, independent program beside the first — the recorded
runs report **8,674 checks, all passing**.

A program takes the object exhibited in its paper as an *input* and derives everything else;
only the derived quantities are checks. Where a program cannot verify a step — typically a
lemma its paper imports by citation — it says so in its own output, and `REVIEW_NOTE.md` quotes that
disclosure rather than paraphrasing it.

The transcripts were produced on machines other than the ones that found the results. Each
opens with the SHA-256 of the program that produced it, so the pairing is checkable:

```sh
shasum -a 256 verify.py    # equals the sha256: line at the top of verify.output.txt
```

## Provenance of this layout

Until 25 August 2026 this repository held one `.tex` and one `.pdf` per paper and no code. That state is
preserved on the branch [`papers-only-2026-08-24`](https://github.com/owenshens/Counterexamples/tree/papers-only-2026-08-24). Nothing in the mathematics changed in the move;
the papers whose sources differ are those where a sentence about the artifact had become
false once the artifact shipped beside it.

## Licence

The verification programs and their recorded output are under the **MIT Licence** ([LICENSE](LICENSE)).
The papers and the per-paper review notes are under **CC BY 4.0** ([LICENSE-PAPERS](LICENSE-PAPERS)).

Cite the individual paper by its title.
