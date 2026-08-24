# Counterexamples

Output of an automated mathematics solver framework we developed. The framework mines the
literature for open conjectures and open existence questions, attacks them, and verifies
what survives. This repository holds the results: 117 short papers.

**[OPEN_RESULTS.md](OPEN_RESULTS.md)** — the list: **84 papers refuting a published
conjecture or question**, **13 exhibiting an object whose existence was open**, and
**20 proving a conjecture or recording a note**. The 117 papers act on
**115 distinct statements**: two statements each carry two independent papers, which
exhibit different objects and are listed separately.

| directory | n | contents |
|---|---|---|
| [`counterexamples/`](counterexamples/) | 84 | a published conjecture or question is refuted |
| [`constructions/`](constructions/) | 13 | an object whose existence was open is exhibited |
| [`notes/`](notes/) | 20 | conjectures proved, expository notes, and re-proofs of results due to others |

Every paper is checkable from the document alone: the object is printed in full, and any
computation is specified in the text. Each was checked against the cited source's own
definitions, and an independent reviewer then tried to break the conclusion.

```sh
tectonic -X compile counterexamples/<name>.tex
```

A compiled PDF sits beside every source. Cite the individual paper by its title.
