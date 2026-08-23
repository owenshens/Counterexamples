# Counterexamples

Output of an automated mathematics solver framework we developed. The framework mines the
literature for open conjectures and open existence questions, attacks them, and verifies
what survives. This repository holds the results: 88 short papers.

**[OPEN_RESULTS.md](OPEN_RESULTS.md)** — the list: **63 refuted conjectures**,
**13 settled existence questions**, and **4 conjectures or open
problems settled affirmatively**.

| directory | n | contents |
|---|---|---|
| [`counterexamples/`](counterexamples/) | 66 | a published conjecture or question is refuted |
| [`constructions/`](constructions/) | 14 | an object whose existence was open is exhibited |
| [`notes/`](notes/) | 8 | conjectures proved, expository notes, and re-proofs of results due to others |

Every paper is checkable from the document alone: the object is printed in full, and any
computation is specified in the text. Each was checked against the cited source's own
definitions, and an independent reviewer then tried to break the conclusion.

Where a result rests on a computation too long to do by hand, the paper specifies that
computation precisely enough to be repeated, and states what it does and does not
establish. Verification code and run transcripts are held separately and are available
on request.

```sh
tectonic -X compile counterexamples/<name>.tex
```

A compiled PDF sits beside every source. Cite the individual paper by its title.
