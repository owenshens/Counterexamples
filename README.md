# Counterexamples and Constructions in Combinatorics and Algebra

60 short, self-contained papers. Each refutes a published conjecture or settles an open
existence question by exhibiting an explicit object.

Every paper is **checkable from the document alone**: the refuting or settling object is printed in
full, and any computation is specified in the text in enough detail to be re-run. No paper depends
on files outside this repository.

| directory | n | contents |
|---|---|---|
| [`counterexamples/`](counterexamples/) | 42 | a published conjecture or question is refuted |
| [`constructions/`](constructions/) | 14 | an object whose existence was open is exhibited |
| [`notes/`](notes/) | 4 | expository notes, and re-proofs of results due to others |

**[OPEN_RESULTS.md](OPEN_RESULTS.md)** lists every paper and what it establishes:
**39 refuted conjectures** and **13 settled existence questions**,
reported separately and never added together.

## Building

Each paper is a standalone `amsart` document with no local dependencies:

```sh
tectonic -X compile counterexamples/<name>.tex
```

A compiled PDF is committed beside every source.

## How these were checked

Each paper was read against **the cited source's own text**, not against its own quotation of that
source. For every load-bearing term the source's definition was compared with the paper's, so that a
counterexample cannot pass by relying on a different reading. Each object was then re-derived from the
printed data and re-tested under the source's definitions. Finally an independent reviewer was asked
to break the first reviewer's conclusion; where the two disagreed, the more sceptical verdict was
taken and the paper is not in this collection.

Some results are computer-assisted. Where that is so, the paper specifies the computation rather than
shipping code, and says plainly which parts are checkable by hand.

## Citing

Cite the individual paper by its title.
