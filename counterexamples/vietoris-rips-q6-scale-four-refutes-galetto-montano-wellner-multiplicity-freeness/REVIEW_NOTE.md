# A Negative Answer to Problem 8.8 of Galetto, Montaño, and Wellner

`vietoris-rips-q6-scale-four-refutes-galetto-montano-wellner-multiplicity-freeness`

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
passes. The recorded run reports **37 checks, all passing**:

    VERDICT: ALL 37 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    9630ec10fd54ae5cccbc6f39fd1ef0a4afbecffd0a4e71724b824508be60ae87

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> LIMITATIONS  (what this file does NOT establish)
> | 1. Lemma 2, Psi(g) = I(F_g;-1), is ASSUMED here. The paper proves it
> |    in full (page 2, by the Hopf trace formula); this program does not
> |    re-derive that proof. It computes only the right-hand side; it
> |    never builds a chain complex and never computes homology.
> | 2. The individual dimensions 239 (degree 7) and 14 (degree 15) are
> |    NOT derived here. They are cited by the paper from GMW Appendix
> |    C.3. C9 tests only their SUM, 239 + 14 = 253, against the derived
> |    Psi at the identity class.
> | 3. The claim that reduced H_i(X^{6,4};Q) VANISHES outside {7,15} is
> |    likewise external and unverified here. The averaging step
> |    (1/|G|) sum_g Psi(g) = d7 + d15 needs BOTH that vanishing AND the
> |    fact that 7 and 15 are odd (so both signs are +1). A nonzero
> |    EVEN-degree group in GMW's computation would enter with the
> |    opposite sign and break the identification -- this file cannot
> |    close that gap. Doing so would mean building the rational chain
> |    complex of VR(Q_6;4) with faces up to dimension 15 on 64 vertices.
> | 4. Reduced vs unreduced homology: the theorem is stated with
> |    unreduced H_7, H_15 while the computation is with reduced homology.
> |    They agree in positive degrees, and 7, 15 > 0, so the step is
> |    sound -- but it is a step, not an identity of definitions.
> | 5. WHICH of degree 7 or 15 carries the repeated trivial constituent
> |    is NOT determined (C10). Only 'at least one of H_7, H_15 is not
> |    multiplicity-free' follows from d7 + d15 = 3.
> | 6. Convention: the action is taken as (pi x)_i = x_{pi^{-1}(i)}. The
> |    other convention yields a conjugate of g^{-1}; orbits of g and of
> |    g^{-1} coincide and conjugation by a coordinate permutation is a
> |    Hamming isometry, so every quantity above is unchanged.
> | 7. The theorem is stated over EVERY field of characteristic zero. The
> |    base-change step -- H_i(K;k) = H_i(K;Q) tensor k, and invariants
> |    commuting with it via the Reynolds idempotent |G|^{-1} sum_g g --
> |    is proved as prose in the paper (proof of Theorem 1, page 4) and
> |    is NOT modelled here. Everything this file computes is over Z and
> |    concerns Q-coefficients only. The step is standard and correct,
> |    but here it is assumed, not verified.
> | 8. C6's enumerative brute force cannot reach the classes with
> |    |O_g| > 20. Those are covered instead by C6b, a second exact
> |    evaluator using a different branch rule and no component
> |    multiplicativity, under a node budget; any class the budget did
> |    not reach is named in the C6b output rather than silently passing.
> | Note on items 1-7: equation (1) -- items 2 and 3 -- is the ONE input
> | the paper itself declares as taken on trust, in its abstract and in
> | Section 1. Items 1, 4 and 7 are steps the paper PROVES (Lemma 2 on
> | page 2; both others in the proof of Theorem 1 on page 4); they are
> | 'assumed' here only in the sense that this program does not re-verify
> | them. This program's assumption set is therefore strictly LARGER than
> | the paper's, and no item in this block is a gap in the paper.
> VERDICT
