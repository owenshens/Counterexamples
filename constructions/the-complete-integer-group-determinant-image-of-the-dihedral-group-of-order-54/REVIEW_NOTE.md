# Pure powers 3^7 and 3^8 in the integer group determinant image of D_54

`pure-powers-3-7-and-3-8-in-the-integer-group-determinant-image-of-d-54`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## What is claimed

Boerkoel and Pinner (*Acta Arith.* **186** (2018), 377--395; arXiv:1802.07336), Theorem 5.4,
determine the integer group determinants of `D_{2p^k}` outside the window `2k+1 <= b <= 3k-1`
in the exponent `b = v_p`. The paper exhibits pure powers filling that window at `p = 3`,
`k = 3`, `G = D_54`, where it is `b in {7,8}`; no priority or novelty claim is made in the
paper, and no literature search beyond the source was carried out there:

1. **Two witnesses (Theorem 1 of the paper's Section 2).** `(x+x^4) + (x+x^4-1)s` and
 `(1+x+x^2+x^4+x^5) + (x+x^2+x^4+x^5)s` in `Z[D_54]` have group determinants exactly
 `3^7 = 2187` and `3^8 = 6561`. The point is not existence but **purity**: both values have
 cofactor `1`, i.e. `v_2 = 0` and nothing outside `3` divides them.
2. **Completeness (Theorem 3), granting [BP].** Granting the necessity and achievability
 statements of Theorem 5.4 of the source at `p = 3`, `k = 3`, which the paper quotes and does
 not reprove,
 `S(D_54) = { 2^a 3^b m : m in Z, gcd(m,6) = 1, a = 0 or a >= 2, b = 0 or b >= 7 }`. It
 determines the single group `D_54`. Two conventions are fixed in Section 1 of the paper so that this is an
 exact set equality rather than a description: `S(G)` is the set of **nonzero** determinants
 (`0 = D_G(0)` is attained for every group and is excluded), and the cofactor `m` ranges over
 **all** of `Z`, not just the positive integers -- `S(D_54)` does contain negative values,
 which the paper's polarity Remark exhibits and the completeness proof now routes through
 explicitly.
The quoted statement of Theorem 5.4 and its page numbers are taken from the version of record,
as cited in Section 1 of the paper.

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package and no external data file.
The program prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **65 checks, all passing**:

 VERDICT: ALL 65 CHECKS PASS

It reads the objects exhibited in the paper -- the four coefficient vectors of Theorem 1 and the
published control witness -- and re-derives every quantity the paper asserts about them. Nothing is read from disk. All arithmetic is exact integer arithmetic
(fraction-free Bareiss determinants, Sylvester resultants, asserted-exact polynomial division);
no decision anywhere depends on a floating-point value.

What it checks, in the order it prints:

1. **The group.** `D_54` is built from the presentation `r^27 = s^2 = 1`, `srs = r^{-1}`;
 order, closure, identity, two-sided inverses, associativity on all `54^3 = 157{,}464`
 triples, the dihedral relation, and the orders of `r` and `s`.
2. **Wiring controls.** A spot value of the `54x54` group determinant that is **negative**
 (`-53`), so the determinant is not being computed up to sign; and multiplicativity
 `D(A)D(B) = D(AB)` on three deterministic pairs of group-ring elements, each with a
 **nonzero** product (up to 54 digits) so that `0 == 0` cannot satisfy the identity.
3. **The witnesses, by four independent routes.** For each of W1 and W2: `h = f f* - g g*`
 recomputed from `f` and `g` and compared with the printed `h`; the shortcut `h = f + f* - 1`
 valid because `g = f - 1`; symmetry of `h` and parity of `h(1)`; the realization stated in
 the paper; then the determinant by (a) the full `54x54` determinant `det(a_{uv^{-1}})` built
 from the multiplication table, which is the **definition**, (b) the `27x27` circulant of `h`,
 (c) `Res(x^27-1, h)`, and (d) `h(1) Res(Phi_3,h) Res(Phi_9,h) Res(Phi_27,h)`. Then that each
 `Res(Phi_{3^m}, h)` is the perfect square `P_m^2 = 9` with `v_3(P_m) = 1`, the valuation
 formula of Lemma 2, `v_3(D)`, `v_2(D) = 0`, cofactor `1`, `D = 3^7` resp. `3^8` exactly, the
 hand identity `P_1 = h(1) - 3H_1` of the paper's Remark, and the sign flip under swapping
 `f` and `g`.
4. **A forced positive control on a published value.** The `p^5` witness printed in the source
 returns the published `243` for `D_18` by all four routes; read instead at `n = 27` it
 returns `6143283 = 3^7 * 53^2`, whose cofactor is *not* `1` -- which is exactly why the pure
 powers of Theorem 1 are load-bearing and the source's own construction does not suffice.
5. **A negative control.** Over a deterministic family of 243 symmetric `h`, no `v_3` lands in
 the band `1..6` that Theorem 5.4's necessity direction forbids, and no `v_2` equals `1`;
 both `v_3 = 7` and `v_3 = 8` do occur in the family, so the test is not vacuous.
6. **The exponent semigroup** of Theorem 3: `{0} u {7,8} u [9,inf)` closed under addition is
 `{0} u [7,inf)`, and without `b = 7, 8` the source alone reaches only `b = 0` or `b >= 9`.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

 3a58b4b6eb0bd51dd28800fe03ab28f33d79016f41b19923a8247ee432130120

The recorded exit status is `0`, on Python 3.9.6; the whole run takes about two seconds on one
core.

**Relation to the original computation.** `verify.py` is not a re-run of the program that first
produced these values. That earlier program was a different file, dispatched to a single-CPU
remote worker, exact integer arithmetic and standard library only, seeded where it used
randomness; its stdout was captured complete and its job exited `0`. It is not part of this
folder, and no claim in the paper depends on it: every input `verify.py` consumes is printed in
the paper, and every value the paper states is re-derived here from those printed inputs.

Two limits of that earlier record should be stated rather than glossed. (i) A second copy of it,
run with several loop bounds cut, **exited nonzero** on an indexing error after the paper's
witness rows and controls had already printed; that job is incomplete, not negative, and nothing
above rests on it. (ii) Two auxiliary scripts from the same record -- a discovery scan and an
exhaustive box census over roughly `6.4` million candidates -- have **no captured output at
all**, so their results are not verifiable from that record. Accordingly the paper makes no
minimality, uniqueness or exhaustiveness claim of any kind, and `verify.py` runs no census.

No other reproduction claim is made here. `python3 verify.py` is the whole of it.

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> NOTE SCOPE: every claim re-derived above is a claim the note makes about the objects printed in it. NOT RE-RUN, and not covered by any check here: (1) the achievability results of [BP] that the completeness theorem consumes -- every 2^a m with (m,6)=1 and a = 0 or a >= 2 (the b = 0 members), and every pure power 3^b with b >= 9 -- which are quoted from the source and not verified; (2) the necessity statement of [BP] Theorem 5.4 / Lemma 4.4 in general, including the p = 2 clauses, which is likewise quoted and not verified: Step 5 is a finite deterministic CONTROL over 243 symmetric h with support in {0,+-1,+-2,+-3,+-4} and coefficients in {-1,0,1}, not a proof; (3) the measure set of any group other than D_54 -- nothing here determines S(D_162), S(D_250) or S of any p > 3; (4) minimality, uniqueness or exhaustiveness of W1 and W2, on which the note makes no claim and this program runs no census.

Beyond that, and beyond Section 5 of the paper, three things a referee should weigh:

- **The completeness statement is only as strong as the source it consumes.** Theorem 3's
 `subseteq` direction is Theorem 5.4 of the source and its `supseteq` direction uses the
 source's achievability for `b = 0` and for `b >= 9`. Neither was re-proved or re-computed
 here; if either is wrong, Theorem 3 falls and Theorem 1 stands.
- **Existence versus purity.** *Mere existence* of measures with `v_3 = 7` and `v_3 = 8` for
 `D_54` follows in a line from constructions already printed in the source, with nontrivial
 cofactors such as `53^2`. Those cofactors cannot be removed, because `S(G)` is closed under
 multiplication and not under division. The new content is that the two powers are attained
 purely; the paper says so in the Remark closing Section 4 and the program exhibits the impure
 value in Step 4.
- **Prior art.** The paper makes no priority or novelty claim and states that no literature
 search beyond the source was carried out. For the record, a search over arXiv, zbMATH,
 Crossref, OpenCitations and Semantic Scholar found no treatment of `D_{2p^k}` with `k >= 3`
 or of order 54, and the source's own most recent successor paper lists `D_{2p^k}`, `k >= 2`,
 nowhere as complete. Two channels were **not** consulted: OpenAlex (repeated HTTP 429 and
 timeouts) and MathSciNet (no access), and no erratum or correspondence channel was checked.
