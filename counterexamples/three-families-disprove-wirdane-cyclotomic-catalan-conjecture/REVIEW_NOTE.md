# Counterexamples to Wirdane's Cyclotomic Catalan Conjecture

`three-families-disprove-wirdane-cyclotomic-catalan-conjecture`

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
passes. The recorded run reports **52 checks, all passing**:

    VERDICT: ALL 52 CHECKS PASS

It reads the object exhibited in the paper as input, and derives the quantities it compares
against the paper's statements.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program
that produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    dc87789a14e22ca7a474d3ecd4d57d2eaa19332cfadea0c15f70cc7954732604

## Scope

The program's closing statement of what it does not cover, quoted from its output:

> --- GAPS: what passing these checks does NOT establish -----------
> G-A. Conjecture 6.6's left-hand side is C^{(mu)}_{n,n}, not
> C_n(zeta_mu).  The paper bridges them with its eq. (3),
> whose last step is the [FH] interpretation cited from a
> paywalled 1985 paper.  This program verifies that bridge
> only for 0 <= n <= 8 with mu in {2,3,4,5} (check G7), and
> univariately for n <= 11 (check H1 with [Wirdane, Thm 3.3]).
> The paper's headline counterexample (n,mu) = (2,3) lies
> INSIDE that range, so the headline result is fully certified.
> But the 80-failure sweep and the family members with n > 11
> or mu > 5 are certified here only as statements about
> C_n(zeta_mu); their transfer to Conjecture 6.6 rests on a
> citation this program cannot check.
> G-B. The three families are INFINITE.  Only m <= 14, r <= 9 and
> p <= 31 are tested.  The paper's own proofs of the general
> statements -- F(z)+F(-z)=2 at order two, the G(z) identity
> behind eq. (7), and the ideal argument mod (1-zeta_p) --
> are NOT machine-checked; only their consequences are, and
> only in the finite ranges above.
> G-C. Source fidelity is a document fact.  Beyond the statement
> -> number map printed below, the TRANSCRIPTIONS of Def 6.1,
> Def 6.2, Rem 6.7 and of the statement of Conjecture 6.6
> itself are unverifiable without the e-print.  One mitigation
> is recorded: eq. (3) is insensitive to which endpoint of
> an inversion assigns its residue class, since any assignment
> of each inversion to exactly one class gives sum_i inv_i
> inv, so a slip in Def 6.1 of that kind would not affect it.
> Exactly ONE datum in this program cross-validates the
> transcription at all: check H3, where the C_{3,3}(q)
> tabulated in [Wirdane] is reproduced independently, both
> from the triangle recurrence and from the enumeration of
> S'_{3,3}(312).  Every other quoted string -- Conj 6.6,
> Rem 6.7, Def 6.1, Def 6.2, and the attributions of Prop 2.4
> and Rem 6.4 -- is uncorroborated here, and the recipe below
> is the only route to checking it.
> G-D. Some checks assert MORE than the paper does: B1 (a table
> fixed here, not printed by the paper), C1 and C2 (all
> n <= 30, the paper says n = 3r+2), F4 (n <= 29, the paper
> says n <= 11), F5 (all mu >= 3 pairs, the paper says (2,3))
> and E7.  A FAIL in those is a defect in a target fixed here
> or in an unproved extrapolation, not necessarily in the paper.
> G-E. Four checks are structural and carry no evidence about the
> paper: A2 and A3 (eq. (1) reindexed), G2 (definitional
> emptiness), G5 (a partition of the same loop).  G7 is forced
> by G5 together with G6.  Discount them: of the 52 checks,
> 4 are implementation checks.
> --- source fidelity (a document fact, NOT one of the checks) ------
> The statement -> number map the paper relies on cannot be
> decided by computation; a reader confirms it against
> the numbered statements of arXiv:2605.14682v1, either by
> reading them off the abs page / PDF at
> https://arxiv.org/abs/2605.14682v1
> or from the LaTeX source, whose main file name is not known
> to this program and must be located rather than assumed:
> curl -sL https://arxiv.org/e-print/2605.14682v1 -o w.tar.gz
> mkdir w && tar xzf w.tar.gz -C w
> grep -l documentclass w/*.tex        # the main file, MAIN
> <any LaTeX engine> MAIN.tex          # then read off the
> numbering below
> the cyclotomic Catalan conjecture            expected: Conjecture 6.6
> the definition of inv_i                      expected: Definition 6.1
> the definition of C_{n,k}(q_1,...,q_mu)      expected: Definition 6.2
> C_{n,k}(q) = sum over S'_{n,k}(312) of q^inv expected: Theorem 3.3
> the [FH] interpretation of C_n(q)            expected: Remark 3.4
> C_{n,n}(q) = C_n(q), its proof deferred there expected: Proposition 2.4
> specialization of all q_i to a common q      expected: Corollary 6.5
> the mu = 2 sign rule contradicted here       expected: Remark 6.7
