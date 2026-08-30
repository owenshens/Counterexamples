# Depth-Four Hole Refill Is Dry on a 2969-Word (23,6,10) Constant-Weight Code

`depth-four-hole-refill-is-dry-on-a-2969-word-23-6-10-constant-weight-code`

Supporting material for this paper: the program that checks its computational claims, and a
record of that program being run. The result is a closure, not a counterexample — an open
problem stated by Lysenstoen (arXiv:2607.19550v1) is answered affirmatively.

## Contents

| file | |
|---|---|
| `paper.tex`, `paper.pdf` | the paper, including the full listing of the code `a23.6.10.2969H` in Appendix A |
| `verify.py` | verification program |
| `verify.output.txt` | recorded run |

## The verification program

```sh
python3 verify.py
```

Python 3.9 or later, standard library only: no third-party package, no network and no external
data file. Its only input is the listing of Appendix A of the paper, reproduced verbatim inside
the program as `C0_LISTING`; the two copies are byte-identical. From those 2969 words it
re-derives every integer the paper prints, including the eight that the paper takes from its
source. It prints one line per check and a closing verdict, and exits 0 only if every check
passes. The recorded run reports **39 checks, all passing**, in 87 seconds on one core:

    VERDICT: ALL 39 CHECKS PASS

The checks fall into six groups. **A** — the listing parses, the 2969 words are distinct, all of
weight 10, and their maximum pairwise intersection is exactly 7, so the code has minimum
distance exactly 6; the SHA-256 identity matches. **B** — the blocker census over all
1,141,097 candidates, built twice by disjoint methods (2-swap neighbourhoods; an 8-subset
containment index sweeping all C(23,10) masks) which agree on all 30,247 candidates with at
most four blockers, and which independently establish that no candidate has an empty blocker
set. **C** — the source paper's own eight published integers about this code, reproduced,
together with the disambiguation of the two different pair-counting conventions it uses. **D** —
the deletion space closes to 70 + 3567 + 143,112 + 4,948,069 = 5,094,818, and is verified to be
the least fixed point (every member is generated, and adjoining any further blocker set to any
member of size at most three produces nothing new). **E** — the sweep: all 5,094,818 deletion
sets, no shard and no cap; the 14,447 whose candidate set is larger than the deletion set each
receive an exact maximum-clique computation, and the remaining 5,080,371 are settled by that
count alone, a clique being unable to exceed the number of candidates. Both histograms and their
cross-agreement are reproduced and no witness exists. **F** — the break-even exchanges
at k = 1, 2, 3, 4, each rechecked as a full O(n^2) validity test on the refilled 2969-word code.

## Provenance

`verify.output.txt` holds the program's output, preceded by a short header and followed by an
exit status, both written by the run harness. The header records the SHA-256 of the program that
produced the output, so the two files can be matched:

```sh
shasum -a 256 verify.py
```

    e09bb7f631b32f18e02d10f14c16dd3b22b276e415bb1dc660fa0ece71bbd600

The run recorded there was made on this control plane under Python 3.9.6, not on a remote
worker: the program needs about 90 seconds of one core and a measured peak of 1.29 GiB, which is
why it can be re-run by a referee unaided.

## Scope

The program **re-derives the paper's claim**; it does not merely confirm the exhibited object.
Group E decides every one of the 5,094,818 deletion sets of the reduced space, and by Lemmas 1
and 2 of the paper that space carries the whole depth-at-most-4 question, so groups A–E together
are a proof of Theorems 1 and 2 from the listing alone. Group F is the separate, and much
smaller, confirmation of the exhibited break-even exchanges. Corrupting the exhibited object
makes the program fail: substituting the conflicting word `4A3705` for `4E3305` in the claimed
maximum clique fails checks F4 and F5, and perturbing one hex digit of an inserted word fails
F6, each with exit status 1.

The program's own closing statements of what it does not cover, quoted from its output:

> NOT RE-RUN: Section XII of Brouwer-Shearer-Sloane-Smith 1990, which the source paper describes
> as certifying k-optimality of this move for k = 2..5, is paywalled and was never read. If it
> certifies k-optimality for the 2969-word code specifically, the |D| = 4 layer of Theorem 1 is
> not new. Nothing in this program can settle that.

> NOT RE-RUN: the |D| <= 3 layer of the sweep (70 + 3567 + 143112 = 146749 sets) RE-PROVES a
> theorem already published in the source paper; only the |D| = 4 layer is new, and within it
> the source's own depth-unbounded theorem already excludes every S all of whose members carry
> at most two blockers.

> NOT RE-RUN: in check E7 the -2, -1 and 0 buckets are UPPER-BOUND-DERIVED for the 5080371 sets
> with |Cand(D)| <= |D|: the recorded margin there is min(|Cand(D)|,|D|) - |D| and no clique
> search was run, so a set counted in one of those buckets may truly sit lower. Only the -3
> bucket and the 14447 live margins are measured. The conclusion "+1 or more is empty" is
> untouched, a clique being unable to exceed |Cand(D)|.

> NOT RE-RUN: depth 5 and above; any code other than C_0; and any lower bound on A(23,6,10).
> C_0 is superseded -- the source paper itself proves A(23,6,10) >= 2979 and Brouwer's table
> reads 2992 -- so nothing here moves a bound.

> NOT RE-RUN: the two forced-positive controls on modified codes (2968 and 2961 words) that were
> used during the search to show the decision procedure can return YES. Their role is replaced
> here by check F6, which exhibits genuine break-even exchanges at every k in 1..4, and by check
> B7, which compares two disjoint censuses.

Two further limits belong to the paper rather than to the program, and §5 of the paper states
both. The first is attribution. Section XII of Brouwer–Shearer–Sloane–Smith 1990 has since been
located in open access (`https://neilsloane.com/doc/Me153.pdf`); §5 of the paper quotes it, and
it contains both the reduction used here and the statement that k-optimality was achieved for
k = 2..5 in the range of that paper's tables. Whether that certification covers this particular
2969-word listing is not determined, so the paper claims no priority for any layer of Theorem 1.
The program's scope note above, written before the text was located, still describes it as
paywalled and unread; the recorded run is left exactly as it was made. The second limit is that
the code `a23.6.10.2969H` has been superseded twice, so the theorem is about an object that is
no longer a record holder and no entry of Brouwer's table changes.

Finally, the attainment claims of §4 of the paper — the deletion set with five insertable
candidates, its ten pairwise intersections, and the break-even exchanges at k = 1, 2, 3, 4 —
need no program at all: every object involved is printed in the paper and the arithmetic is on
23-bit masks. Only Theorem 1 and Theorem 2 require the census and the sweep.
