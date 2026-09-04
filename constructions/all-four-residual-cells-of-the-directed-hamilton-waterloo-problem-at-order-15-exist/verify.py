#!/usr/bin/env python3
"""Verification program for the note

    All Four Residual Cells of the Directed Hamilton--Waterloo Problem at Order 15 Exist

It re-derives every quantity that note claims, from the objects PRINTED IN THE NOTE.
The listing in WITNESSES below is character-for-character the listing typeset in
Section 3 of paper.tex, and it is the only input: no external file, no third-party
module, no solver, no floating point.  Python 3.9+, standard library only.

What it checks, for each of the four cells (m,n,r,s) named in the exception clause of
printed Lemma 4.2 of Yetgin-Odabasi-Ozkan:

  * the cell's necessary conditions (m | 15, n | 15, r + s = 14, 15(r+s) = 210 arcs);
  * exactly 14 factors are printed;
  * every printed factor is a set of disjoint cycles covering Z_15 exactly once,
    i.e. a permutation of Z_15;
  * every printed factor is fixed-point-free (no cycle of length 1);
  * every factor's cycle type equals its own printed label m^(15/m);
  * the census of labels is r factors of type m and s of type n;
  * the 210 arcs (x, sigma_i(x)) are EXACTLY the 210 arcs of K*_15, as a set identity
    (none missing, none repeated) -- this is the decisive condition;
  * adjoining the identity as row 0 gives a Latin square of order 15 (every column is a
    permutation of Z_15) -- the equivalent reformulation, computed independently;
  * i |-> sigma_i(0) is a bijection onto {1,...,14} (the "second entry" diagnostic).

and, globally: that the four cells listed are exactly the four cells the quoted
exception clause names, and that the four witnesses are pairwise distinct.

It closes with controls of both polarities, because a checker that accepts everything
proves nothing: one positive control on the machinery itself (the 14 translations of
Z_15 partition the arcs of K*_15, with the predicted cycle types), and five negative
controls in which a witness is deliberately damaged -- a factor reversed, a factor
duplicated, a factor deleted, a cycle type mislabelled, a printed opening parenthesis
dropped -- each of which must be REJECTED by the same checker.

Exit status is 0 if and only if every check passed.
"""

import re
import sys
from collections import Counter
from fractions import Fraction   # exact arithmetic only; no float decision is taken

V = 15
FULL_ARCS = frozenset((a, b) for a in range(V) for b in range(V) if a != b)

# The exception clause of printed Lemma 4.2, as a predicate rather than a list:
# "except possibly for r in {11,12,13} when (m,n)=(3,5) and for r=13 when (m,n)=(3,15)",
# together with r + s = 14.
CLAUSE = [(3, 5, r, 14 - r) for r in (11, 12, 13)] + [(3, 15, 13, 1)]

# ----------------------------------------------------------------------------------
# The objects, exactly as Section 3 of paper.tex prints them.
# ----------------------------------------------------------------------------------
WITNESSES = r"""
HWP*(15; 3^11, 5^3)
F1   3^5   (0,1,2)(3,4,5)(6,7,8)(9,10,11)(12,13,14)
F2   3^5   (0,2,5)(1,6,13)(3,8,11)(4,7,12)(9,14,10)
F3   3^5   (0,3,9)(1,14,4)(2,7,6)(5,13,11)(8,12,10)
F4   3^5   (0,4,12)(1,11,6)(2,8,9)(3,5,14)(7,10,13)
F5   3^5   (0,5,11)(1,3,12)(2,6,4)(7,9,8)(10,14,13)
F6   3^5   (0,6,8)(1,7,5)(2,4,14)(3,13,9)(10,12,11)
F7   3^5   (0,7,14)(1,13,3)(2,11,8)(4,6,10)(5,9,12)
F8   3^5   (0,8,13)(1,4,9)(2,10,5)(3,6,12)(7,11,14)
F9   3^5   (0,9,7)(1,10,2)(3,14,8)(4,11,13)(5,12,6)
F10  3^5   (0,10,6)(1,5,8)(2,9,13)(3,7,4)(11,12,14)
F11  3^5   (0,11,1)(2,13,12)(3,10,7)(4,8,5)(6,14,9)
F12  5^3   (0,12,7,2,3)(1,8,14,5,10)(4,13,6,9,11)
F13  5^3   (0,13,8,4,10)(1,12,9,5,7)(2,14,6,3,11)
F14  5^3   (0,14,1,9,4)(2,12,8,10,3)(5,6,11,7,13)

HWP*(15; 3^12, 5^2)
F1   3^5   (0,1,2)(3,4,5)(6,7,8)(9,10,11)(12,13,14)
F2   3^5   (0,2,7)(1,3,14)(4,12,9)(5,6,10)(8,11,13)
F3   3^5   (0,3,6)(1,9,8)(2,4,13)(5,7,14)(10,12,11)
F4   3^5   (0,4,9)(1,10,2)(3,11,6)(5,8,13)(7,12,14)
F5   3^5   (0,5,10)(1,7,11)(2,12,3)(4,14,8)(6,13,9)
F6   3^5   (0,7,1)(2,8,9)(3,5,13)(4,11,12)(6,14,10)
F7   3^5   (0,8,14)(1,11,4)(2,6,5)(3,9,7)(10,13,12)
F8   3^5   (0,9,13)(1,4,6)(2,11,8)(3,10,14)(5,12,7)
F9   3^5   (0,10,4)(1,8,12)(2,13,6)(3,7,9)(5,14,11)
F10  3^5   (0,11,3)(1,14,13)(2,5,4)(6,9,12)(7,10,8)
F11  3^5   (0,12,8)(1,5,9)(2,10,7)(3,13,4)(6,11,14)
F12  3^5   (0,13,11)(1,12,5)(2,9,14)(3,8,10)(4,7,6)
F13  5^3   (0,6,8,3,12)(1,13,7,4,10)(2,14,9,5,11)
F14  5^3   (0,14,4,8,5)(1,6,12,2,3)(7,13,10,9,11)

HWP*(15; 3^13, 5^1)
F1   3^5   (0,1,2)(3,4,5)(6,7,8)(9,10,11)(12,13,14)
F2   3^5   (0,2,4)(1,9,7)(3,11,14)(5,6,13)(8,10,12)
F3   3^5   (0,4,8)(1,6,2)(3,13,9)(5,7,12)(10,14,11)
F4   3^5   (0,5,10)(1,12,6)(2,9,13)(3,14,4)(7,11,8)
F5   3^5   (0,6,5)(1,14,8)(2,11,3)(4,9,12)(7,13,10)
F6   3^5   (0,7,3)(1,5,12)(2,10,9)(4,13,11)(6,8,14)
F7   3^5   (0,8,9)(1,13,4)(2,5,11)(3,12,10)(6,14,7)
F8   3^5   (0,9,1)(2,14,10)(3,5,8)(4,12,7)(6,11,13)
F9   3^5   (0,10,6)(1,11,5)(2,12,14)(3,8,13)(4,7,9)
F10  3^5   (0,11,7)(1,10,13)(2,8,4)(3,6,12)(5,14,9)
F11  3^5   (0,12,11)(1,7,14)(2,13,8)(3,9,6)(4,10,5)
F12  3^5   (0,13,12)(1,8,11)(2,3,7)(4,6,10)(5,9,14)
F13  3^5   (0,14,13)(1,3,10)(2,7,5)(4,11,6)(8,12,9)
F14  5^3   (0,3,1,4,14)(2,6,9,11,12)(5,13,7,10,8)

HWP*(15; 3^13, 15^1)
F1   3^5   (0,1,2)(3,4,5)(6,7,8)(9,10,11)(12,13,14)
F2   3^5   (0,2,9)(1,7,4)(3,13,10)(5,14,11)(6,8,12)
F3   3^5   (0,4,11)(1,12,3)(2,7,14)(5,8,9)(6,10,13)
F4   3^5   (0,5,13)(1,8,10)(2,3,11)(4,6,14)(7,9,12)
F5   3^5   (0,6,12)(1,11,13)(2,4,10)(3,14,8)(5,9,7)
F6   3^5   (0,7,3)(1,9,14)(2,10,8)(4,12,11)(5,6,13)
F7   3^5   (0,8,1)(2,5,12)(3,10,6)(4,13,9)(7,11,14)
F8   3^5   (0,9,6)(1,5,7)(2,11,3)(4,8,13)(10,12,14)
F9   3^5   (0,10,14)(1,3,12)(2,13,7)(4,9,8)(5,11,6)
F10  3^5   (0,11,7)(1,10,5)(2,6,4)(3,8,14)(9,13,12)
F11  3^5   (0,12,4)(1,14,6)(2,8,5)(3,7,13)(9,11,10)
F12  3^5   (0,13,8)(1,6,11)(2,14,9)(3,5,4)(7,12,10)
F13  3^5   (0,14,5)(1,13,2)(3,6,9)(4,7,10)(8,11,12)
F14  15^1  (0,3,9,1,4,14,13,11,8,7,6,2,12,5,10)
"""

# ----------------------------------------------------------------------------------
# check bookkeeping
# ----------------------------------------------------------------------------------
_n_pass = 0
_n_fail = 0


def check(ok, name, detail=''):
    global _n_pass, _n_fail
    if ok:
        _n_pass += 1
        print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))
    else:
        _n_fail += 1
        print('FAIL %s%s' % (name, (' [%s]' % detail) if detail else ''))
    return ok


def note(text):
    print('NOTE %s' % text)


# ----------------------------------------------------------------------------------
# parsing the printed listing
# ----------------------------------------------------------------------------------
HEAD_RE = re.compile(r'^HWP\*\(15;\s*(\d+)\^(\d+),\s*(\d+)\^(\d+)\)$')
FACT_RE = re.compile(r'^F(\d+)\s+(\d+)\^(\d+)\s+(\(.*\))$')


def parse_listing(text):
    """-> [ {'spec': (m,n,r,s), 'factors': [(index, m_i, k_i, cycletext), ...]} ]"""
    cells = []
    for raw in text.strip().split('\n'):
        line = raw.strip()
        if not line:
            continue
        m = HEAD_RE.match(line)
        if m:
            cells.append({'spec': tuple(int(g) for g in
                                        (m.group(1), m.group(3), m.group(2), m.group(4))),
                          'factors': []})
            continue
        m = FACT_RE.match(line)
        if m:
            if not cells:
                raise SystemExit('factor line before any cell header: %r' % line)
            cells[-1]['factors'].append((int(m.group(1)), int(m.group(2)),
                                         int(m.group(3)), m.group(4)))
            continue
        raise SystemExit('listing line parses as neither a header nor a factor: %r' % line)
    return cells


def cycles_of(cycletext):
    return [[int(t) for t in grp.split(',')]
            for grp in re.findall(r'\(([^)]*)\)', cycletext)]


def tag(spec):
    m, n, r, s = spec
    return '%dp%d_%dp%d' % (m, r, n, s)


def pretty(spec):
    m, n, r, s = spec
    return 'HWP*(15; %d^%d, %d^%d)' % (m, r, n, s)


# ----------------------------------------------------------------------------------
# the checker proper: one function, reused by the target cells and by the controls
# ----------------------------------------------------------------------------------
def evaluate(cell):
    """Every condition of the criterion of the note, on one printed listing.

    Returns a dict of booleans plus detail strings. It never raises on damaged input:
    the negative controls feed it damaged input on purpose.
    """
    m, n, r, s = cell['spec']
    facs = cell['factors']
    out = {}

    out['factor_count'] = (len(facs) == V - 1)
    out['factor_count_detail'] = '%d factors printed, %d required' % (len(facs), V - 1)

    perms_ok, fpf_ok, types_ok = True, True, True
    type_problems, perm_problems = [], []
    arcs = Counter()
    census = Counter()
    perms = []

    for (idx, fm, fk, cycletext) in facs:
        cycs = cycles_of(cycletext)
        pts = [p for c in cycs for p in c]
        if sorted(pts) != list(range(V)):
            perms_ok = False
            perm_problems.append('F%d covers %r' % (idx, sorted(pts)))
        sigma = {}
        for c in cycs:
            if len(c) == 1:
                fpf_ok = False
            for j, x in enumerate(c):
                y = c[(j + 1) % len(c)]
                if x == y:
                    fpf_ok = False
                sigma[x] = y
                arcs[(x, y)] += 1
        if len(sigma) != V:
            perms_ok = False
            perm_problems.append('F%d defines %d of %d images' % (idx, len(sigma), V))
        if fm * fk != V or sorted(len(c) for c in cycs) != [fm] * fk:
            types_ok = False
            type_problems.append('F%d labelled %d^%d but has cycle lengths %r'
                                 % (idx, fm, fk, sorted(len(c) for c in cycs)))
        census[fm] += 1
        perms.append(sigma)

    out['perms'] = perms_ok
    out['perms_detail'] = '; '.join(perm_problems) or 'all factors are permutations of Z_15'
    out['fpf'] = fpf_ok
    out['types'] = types_ok
    out['types_detail'] = '; '.join(type_problems) or 'every cycle type equals its printed label'

    out['census'] = (census == Counter({m: r, n: s}))
    out['census_detail'] = 'labels seen %s, cell wants {%d: %d, %d: %d}' % (
        dict(sorted(census.items())), m, r, n, s)

    missing = sorted(FULL_ARCS - set(arcs))
    repeated = sorted(a for a, c in arcs.items() if c > 1)
    alien = sorted(set(arcs) - FULL_ARCS)
    out['arcs'] = (not missing and not repeated and not alien
                   and sum(arcs.values()) == len(FULL_ARCS))
    out['arcs_detail'] = ('%d arcs, each exactly once, and the set equals A(K*_15)'
                          % sum(arcs.values())) if out['arcs'] else (
        '%d arcs printed; missing %d %s; repeated %d %s; not arcs of K*_15: %d %s'
        % (sum(arcs.values()), len(missing), missing[:4], len(repeated), repeated[:4],
           len(alien), alien[:4]))

    complete = all(len(p) == V for p in perms)
    if complete and perms:
        cols_ok = all(sorted([x] + [p[x] for p in perms]) == list(range(V))
                      for x in range(V))
        firsts = sorted(p[0] for p in perms)
        bij_ok = (firsts == list(range(1, len(perms) + 1)))
    else:
        cols_ok, bij_ok, firsts = False, False, []
    out['latin'] = cols_ok and len(perms) == V - 1
    out['bijection'] = bij_ok and len(perms) == V - 1
    out['bijection_detail'] = 'sigma_i(0) over i = %s' % (firsts if firsts else 'undefined')
    out['perm_maps'] = perms
    return out


def all_ok(ev):
    return all(ev[k] for k in ('factor_count', 'perms', 'fpf', 'types', 'census',
                               'arcs', 'latin', 'bijection'))


def clone(cell):
    return {'spec': cell['spec'], 'factors': list(cell['factors'])}


# ----------------------------------------------------------------------------------
def main():
    print('verification of the note: all four residual cells of the directed')
    print('Hamilton-Waterloo problem at order 15 exist -- HWP*(15; 3^11,5^3),')
    print('HWP*(15; 3^12,5^2), HWP*(15; 3^13,5^1), HWP*(15; 3^13,15^1)')
    print('python %s, exact integer arithmetic only' % sys.version.split()[0])
    print()

    print('=== Step 0: the ground set and the exception clause')
    check(len(FULL_ARCS) == 210 and V * (V - 1) == 210,
          'K_star_15_has_210_arcs', '15*14 = %d = |A(K*_15)|' % len(FULL_ARCS))
    check(Fraction(210, V) == 14,
          'a_solution_needs_exactly_14_spanning_factors', '210/15 = 14 = v-1')

    cells = parse_listing(WITNESSES)
    check(len(cells) == 4, 'listing_parses_as_four_labelled_cells',
          '%d cell headers, %s' % (len(cells), [pretty(c['spec']) for c in cells]))
    check(sorted(c['spec'] for c in cells) == sorted(CLAUSE),
          'the_four_listed_cells_are_exactly_the_cells_of_the_exception_clause',
          'clause r in {11,12,13} at (3,5) and r=13 at (3,15), with r+s=14, gives %s'
          % sorted(CLAUSE))
    note('the clause is re-derived from its own wording, not copied from the listing')

    for cell in cells:
        m, n, r, s = cell['spec']
        ok = (V % m == 0 and V % n == 0 and r + s == V - 1 and V * (r + s) == 210)
        check(ok, 'cell_%s_necessary_conditions' % tag(cell['spec']),
              '%d | 15, %d | 15, r+s = %d+%d = 14, arcs 15*14 = 210' % (m, n, r, s))

    results = {}
    for cell in cells:
        spec = cell['spec']
        t = tag(spec)
        print()
        print('=== %s' % pretty(spec))
        ev = evaluate(cell)
        results[spec] = ev
        check(ev['factor_count'], 'cell_%s_prints_14_factors' % t, ev['factor_count_detail'])
        check(ev['perms'], 'cell_%s_every_factor_is_a_permutation_of_Z_15' % t, ev['perms_detail'])
        check(ev['fpf'], 'cell_%s_every_factor_is_fixed_point_free' % t,
              'no cycle of length 1 and no loop')
        check(ev['types'], 'cell_%s_cycle_types_match_the_printed_labels' % t, ev['types_detail'])
        check(ev['census'], 'cell_%s_census_of_labels_matches_the_cell' % t, ev['census_detail'])
        check(ev['arcs'], 'cell_%s_arc_set_is_exactly_A_of_K_star_15' % t, ev['arcs_detail'])
        check(ev['latin'], 'cell_%s_identity_adjoined_array_is_a_latin_square' % t,
              'all 15 columns are permutations of Z_15')
        check(ev['bijection'], 'cell_%s_sigma_i_of_0_is_a_bijection_onto_1_to_14' % t,
              ev['bijection_detail'])
        check(all_ok(ev), 'cell_%s_satisfies_every_checked_condition_of_the_criterion' % t,
              'factor count, permutation, fixed-point-freeness, cycle types, label '
              'census, arc-set identity, Latin square and second-entry bijection all hold')

    print()
    print('=== Step 2: the four witnesses are distinct objects')
    sets = [frozenset(tuple(sorted(p.items())) for p in results[c['spec']]['perm_maps'])
            for c in cells]
    check(len(set(sets)) == 4, 'the_four_witnesses_are_pairwise_distinct_factor_sets',
          'four distinct sets of 14 permutations')

    print()
    print('=== Step 3: controls, both polarities')

    # positive control on the machinery: the translations of Z_15 partition A(K*_15)
    trans = [{x: (x + c) % V for x in range(V)} for c in range(1, V)]
    tarcs = Counter((x, p[x]) for p in trans for x in range(V))
    check(set(tarcs) == FULL_ARCS and all(v == 1 for v in tarcs.values()),
          'control_positive_the_14_translations_of_Z_15_partition_A_of_K_star_15',
          '%d arcs, each exactly once' % sum(tarcs.values()))

    def cyclens(p):
        seen, lens = set(), []
        for x in range(V):
            if x in seen:
                continue
            L, y = 0, x
            while y not in seen:
                seen.add(y)
                y = p[y]
                L += 1
            lens.append(L)
        return sorted(lens)

    from math import gcd
    predicted = {c: [V // gcd(c, V)] * gcd(c, V) for c in range(1, V)}
    check(all(cyclens(trans[c - 1]) == sorted(predicted[c]) for c in range(1, V)),
          'control_positive_translation_cycle_types_are_as_predicted',
          'x+c has gcd(c,15) cycles of length 15/gcd(c,15): 3^5 for c in {5,10}, '
          '5^3 for c in {3,6,9,12}, 15^1 for gcd(c,15)=1')

    base = cells[0]
    tb = tag(base['spec'])

    # negative control 1: reverse one factor
    dmg = clone(base)
    idx, fm, fk, txt = dmg['factors'][1]
    rev = ''.join('(%s)' % ','.join(str(x) for x in ([c[0]] + c[1:][::-1]))
                  for c in cycles_of(txt))
    dmg['factors'][1] = (idx, fm, fk, rev)
    ev = evaluate(dmg)
    check(not ev['arcs'] and not all_ok(ev),
          'control_negative_reversing_one_factor_of_%s_is_rejected' % tb,
          'arc set identity fails: %s' % ev['arcs_detail'][:90])
    check(not ev['bijection'],
          'control_negative_the_reversal_is_also_convicted_by_the_second_entry_diagnostic',
          ev['bijection_detail'])

    # negative control 2: duplicate a factor
    dmg = clone(base)
    dmg['factors'][1] = dmg['factors'][2]
    ev = evaluate(dmg)
    check(not ev['arcs'] and not all_ok(ev),
          'control_negative_duplicating_a_factor_of_%s_is_rejected' % tb,
          'arc set identity fails: %s' % ev['arcs_detail'][:90])

    # negative control 3: delete a factor
    dmg = clone(base)
    del dmg['factors'][3]
    ev = evaluate(dmg)
    check(not ev['factor_count'] and not ev['arcs'] and not all_ok(ev),
          'control_negative_deleting_a_factor_of_%s_is_rejected' % tb,
          '%s; %s' % (ev['factor_count_detail'], ev['arcs_detail'][:70]))

    # negative control 4: mislabel a cycle type
    dmg = clone(base)
    idx, fm, fk, txt = dmg['factors'][-1]
    dmg['factors'][-1] = (idx, 3, 5, txt)
    ev = evaluate(dmg)
    check(not ev['types'] and not all_ok(ev),
          'control_negative_mislabelling_a_cycle_type_in_%s_is_rejected' % tb,
          ev['types_detail'][:110])

    # negative control 5: drop one printed opening parenthesis
    dmg = clone(base)
    idx, fm, fk, txt = dmg['factors'][0]
    dmg['factors'][0] = (idx, fm, fk, txt.replace('(', '', 1))
    ev = evaluate(dmg)
    check(not ev['perms'] and not all_ok(ev),
          'control_negative_dropping_one_printed_parenthesis_in_%s_is_rejected' % tb,
          ev['perms_detail'][:110])

    print()
    note('SCOPE: this program checks the four exhibited factorizations of K*_15 and '
         'nothing else. It re-derives that each is a solution of its cell, hence that '
         'the four cells of the exception clause of printed Lemma 4.2 of '
         'Yetgin-Odabasi-Ozkan all exist.')
    note('NOT RE-RUN: the searches that FOUND these objects (their logs were not '
         'preserved; a re-run would return different objects, and the objects printed '
         'in the note are the claim); and the cases of printed Lemma 4.2 outside the '
         'four cells, which are the source authors\' own. No non-existence, census or '
         'exhaustion claim is checked or made here: every verdict above is the '
         'existence of an exhibited object.')

    print()
    if _n_fail:
        print('VERDICT: %d of %d CHECKS PASS -- %d FAILED'
              % (_n_pass, _n_pass + _n_fail, _n_fail))
        return 1
    print('VERDICT: ALL %d CHECKS PASS' % _n_pass)
    return 0


if __name__ == '__main__':
    sys.exit(main())
