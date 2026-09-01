#!/usr/bin/env python3
"""Verification program for

    "A Hash-and-Separable Matrix of Order 9"

(the paper's earlier title carried "D'yachkov's Open Problem 1 Answers Yes"; the paper as revised
claims only that a C_HS(q,4)-matrix with q < 13 exists, which is what Open problem 1 asks for, and
claims no priority for it.  Nothing printed below should be read as settling that problem.)

Python 3.9+, STANDARD LIBRARY ONLY, exact integer arithmetic, no floating point anywhere.

It reads the objects PRINTED IN THE PAPER -- the 9x9 matrix M of Section 2 and the symbol table of
Section 2 -- and re-derives every quantity the paper claims.  Nothing is read from disk.

The published control matrices are NOT printed in the paper (Examples 7, 8, 9 and 10 of the source,
including the C_H(8,4) and the C_HS(13,4) of Example 10, of which the paper quotes only the two
values).  They are transcribed here from the source's e-print, with the source's own labels, purely
as CONTROLS: they are what makes a negative answer from these checkers believable.  Their line
numbers in the source's e-print `lect_dse.tex` are given beside each one.
"""

import itertools
import sys

STAR = 0                                   # internal code for the source's `*`

# ----------------------------------------------------------------------------------------------
# check harness
# ----------------------------------------------------------------------------------------------
_n = 0
_bad = []


def PASS(name, detail=''):
    global _n
    _n += 1
    print('PASS %s%s' % (name, (' [%s]' % detail) if detail else ''))


def CHECK(name, cond, detail=''):
    if cond:
        PASS(name, detail)
    else:
        _bad.append(name)
        print('FAIL %s [%s]' % (name, detail))


def NOTE(s):
    print('NOTE %s' % s)


def HEAD(s):
    print('')
    print('=== %s' % s)


# ----------------------------------------------------------------------------------------------
# the object, exactly as printed in Section 2 of the paper
# ----------------------------------------------------------------------------------------------
PAPER_MATRIX = """
1 2 3 4 * * * * *
5 * * * 9 * 4 * 8
6 * * * * 8 3 2 *
7 * 9 * * 4 * * 2
* 8 6 * 1 * * 4 *
* * * 9 * 1 * 7 3
* 9 * * * * 1 5 6
* 7 * 6 3 5 * * *
* * 5 8 2 * 7 * *
"""

# the table of Section 2: symbol -> (rows carrying it, columns carrying it), 1-indexed
PAPER_SYMBOL_TABLE = {
    1: ((1, 5, 6, 7), (1, 5, 6, 7)),
    2: ((1, 3, 4, 9), (2, 5, 8, 9)),
    3: ((1, 3, 6, 8), (3, 5, 7, 9)),
    4: ((1, 2, 4, 5), (4, 6, 7, 8)),
    5: ((2, 7, 8, 9), (1, 3, 6, 8)),
    6: ((3, 5, 7, 8), (1, 3, 4, 9)),
    7: ((4, 6, 8, 9), (1, 2, 7, 8)),
    8: ((2, 3, 5, 9), (2, 4, 6, 9)),
    9: ((2, 4, 6, 7), (2, 3, 4, 5)),
}

# the negative control at the witness's own k = 4: D'yachkov's C_H(8,4), Example 10,
# `lect_dse.tex` lines 3764--3772.  It is hash but NOT separable.  The paper quotes Example 10's
# two values but does not print this matrix, so it is transcribed here from the source.
PAPER_CH_8_4 = """
* * * 1 * 2 3 5
* * 1 * 2 * 4 6
* 1 * * 3 4 * 7
1 * * * 5 6 7 *
* 2 3 4 * * * 8
2 * 5 6 * * 8 *
3 5 * 7 * 8 * *
4 6 7 * 8 * * *
"""


def parse(block):
    """A printed matrix -> list of rows of ints, with 0 for `*`. Rejects anything else."""
    rows = []
    for line in block.strip().split('\n'):
        toks = line.split()
        if not toks:
            continue
        rows.append([STAR if t == '*' else int(t) for t in toks])
    return rows


def show(C):
    return ' / '.join(' '.join('*' if v == STAR else str(v) for v in r) for r in C)


# ----------------------------------------------------------------------------------------------
# the three propositions of the source, implemented literally
# ----------------------------------------------------------------------------------------------
def viol_P5(C, q, k):
    """Proposition 5 (`lect_dse.tex` 3616--3629): every symbol of [q] occurs exactly k times;
    every row and every column holds exactly q-k stars; no symbol repeats in a row or column."""
    bad = []
    cnt = {}
    for i in range(q):
        for j in range(q):
            v = C[i][j]
            if v != STAR:
                cnt[v] = cnt.get(v, 0) + 1
    for a in range(1, q + 1):
        if cnt.get(a, 0) != k:
            bad.append(('symbol-count', a, cnt.get(a, 0)))
    for i in range(q):
        if sum(1 for j in range(q) if C[i][j] == STAR) != q - k:
            bad.append(('row-stars', i + 1))
    for j in range(q):
        if sum(1 for i in range(q) if C[i][j] == STAR) != q - k:
            bad.append(('col-stars', j + 1))
    for i in range(q):
        vs = [C[i][j] for j in range(q) if C[i][j] != STAR]
        if len(set(vs)) != len(vs):
            bad.append(('row-repeat', i + 1))
    for j in range(q):
        vs = [C[i][j] for i in range(q) if C[i][j] != STAR]
        if len(set(vs)) != len(vs):
            bad.append(('col-repeat', j + 1))
    return bad


def viol_P6(C, q):
    """Proposition 6 (3643--3654), the hash clause, by exhaustive search over all (i,j,m,n):
    if c_i(j) = c_m(n) = a != * with i != m and j != n then c_i(n) = c_m(j) = *."""
    bad = []
    for i in range(q):
        for j in range(q):
            a = C[i][j]
            if a == STAR:
                continue
            for m in range(q):
                if m == i:
                    continue
                for n in range(q):
                    if n == j:
                        continue
                    if C[m][n] == a and (C[i][n] != STAR or C[m][j] != STAR):
                        bad.append(('P6', a, i + 1, j + 1, m + 1, n + 1))
    return bad


# The six forbidden 3x3 forms of Proposition 7 (3660--3685).  The source records them as "the
# permutations of the same three columns"; they are generated here from one base form so that
# claim is exercised rather than assumed.
_BASE = [[None, 'a', 'c'],
         ['a', None, 'b'],
         ['c', 'b', None]]
PATTERNS = [[[_BASE[r][p[c]] for c in range(3)] for r in range(3)]
            for p in itertools.permutations(range(3))]


def viol_P7_printed(C, q):
    """Proposition 7 as printed: no 3x3 submatrix on rows i1<i2<i3 and columns j1<j2<j3 equals any
    of the six forms, under a symbol assignment with a, b, c pairwise distinct."""
    bad = []
    for tri_i in itertools.combinations(range(q), 3):
        for tri_j in itertools.combinations(range(q), 3):
            sub = [[C[i][j] for j in tri_j] for i in tri_i]
            for P in PATTERNS:
                asg = {}
                ok = True
                for r in range(3):
                    for c in range(3):
                        want, got = P[r][c], sub[r][c]
                        if want is None:
                            if got != STAR:
                                ok = False
                        elif got == STAR:
                            ok = False
                        elif want in asg:
                            if asg[want] != got:
                                ok = False
                        else:
                            asg[want] = got
                        if not ok:
                            break
                    if not ok:
                        break
                if ok and len(set(asg.values())) == 3:
                    bad.append((tuple(i + 1 for i in tri_i), tuple(j + 1 for j in tri_j),
                                tuple(asg[x] for x in 'abc')))
    return bad


def viol_P7_ordered(C, q):
    """The same condition written out over ORDERED row and column triples against the single base
    form -- a second encoding of Proposition 7, with no pattern list at all.  Per Section 4 of the
    paper it is an encoding inside this one program and shares this program's one transcription of
    M, so its agreement with the printed form bears on the encodings and is not an independent
    confirmation of the matrix."""
    bad = []
    for i1, i2, i3 in itertools.permutations(range(q), 3):
        for j1, j2, j3 in itertools.permutations(range(q), 3):
            if C[i1][j1] != STAR or C[i2][j2] != STAR or C[i3][j3] != STAR:
                continue
            a, b, c = C[i1][j2], C[i2][j3], C[i1][j3]
            if STAR in (a, b, c) or len({a, b, c}) != 3:
                continue
            if C[i2][j1] == a and C[i3][j2] == b and C[i3][j1] == c:
                bad.append((tuple(x + 1 for x in (i1, i2, i3)),
                            tuple(x + 1 for x in (j1, j2, j3)), (a, b, c)))
    return bad


# ----------------------------------------------------------------------------------------------
# the raw definitions of the source (Definitions 2, 4 and 5), a second and independent route
# ----------------------------------------------------------------------------------------------
def induced_code(C, q):
    """Definition 2's induced code: one codeword (a, i, j) per non-star cell, 1-indexed."""
    return [(C[i][j], i + 1, j + 1) for i in range(q) for j in range(q) if C[i][j] != STAR]


def viol_homogeneous(B, q, k):
    bad = []
    for c in range(3):
        cnt = {}
        for w in B:
            cnt[w[c]] = cnt.get(w[c], 0) + 1
        for a in range(1, q + 1):
            if cnt.get(a, 0) != k:
                bad.append(('coord%d' % (c + 1), a, cnt.get(a, 0)))
    return bad


def viol_distance(B):
    """Two codewords of length 3 may agree in at most one coordinate (Hamming distance >= 2)."""
    return [(x, y) for x, y in itertools.combinations(B, 2)
            if sum(1 for c in range(3) if x[c] == y[c]) >= 2]


def viol_hash3(B):
    """Definition 5 with s = 3: every 3-subset has a coordinate on which all three differ."""
    return [e for e in itertools.combinations(B, 3)
            if not any(len({w[c] for w in e}) == 3 for c in range(3))]


def viol_sep3(B):
    """Definition 4 with s = 3: the coordinate-wise unions of DISTINCT subsets of size <= 3 are
    themselves distinct."""
    seen, bad = {}, []
    for n in (1, 2, 3):
        for e in itertools.combinations(B, n):
            key = tuple(frozenset(w[c] for w in e) for c in range(3))
            if key in seen:
                bad.append((seen[key], e))
            else:
                seen[key] = e
    return bad


def n_sep3_subsets(m):
    return m + m * (m - 1) // 2 + m * (m - 1) * (m - 2) // 6


# ==============================================================================================
# Step 1 -- the object as printed
# ==============================================================================================
HEAD("Step 1: the object printed in Section 2 of the paper")

M = parse(PAPER_MATRIX)
q, k = 9, 4
CHECK('matrix_is_9_by_9', len(M) == 9 and all(len(r) == 9 for r in M),
      'shape %dx%s' % (len(M), sorted({len(r) for r in M})))
CHECK('every_entry_is_a_star_or_a_symbol_in_1_to_9',
      all(v == STAR or 1 <= v <= q for r in M for v in r), 'alphabet {*} u [9]')
nonstar = [(i, j) for i in range(q) for j in range(q) if M[i][j] != STAR]
CHECK('nonstar_cells_number_36_equals_k_times_q', len(nonstar) == k * q,
      '%d = %d*%d' % (len(nonstar), k, q))
NOTE('M = ' + show(M))

# ==============================================================================================
# Step 2 -- Proposition 5
# ==============================================================================================
HEAD("Step 2: Proposition 5 (source lines 3616--3629)")

cnt = {a: sum(1 for i in range(q) for j in range(q) if M[i][j] == a) for a in range(1, q + 1)}
CHECK('each_of_the_9_symbols_occurs_exactly_4_times',
      all(cnt[a] == k for a in cnt), 'counts %s' % sorted(cnt.values()))
rs = [sum(1 for j in range(q) if M[i][j] == STAR) for i in range(q)]
cs = [sum(1 for i in range(q) if M[i][j] == STAR) for j in range(q)]
CHECK('every_row_holds_exactly_q_minus_k_equals_5_stars', all(x == q - k for x in rs), 'rows %s' % rs)
CHECK('every_column_holds_exactly_q_minus_k_equals_5_stars', all(x == q - k for x in cs), 'cols %s' % cs)
CHECK('no_symbol_repeats_within_a_row',
      all(len({M[i][j] for j in range(q) if M[i][j] != STAR}) == k for i in range(q)))
CHECK('no_symbol_repeats_within_a_column',
      all(len({M[i][j] for i in range(q) if M[i][j] != STAR}) == k for j in range(q)))
v5 = viol_P5(M, q, k)
CHECK('proposition_5_holds_zero_violations', not v5, '%d violations' % len(v5))

# ==============================================================================================
# Step 3 -- Proposition 6, twice: literally, and by the paper's block criterion
# ==============================================================================================
HEAD("Step 3: Proposition 6 (3643--3654), literally and by the block criterion")

v6 = viol_P6(M, q)
CHECK('proposition_6_holds_by_exhaustive_search_over_all_i_j_m_n', not v6,
      '%d violations; %d (i,j,m,n) quadruples examined' % (len(v6), k * q * (q - 1) * (q - 1)))

rows_of = {a: sorted(i + 1 for i in range(q) for j in range(q) if M[i][j] == a) for a in range(1, q + 1)}
cols_of = {a: sorted(j + 1 for i in range(q) for j in range(q) if M[i][j] == a) for a in range(1, q + 1)}
CHECK('each_symbol_occupies_4_distinct_rows_and_4_distinct_columns',
      all(len(rows_of[a]) == k and len(cols_of[a]) == k for a in rows_of))
CHECK('the_symbol_table_of_section_2_agrees_with_the_matrix',
      all(tuple(rows_of[a]) == PAPER_SYMBOL_TABLE[a][0] and tuple(cols_of[a]) == PAPER_SYMBOL_TABLE[a][1]
          for a in range(1, q + 1)),
      'all 9 rows of the printed table re-derived from M')

blocks_ok, block_sizes = True, []
for a in range(1, q + 1):
    R = [i - 1 for i in rows_of[a]]
    Cc = [j - 1 for j in cols_of[a]]
    inside = [(i, j) for i in R for j in Cc if M[i][j] != STAR]
    block_sizes.append(len(inside))
    if len(inside) != k or any(M[i][j] != a for i, j in inside):
        blocks_ok = False
CHECK('the_nine_4x4_blocks_R_a_x_C_a_contain_only_the_four_cells_of_a', blocks_ok,
      'non-star cells per block %s, all equal to a' % sorted(set(block_sizes)))
NOTE('so each block is a 4x4 permutation matrix over {a, *}, which is the paper\'s '
     'one-line equivalent of Proposition 6')

# the criterion is only worth anything if it can also say NO; test it against a matrix that
# satisfies Proposition 5 and violates Proposition 6.
CTRL_C_6_4 = parse("""
1 2 4 5 * *
* 1 2 4 5 *
* * 1 3 4 6
6 * * 1 3 4
5 6 * * 2 3
3 5 6 * * 2
""")                                        # the plain C(6,4) of Example 7 of the source


def block_criterion_fails(C, qq, kk):
    for a in range(1, qq + 1):
        R = [i for i in range(qq) for j in range(qq) if C[i][j] == a]
        Cc = [j for i in range(qq) for j in range(qq) if C[i][j] == a]
        ins = [(i, j) for i in set(R) for j in set(Cc) if C[i][j] != STAR]
        if len(ins) != kk or any(C[i][j] != a for i, j in ins):
            return True
    return False


CHECK('block_criterion_and_proposition_6_agree_on_the_control_C_6_4',
      block_criterion_fails(CTRL_C_6_4, 6, 4) and len(viol_P6(CTRL_C_6_4, 6)) == 68,
      'control C(6,4) of Example 7: both routes reject; %d P6 violations'
      % len(viol_P6(CTRL_C_6_4, 6)))
CHECK('control_C_6_4_satisfies_proposition_5_so_the_rejection_is_p6_alone',
      not viol_P5(CTRL_C_6_4, 6, 4), '0 P5 violations')

# ==============================================================================================
# Step 4 -- Proposition 7, in two encodings, plus the paper's reduction
# ==============================================================================================
HEAD("Step 4: Proposition 7 (3660--3685), two encodings and the reduction of Section 2")

CHECK('the_six_forbidden_forms_are_the_six_column_permutations_of_one',
      len(PATTERNS) == 6 and len({tuple(map(tuple, P)) for P in PATTERNS}) == 6,
      '6 distinct forms generated from the base form')

p7p = viol_P7_printed(M, q)
n_tri = len(list(itertools.combinations(range(q), 3)))
CHECK('proposition_7_holds_printed_form_no_forbidden_3x3_submatrix', not p7p,
      '%d hits over %d x %d increasing triples x 6 forms = %d pattern tests'
      % (len(p7p), n_tri, n_tri, n_tri * n_tri * 6))
p7o = viol_P7_ordered(M, q)
CHECK('proposition_7_holds_ordered_form_no_forbidden_configuration', not p7o,
      '%d hits over %d x %d ordered triples' % (len(p7o), q * (q - 1) * (q - 2), q * (q - 1) * (q - 2)))
CHECK('the_two_proposition_7_encodings_agree_on_M', (not p7p) == (not p7o), 'both report zero')

# The paper's reduction: a violation needs three rows pairwise sharing symbols, three pairwise
# distinct symbols a, b, c placed as in the base form, and three distinct columns.  Given
# Proposition 6 the three star cells on the diagonal are IMPLIED, so only the column equations
# remain.
def share(C, qq, i, m):
    A = {C[i][j] for j in range(qq) if C[i][j] != STAR}
    B = {C[m][j] for j in range(qq) if C[m][j] != STAR}
    return sorted(A & B)


def reduction(C, qq):
    """-> (candidates, violations). A candidate is a row triple plus a distinct symbol triple in
    the right sharing pattern; a violation is a candidate meeting all three column equations."""
    colof = {(C[i][j], i): j for i in range(qq) for j in range(qq) if C[i][j] != STAR}
    cand, bad = 0, []
    for i1, i2, i3 in itertools.combinations(range(qq), 3):
        for a in share(C, qq, i1, i2):
            for b in share(C, qq, i2, i3):
                for c in share(C, qq, i1, i3):
                    if len({a, b, c}) != 3:
                        continue
                    cand += 1
                    j2, j1, j3 = colof[(a, i1)], colof[(a, i2)], colof[(b, i2)]
                    if (colof[(b, i3)] == j2 and colof[(c, i3)] == j1
                            and colof[(c, i1)] == j3 and len({j1, j2, j3}) == 3):
                        bad.append((i1 + 1, i2 + 1, i3 + 1, a, b, c))
    return cand, bad


pairs_share = {(i, m): share(M, q, i, m) for i, m in itertools.combinations(range(q), 2)}
CHECK('all_84_row_triples_of_M_pairwise_share_a_symbol_so_no_triple_is_free',
      all(len(v) >= 1 for v in pairs_share.values()) and n_tri == 84,
      '%d of %d row pairs share at least one symbol' % (sum(1 for v in pairs_share.values() if v),
                                                        len(pairs_share)))
cand, red_bad = reduction(M, q)
CHECK('the_reduction_leaves_exactly_162_candidate_configurations', cand == 162,
      '%d candidates' % cand)
CHECK('every_candidate_fails_one_of_the_three_column_equations', not red_bad,
      '%d of %d survive' % (len(red_bad), cand))
CHECK('the_reduction_and_the_two_brute_force_encodings_agree_on_M',
      (not red_bad) and (not p7p) and (not p7o),
      'three encodings inside one program, all zero; per Section 4 of the paper they share one '
      'transcription of M, so their agreement bears on the encodings and is NOT an independent '
      'confirmation')

# ==============================================================================================
# Step 5 -- the raw definitions, not the propositions
# ==============================================================================================
HEAD("Step 5: Definitions 2, 4 and 5 evaluated directly on the induced code")

B = induced_code(M, q)
CHECK('the_induced_code_has_36_codewords_of_length_3', len(B) == 36 and all(len(w) == 3 for w in B),
      '|B| = %d' % len(B))
CHECK('the_induced_code_is_9_4_3_homogeneous', not viol_homogeneous(B, q, k),
      'every symbol appears exactly 4 times in each of the 3 coordinates')
CHECK('any_two_codewords_agree_in_at_most_one_coordinate', not viol_distance(B),
      '0 pairs at Hamming distance < 2 over %d pairs' % (len(B) * (len(B) - 1) // 2))
h3 = viol_hash3(B)
CHECK('the_induced_code_is_3_hash_by_definition_5', not h3,
      '%d bad triples of %d' % (len(h3), len(list(itertools.combinations(range(len(B)), 3)))))
s3 = viol_sep3(B)
CHECK('the_induced_code_is_3_separable_by_definition_4', not s3,
      '%d collisions over all %d subsets of size <= 3' % (len(s3), n_sep3_subsets(len(B))))
CHECK('the_definitional_route_and_the_propositional_route_agree_on_M',
      (not h3) and (not s3) and (not v5) and (not v6) and (not p7p),
      'C_HS(9,4) by Propositions 5/6/7 AND C_HS by Definitions 2/4/5')

# ==============================================================================================
# Step 6 -- controls: both polarities, on the source's own labelled matrices
# ==============================================================================================
HEAD("Step 6: controls on the source's own labelled matrices (Examples 8, 9, 10)")

# Published matrices, transcribed from the source's e-print `lect_dse.tex`.
CTRL = {
    'C_S(3,2)  Example 8  (separable, not hash)': (3, 2, parse("1 2 *\n* 1 3\n3 * 2")),
    'C_H(3,2)  Example 8  (hash, not separable)': (3, 2, parse("* 3 1\n1 2 *\n3 * 2")),
    'C_HS(4,2) Example 8': (4, 2, parse("1 * 3 *\n* 1 * 3\n4 * 2 *\n* 4 * 2")),
    'C_H(6,3)  Example 9  (hash, not separable)':
        (6, 3, parse("* * 1 2 3 *\n* 1 * 5 * 3\n1 * * * 5 4\n* 2 5 * * 6\n2 * 4 * 6 *\n3 4 * 6 * *")),
    'C_HS(7,3) Example 9':
        (7, 3, parse("* * 1 2 3 * *\n* 1 * 5 7 * *\n1 * * * * 7 3\n* 2 * * * 5 4\n"
                     "2 * 7 * * * 6\n* 3 * 4 * 6 *\n5 * 4 * 6 * *")),
    'C_H(8,4)  Example 10 (hash, not separable)': (8, 4, parse(PAPER_CH_8_4)),
    'C_HS(13,4) Example 10 (the published record)': (13, 4, parse("""
* * * * * * 4 2 3 * 8 * *
* * * * * * * 7 13 5 1 * *
* * * * * * * 9 10 * * 1 8
* * * * * * 6 12 11 * * 5 *
* * * * * * * * * 10 11 13 3
* * * * * * * * * 9 12 7 2
1 * * 13 * * * * * 6 * * 4
2 7 9 12 * * * * * * * * *
3 6 11 10 * * * * * * * * *
* 5 * * 11 9 13 * * * * * *
8 4 * * 10 12 * * * * * * *
* * 4 5 6 7 * * * * * * *
* * 8 * 3 2 1 * * * * * *
""")),
}

res = {}
for name, (qq, kk, C) in CTRL.items():
    B2 = induced_code(C, qq)
    res[name] = {'P5': len(viol_P5(C, qq, kk)), 'P6': len(viol_P6(C, qq)),
                 'P7p': len(viol_P7_printed(C, qq)), 'P7o': len(viol_P7_ordered(C, qq)),
                 'hash3': len(viol_hash3(B2)), 'sep3': len(viol_sep3(B2)),
                 'homog': len(viol_homogeneous(B2, qq, kk)), 'dist': len(viol_distance(B2)),
                 'nB': len(B2)}
    NOTE('%-46s |B|=%-3d P5=%-2d P6=%-2d P7printed=%-3d P7ordered=%-3d 3hash=%-2d 3sep=%-3d'
         % (name, len(B2), res[name]['P5'], res[name]['P6'], res[name]['P7p'],
            res[name]['P7o'], res[name]['hash3'], res[name]['sep3']))

r13 = res['C_HS(13,4) Example 10 (the published record)']
CHECK('control_the_published_C_HS_13_4_passes_all_three_propositions',
      r13['P5'] == 0 and r13['P6'] == 0 and r13['P7p'] == 0 and r13['P7o'] == 0,
      'P5/P6/P7 all clean, |B| = %d' % r13['nB'])
CHECK('control_the_published_C_HS_13_4_is_C_HS_by_the_raw_definitions_too',
      r13['hash3'] == 0 and r13['sep3'] == 0 and r13['homog'] == 0 and r13['dist'] == 0,
      'so no corpus error hides in the C_HS(13,4) the source prints at q = 13; the paper claims no '
      'priority over it')

r84 = res['C_H(8,4)  Example 10 (hash, not separable)']
CHECK('control_the_published_C_H_8_4_passes_propositions_5_and_6',
      r84['P5'] == 0 and r84['P6'] == 0, 'as the source labels it: a hash code')
CHECK('control_the_published_C_H_8_4_fails_proposition_7',
      r84['P7p'] == 32 and r84['P7o'] == 192,
      '%d printed-form hits and %d ordered-form hits' % (r84['P7p'], r84['P7o']))
CHECK('control_the_published_C_H_8_4_fails_3_separability_by_definition_4',
      r84['sep3'] == 32 and r84['hash3'] == 0,
      '%d union collisions, 3-hash clean -- the definitional route agrees with Proposition 7'
      % r84['sep3'])
NOTE('this is the anti-vacuity control that matters: it is at the SAME k = 4 as the witness, it '
     'comes from the source itself, and it shows the Proposition 7 clauses do real work rather '
     'than being satisfied by everything that passes Propositions 5 and 6')

r63 = res['C_H(6,3)  Example 9  (hash, not separable)']
CHECK('control_the_published_C_H_6_3_fails_proposition_7',
      r63['P5'] == 0 and r63['P6'] == 0 and r63['P7p'] == 4 and r63['P7o'] == 24,
      '%d printed-form and %d ordered-form hits' % (r63['P7p'], r63['P7o']))
r32h = res['C_H(3,2)  Example 8  (hash, not separable)']
CHECK('control_the_published_C_H_3_2_fails_proposition_7_at_one_configuration',
      r32h['P5'] == 0 and r32h['P6'] == 0 and r32h['P7p'] == 1 and r32h['P7o'] == 6,
      'the single unordered hit the source itself marks')
hit = viol_P7_ordered(CTRL['C_H(3,2)  Example 8  (hash, not separable)'][2], 3)
CHECK('control_C_H_3_2_the_hit_is_at_rows_1_2_3_and_columns_1_3_2',
      ((1, 2, 3), (1, 3, 2)) in {(h[0], h[1]) for h in hit},
      'located configuration %s' % (sorted({(h[0], h[1]) for h in hit})[0],))
r32s = res['C_S(3,2)  Example 8  (separable, not hash)']
CHECK('control_the_published_C_S_3_2_fails_proposition_6_but_is_3_separable',
      r32s['P5'] == 0 and r32s['P6'] == 4 and r32s['sep3'] == 0,
      '%d P6 violations, 0 union collisions' % r32s['P6'])
CHECK('control_the_published_C_S_3_2_fails_3_hash_by_definition_5',
      r32s['hash3'] == 2, '%d non-hashing triples' % r32s['hash3'])
r42 = res['C_HS(4,2) Example 8']
r73 = res['C_HS(7,3) Example 9']
CHECK('control_the_published_C_HS_4_2_and_C_HS_7_3_pass_every_check',
      all(x == 0 for x in (r42['P5'], r42['P6'], r42['P7p'], r42['P7o'], r42['hash3'], r42['sep3'],
                           r73['P5'], r73['P6'], r73['P7p'], r73['P7o'], r73['hash3'], r73['sep3'])),
      'both polarities of every checker are now exercised on published objects')
CHECK('the_reduction_also_reproduces_the_negative_controls',
      reduction(CTRL['C_H(6,3)  Example 9  (hash, not separable)'][2], 6)[1]
      and reduction(parse(PAPER_CH_8_4), 8)[1],
      'C_H(6,3): %d surviving; C_H(8,4): %d surviving -- the reduction can say NO'
      % (len(reduction(CTRL['C_H(6,3)  Example 9  (hash, not separable)'][2], 6)[1]),
         len(reduction(parse(PAPER_CH_8_4), 8)[1])))

# ==============================================================================================
# Step 7 -- the shared-incidence counts of Section 2, a counting bound the paper does not state,
# the hand lemma of Section 3, and the conclusion
# ==============================================================================================
HEAD("Step 7: the shared-incidence counts of Section 2, a counting bound the paper does not state, "
     "the hand lemma of Section 3, and the conclusion")

corow = {}
for a, b in itertools.combinations(range(1, q + 1), 2):
    corow[(a, b)] = len(set(rows_of[a]) & set(rows_of[b]))
two = sum(1 for v in corow.values() if v == 2)
one = sum(1 for v in corow.values() if v == 1)
CHECK('of_the_36_symbol_pairs_18_co_occur_in_two_rows_and_18_in_exactly_one',
      two == 18 and one == 18 and two + one == 36, '%d twice, %d once' % (two, one))
CHECK('the_total_symbol_pair_incidence_is_54_equals_q_times_c_k_2',
      sum(corow.values()) == q * (k * (k - 1) // 2) == 54,
      '%d = 9 * C(4,2)' % sum(corow.values()))
CHECK('so_the_symbol_pairs_of_distinct_rows_are_not_distinct_and_the_counting_bound_lapses',
      two > 0, '%d symbol pairs are reused, so 6q <= C(q,2) is not forced' % two)
lin4 = [qq for qq in range(2, 40) if (k * (k - 1) // 2) * qq <= qq * (qq - 1) // 2]
lin3 = [qq for qq in range(2, 40) if 3 * qq <= qq * (qq - 1) // 2]
CHECK('the_counting_bound_under_that_hypothesis_would_force_q_at_least_13_for_k_4',
      min(lin4) == 13, 'least q with C(4,2)q <= C(q,2) is %d' % min(lin4))
CHECK('and_q_at_least_7_for_k_3_matching_the_published_q_3_hs_equals_7',
      min(lin3) == 7, 'least q with C(3,2)q <= C(q,2) is %d' % min(lin3))
CHECK('13_and_7_are_k_squared_minus_k_plus_1_for_k_4_and_k_3',
      4 * 4 - 4 + 1 == 13 and 3 * 3 - 3 + 1 == 7, '|PG(2,3)| = 13, |PG(2,2)| = 7')
NOTE('SCOPE: of this step, only the 18/18 and 54 = q*C(k,2) shared-incidence counts and the '
     'q >= 2k-1 hand lemma correspond to anything in the paper (the (P7) proof of Section 2 and '
     'the q >= 2k-1 lemma of Section 3).  The counting-bound checks -- that reusing symbol pairs '
     'makes C(k,2)q <= C(q,2) lapse, that the bound would otherwise force q >= 13 at k = 4 and '
     'q >= 7 at '
     'k = 3, and the k^2-k+1 / projective-plane identification -- correspond to NOTHING in the '
     'paper as revised, which draws no counting bound and makes no projective-plane remark.  They '
     'are established here, not claimed there.')

CHECK('the_hand_lemma_arithmetic_a_symbol_forces_q_at_least_2k_minus_1',
      all(k - 1 > qq - k for qq in range(k, 2 * k - 1)) and 2 * k - 1 == 7,
      'for q <= 6 a symbol forces k-1 = 3 stars into a row that may hold only q-k <= 2, so q >= 7')
CHECK('M_witnesses_q_4_hs_at_most_9_which_is_strictly_below_the_published_13',
      q == 9 and 9 < 13 and not (v5 or v6 or p7p or p7o or h3 or s3),
      'q_4^HS <= 9 < 13')
CHECK('a_c_hs_q_4_matrix_with_q_below_13_exists_which_is_what_open_problem_1_asks_for',
      q < 13 and k == 4 and not (v5 or v6 or p7p or p7o),
      'a C_HS(q,4)-matrix with q < 13 exists, exhibited at q = 9; the paper does not claim that '
      'Open problem 1 is thereby settled, and claims no priority')

NOTE('NOT RE-RUN, and not claimed by the paper: (a) NON-EXISTENCE of a C_HS(q,4)-matrix at '
     'q = 7 and q = 8, which in the run behind this paper rests on SAT-solver UNSAT and is NOT '
     'reproduced here -- so this program establishes q_4^HS <= 9 and says nothing about whether '
     'q_4^HS equals 9; (b) the cases q = 4, 5, 6, which the hand lemma above settles by '
     'arithmetic and which are NOT searched exhaustively here; (c) the C_HS(10,4), C_HS(11,4) '
     'and C_HS(12,4) matrices that the same search returned, which are not printed in the paper '
     'and are not examined here; (d) the k = 3 satisfying assignment at q = 6 reported by that '
     'search, whose matrix is not exhibited anywhere in this folder, so q_3^HS <= 6 is NOT '
     'verified here; (e) the source\'s six published minima as MINIMA -- only the objects it '
     'prints are re-checked, never the claim that nothing smaller exists; (f) Open problem 2 of '
     'the source (k >= 5), which is untouched; (g) PRIORITY -- whether a C_HS(q,4)-matrix with '
     'q < 13 is already recorded in the literature (Dyachkov-Rykov, Furedi-Ruszinko, the '
     'Kim-Lebedev and Lebedev tables, or the Blackburn k-plex literature) is not determined here, '
     'and the paper claims none, so nothing above should be read as a claim of first construction '
     'or as settling the source\'s Open problem 1.')

print('')
if _bad:
    print('VERDICT: %d CHECK(S) FAILED: %s' % (len(_bad), ', '.join(_bad)))
    sys.exit(1)
print('VERDICT: ALL %d CHECKS PASS' % _n)
sys.exit(0)
