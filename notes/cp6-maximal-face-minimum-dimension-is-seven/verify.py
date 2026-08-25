#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify.py -- verification artifact for

    "The Minimum Dimension of a Maximal Face of the 6x6 Completely Positive Cone"
    (main theorem: low(6) = 7)

Python 3.9, standard library only.  All arithmetic is exact (int / fractions.Fraction);
there is no floating point anywhere in this file, hence no tolerances.

--------------------------------------------------------------------------------------
TAKEN FROM THE PAPER (data -- transcribed, not derived)
--------------------------------------------------------------------------------------
 1. Table 1 of the paper: 21 rows, each a component label, the integer p, the list of
    minimal-zero supports I_alpha, the list of differences J_alpha \ I_alpha, and the
    certificate the paper prints for that row.  Subsets of {1,...,6} are written in the
    paper as digit strings ("146" = {1,4,6}); that encoding is kept verbatim below.
 2. A SECOND transcription of the source table, Table 1 of
    R. Hildebrand and A. Afonin, "On the structure of the 6x6 copositive cone",
    LAA 693 (2024) 22-38 -- all 22 main components, including the
    O5 row that the paper deliberately does not use.  Recorded in HA's own brace
    notation and parsed by a different parser.  BOTH tables in this file are HAND
    transcriptions entered from the same printed source, so their agreement catches
    divergent typing slips and parser-specific slips only; it does NOT establish
    agreement with the print, because a misreading shared by both transcriptions would
    survive it.  No automated extraction from the published source was performed.
    (HA prints I_alpha once for the split cases 9 and 13; it is duplicated here for
    9.1/9.2 and 13.1/13.2.)
 3. Two integers from the cited literature: the dimension 7 of the maximal face of CP^6
    exhibited by Holmgren-Zhang, and the Kostyukova-Tchemisova range n <= low(n) <= n+3
    for even n >= 6.  The classification of the extreme rays of COP^6 by Afonin,
    Hildebrand and Dickinson is cited, not reproved (38 published pages); this program
    checks what the paper does WITH that classification.

--------------------------------------------------------------------------------------
DERIVED HERE (this is what the checks decide)
--------------------------------------------------------------------------------------
 A. Parsing of both tables; J_alpha = I_alpha u (J_alpha \ I_alpha) reconstructed.
 B. p = |{I_alpha}| for every row, and distinctness of the supports in a row.
 C. Data hygiene: every element lies in {1,...,6}; the printed J_alpha \ I_alpha is
    disjoint from I_alpha (the paper's I_alpha subset J_alpha).
 D. Paper Table 1 vs. the second (brace-notation) transcription of HA Table 1: agreement
    on all 21 shared rows, and that the 22nd HA row is exactly O5.  Both are hand
    transcriptions of the same print (see item 2 above for what that does and does not
    establish).
 E. Every certificate of the paper's table: the pair (beta, alpha) is read out of the
    paper's own printed string, and I_beta subset J_alpha is decided, with both sets
    printed.  Rows with no certificate: p >= n+1 = 7, so Lemma 3(1) applies.
 F. The relation R(beta, alpha) := [I_beta subset J_alpha] must be SYMMETRIC, because
    u_beta^T A u_alpha = u_alpha^T A u_beta for symmetric A.  Full relation matrices
    derived; symmetry decided for all 21 rows.  The O5 row, as transcribed here, is
    derived to FAIL that gate; the paper's own footnote reads that row as understating
    J_alpha for alpha <= 5, and does not use it.  Adding 6 to every J_alpha (justified
    by A e_6 = 0) restores symmetry, after which the row carries a Lemma 3(2)
    certificate.
 G. Robustness sweep: each of HA's 22 components has p >= 7, or at least one unordered
    off-diagonal certificate pair.
 H. The equivalence u_beta^T A u_alpha = 0 <=> I_beta subset J_alpha, checked
    mechanically on nonnegative stand-ins with the prescribed supports and zero sets.
 I. Lemma 3's linear independence, as exact ranks over Q of families of rank-one
    symmetric matrices, including the failure mode when its hypothesis is dropped.
 J. The three nonexceptional face dimensions 15, 20, 15, each as an exact rank of a
    spanning family of rank-one zeros together with a matching upper bound.
 K. The arithmetic bookkeeping: (n^2-5n+8)/2 at n = 6, the interval {6,7} that the
    theorem closes, and the assembled minimum over every branch of the case analysis.

Output: one "PASS "/"FAIL " line per check, then
    VERDICT: ALL n CHECKS PASS   |   VERDICT: m OF n CHECKS FAILED
Exit status 0 if every check passes, 1 otherwise.
"""

import itertools
import sys
from fractions import Fraction

N = 6                      # the paper's ambient dimension
LEMMA3_BOUND = N + 1       # Lemma 3 concludes dim F_A >= n+1
HZ_FACE_DIM = 7            # Holmgren-Zhang: CP^6 has a maximal face of this dimension
KT_RANGE = (N, N + 3)      # Kostyukova-Tchemisova, even n >= 6

# ---------------------------------------------------------------------------------
# DATA 1: Table 1 of the paper, verbatim.
# (label, p as printed, I_alpha string, (J_alpha \ I_alpha) string, certificate string)
# "-" denotes the empty set; "---" in the certificate column means none is printed.
# ---------------------------------------------------------------------------------
PAPER_TABLE = [
    ("1",    6, "12,13,14,25,36,456",              "345,246,23,16,15,-",            "I_2 subset J_1"),
    ("2",    6, "12,13,14,25,356,456",             "345,246,23,16,-,-",             "I_2 subset J_1"),
    ("3",    6, "12,13,14,256,356,456",            "345,24,236,-,-,-",              "I_2 subset J_1"),
    ("4",    6, "12,13,24,345,156,456",            "346,26,15,-,-,-",               "I_2 subset J_1"),
    ("5",    6, "12,13,145,246,346,456",           "35,256,-,-,-,-",                "I_2 subset J_1"),
    ("6",    6, "12,13,245,345,246,356",           "34,256,-,-,-,1",                "I_2 subset J_1"),
    ("7",    6, "15,26,123,234,345,456",           "24,13,6,-,-,-",                 "I_3 subset J_2"),
    ("8",    6, "12,134,135,246,346,256",          "36,5,4,-,-,-",                  "I_3 subset J_2"),
    ("9.1",  6, "12,134,135,246,346,456",          "36,-,-,5,-,2",                  "I_6 subset J_4"),
    ("9.2",  6, "12,134,135,246,346,456",          "356,-,2,-,-,-",                 "I_3 subset J_1"),
    ("10",   6, "12,134,135,246,356,456",          "36,-,-,5,-,2",                  "I_6 subset J_4"),
    ("11",   6, "123,124,125,136,246,346",         "5,5,34,-,-,5",                  "I_3 subset J_2"),
    ("12",   6, "123,124,125,136,246,356",         "-,5,4,-,-,4",                   "I_3 subset J_2"),
    ("13.1", 6, "123,234,345,456,156,126",         "4,15,26,3,-,-",                 "I_2 subset J_1"),
    ("13.2", 6, "123,234,345,456,156,126",         "4,1,6,3,2,5",                   "I_2 subset J_1"),
    ("14",   7, "12,13,14,25,45,36,56",            "345,246,235,146,126,15,234",    "---"),
    ("15",   7, "12,134,135,146,256,356,456",      "34,2,-,-,-,-,-",                "---"),
    ("16",   7, "123,124,125,136,246,346,356",     "-,-,-,-,-,5,4",                 "---"),
    ("17",   7, "123,124,125,136,246,356,456",     "-,-,-,-,-,4,3",                 "---"),
    ("18",   8, "123,234,345,145,125,346,146,126", "-,-,6,6,6,5,5,5",               "---"),
    ("19",   6, "345,145,125,123,156,2346",        "-,6,-,-,4,-",                   "I_5 subset J_2"),
]

# ---------------------------------------------------------------------------------
# DATA 2: Table 1 of Hildebrand-Afonin, LAA 693 (2024) 22-38, all 22 main components,
# in HA's brace notation.  A second HAND transcription of the same print as DATA 1, in a
# different encoding and parsed by a separate parser (see item 2 of the header).
# ---------------------------------------------------------------------------------
HA_TABLE = [
    ("O5",   "{1,2,3},{2,3,4},{3,4,5},{1,4,5},{1,2,5},{6}",
             "{},{},{},{},{},{1,2,3,4,5}"),
    ("1",    "{1,2},{1,3},{1,4},{2,5},{3,6},{4,5,6}",
             "{3,4,5},{2,4,6},{2,3},{1,6},{1,5},{}"),
    ("2",    "{1,2},{1,3},{1,4},{2,5},{3,5,6},{4,5,6}",
             "{3,4,5},{2,4,6},{2,3},{1,6},{},{}"),
    ("3",    "{1,2},{1,3},{1,4},{2,5,6},{3,5,6},{4,5,6}",
             "{3,4,5},{2,4},{2,3,6},{},{},{}"),
    ("4",    "{1,2},{1,3},{2,4},{3,4,5},{1,5,6},{4,5,6}",
             "{3,4,6},{2,6},{1,5},{},{},{}"),
    ("5",    "{1,2},{1,3},{1,4,5},{2,4,6},{3,4,6},{4,5,6}",
             "{3,5},{2,5,6},{},{},{},{}"),
    ("6",    "{1,2},{1,3},{2,4,5},{3,4,5},{2,4,6},{3,5,6}",
             "{3,4},{2,5,6},{},{},{},{1}"),
    ("7",    "{1,5},{2,6},{1,2,3},{2,3,4},{3,4,5},{4,5,6}",
             "{2,4},{1,3},{6},{},{},{}"),
    ("8",    "{1,2},{1,3,4},{1,3,5},{2,4,6},{3,4,6},{2,5,6}",
             "{3,6},{5},{4},{},{},{}"),
    ("9.1",  "{1,2},{1,3,4},{1,3,5},{2,4,6},{3,4,6},{4,5,6}",
             "{3,6},{},{},{5},{},{2}"),
    ("9.2",  "{1,2},{1,3,4},{1,3,5},{2,4,6},{3,4,6},{4,5,6}",
             "{3,5,6},{},{2},{},{},{}"),
    ("10",   "{1,2},{1,3,4},{1,3,5},{2,4,6},{3,5,6},{4,5,6}",
             "{3,6},{},{},{5},{},{2}"),
    ("11",   "{1,2,3},{1,2,4},{1,2,5},{1,3,6},{2,4,6},{3,4,6}",
             "{5},{5},{3,4},{},{},{5}"),
    ("12",   "{1,2,3},{1,2,4},{1,2,5},{1,3,6},{2,4,6},{3,5,6}",
             "{},{5},{4},{},{},{4}"),
    ("13.1", "{1,2,3},{2,3,4},{3,4,5},{4,5,6},{1,5,6},{1,2,6}",
             "{4},{1,5},{2,6},{3},{},{}"),
    ("13.2", "{1,2,3},{2,3,4},{3,4,5},{4,5,6},{1,5,6},{1,2,6}",
             "{4},{1},{6},{3},{2},{5}"),
    ("14",   "{1,2},{1,3},{1,4},{2,5},{4,5},{3,6},{5,6}",
             "{3,4,5},{2,4,6},{2,3,5},{1,4,6},{1,2,6},{1,5},{2,3,4}"),
    ("15",   "{1,2},{1,3,4},{1,3,5},{1,4,6},{2,5,6},{3,5,6},{4,5,6}",
             "{3,4},{2},{},{},{},{},{}"),
    ("16",   "{1,2,3},{1,2,4},{1,2,5},{1,3,6},{2,4,6},{3,4,6},{3,5,6}",
             "{},{},{},{},{},{5},{4}"),
    ("17",   "{1,2,3},{1,2,4},{1,2,5},{1,3,6},{2,4,6},{3,5,6},{4,5,6}",
             "{},{},{},{},{},{4},{3}"),
    ("18",   "{1,2,3},{2,3,4},{3,4,5},{1,4,5},{1,2,5},{3,4,6},{1,4,6},{1,2,6}",
             "{},{},{6},{6},{6},{5},{5},{5}"),
    ("19",   "{3,4,5},{1,4,5},{1,2,5},{1,2,3},{1,5,6},{2,3,4,6}",
             "{},{6},{},{},{4},{}"),
]

# ---------------------------------------------------------------------------------
# check bookkeeping
# ---------------------------------------------------------------------------------
RESULTS = []


def check(label, ok, detail=""):
    """Record one check and print its PASS/FAIL line."""
    RESULTS.append((label, bool(ok)))
    tag = "PASS " if ok else "FAIL "
    line = tag + label
    if detail:
        line += "  [" + detail + "]"
    print(line)
    return bool(ok)


# ---------------------------------------------------------------------------------
# DERIVED A: two independent parsers for the two transcriptions.
# ---------------------------------------------------------------------------------
def parse_paper_sets(field):
    """Paper encoding: comma-separated digit strings, '-' = empty set.
    '146' -> frozenset({1,4,6}).  Returns a list of frozensets."""
    out = []
    for tok in field.split(","):
        tok = tok.strip()
        if tok == "-":
            out.append(frozenset())
            continue
        if not tok or not all(ch.isdigit() for ch in tok):
            raise ValueError("bad paper token %r" % (tok,))
        digits = [int(ch) for ch in tok]
        if len(set(digits)) != len(digits):
            raise ValueError("repeated digit in paper token %r" % (tok,))
        out.append(frozenset(digits))
    return out


def parse_ha_sets(field):
    """HA encoding: '{1,2},{3},{}' -> [frozenset({1,2}), frozenset({3}), frozenset()].
    Deliberately a different parser from parse_paper_sets, so an encoding or parsing slip
    in one table cannot be cancelled by the same slip in the other.  A MISREADING of the
    printed source shared by both hand transcriptions is not caught by this."""
    field = field.strip()
    out = []
    depth = 0
    cur = ""
    for ch in field:
        if ch == "{":
            if depth != 0:
                raise ValueError("nested brace in HA field %r" % (field,))
            depth, cur = 1, ""
        elif ch == "}":
            if depth != 1:
                raise ValueError("unbalanced brace in HA field %r" % (field,))
            depth = 0
            items = [t.strip() for t in cur.split(",") if t.strip() != ""]
            vals = [int(t) for t in items]
            if len(set(vals)) != len(vals):
                raise ValueError("repeated element in HA group {%s}" % (cur,))
            out.append(frozenset(vals))
        elif depth == 1:
            cur += ch
        elif ch not in ", \t":
            raise ValueError("stray character %r in HA field" % (ch,))
    if depth != 0:
        raise ValueError("unterminated brace in HA field %r" % (field,))
    return out


# ---------------------------------------------------------------------------------
# DERIVED A/B: build the rows.  J_alpha = I_alpha u (J_alpha \ I_alpha), and p is the
# cardinality of the derived collection of supports -- never the printed integer.
# ---------------------------------------------------------------------------------
def parse_certificate(text):
    """Read the paper's own certificate string, e.g. 'I_6 subset J_4' -> (6, 4).
    '---' -> None."""
    text = text.strip()
    if text == "---":
        return None
    head, tail = text.split("subset")
    beta = int(head.strip().split("_")[1])
    alpha = int(tail.strip().split("_")[1])
    return (beta, alpha)


def build_rows(table, kind):
    """kind='paper' or kind='ha'.  Returns an ordered list of dicts, each holding the
    parsed I, the parsed J\\I, the DERIVED J, the derived p, and the parsed certificate."""
    parse = parse_paper_sets if kind == "paper" else parse_ha_sets
    rows = []
    for entry in table:
        if kind == "paper":
            label, p_printed, i_field, d_field, cert_field = entry
            cert = parse_certificate(cert_field)
        else:
            label, i_field, d_field = entry
            p_printed, cert = None, None
        isets = parse(i_field)
        dsets = parse(d_field)
        if len(isets) != len(dsets):
            raise ValueError("row %s: %d supports but %d differences"
                             % (label, len(isets), len(dsets)))
        jsets = [isets[k] | dsets[k] for k in range(len(isets))]
        rows.append({
            "label": label,
            "p_printed": p_printed,
            "I": isets,
            "D": dsets,          # J_alpha \ I_alpha exactly as printed
            "J": jsets,          # derived
            "p": len(set(isets)),   # derived: p = |{I_alpha}|
            "n_listed": len(isets),
            "cert": cert,
        })
    return rows


# ---------------------------------------------------------------------------------
# DERIVED B/C: hygiene of the paper's table.
# ---------------------------------------------------------------------------------
def check_data_hygiene(prows):
    ground = frozenset(range(1, N + 1))
    bad_ground, bad_disjoint, bad_p, bad_distinct = [], [], [], []
    print("  derived per-row data (p = |{I_alpha}|, J_alpha = I_alpha u (J_alpha \\ I_alpha)):")
    for r in prows:
        for k in range(r["n_listed"]):
            if not (r["I"][k] <= ground and r["D"][k] <= ground):
                bad_ground.append((r["label"], k + 1))
            if r["I"][k] & r["D"][k]:
                bad_disjoint.append((r["label"], k + 1, sorted(r["I"][k] & r["D"][k])))
        if r["p"] != r["p_printed"]:
            bad_p.append((r["label"], r["p"], r["p_printed"]))
        if r["p"] != r["n_listed"]:
            bad_distinct.append((r["label"], r["n_listed"], r["p"]))
        jtxt = ",".join("".join(str(x) for x in sorted(s)) if s else "-" for s in r["J"])
        print("    %-5s p(derived)=%d p(printed)=%s  J = %s"
              % (r["label"], r["p"], r["p_printed"], jtxt))
    check("paper Table 1: 21 rows transcribed", len(prows) == 21,
          "rows = %d" % len(prows))
    check("paper Table 1: every element of every I_alpha, J_alpha lies in {1,...,6}",
          not bad_ground, "violations = %d" % len(bad_ground))
    check("paper Table 1: printed (J_alpha \\ I_alpha) is disjoint from I_alpha "
          "(consistent with I_alpha subset J_alpha)",
          not bad_disjoint, "violations = %s" % (bad_disjoint if bad_disjoint else "none"))
    check("paper Table 1: minimal-zero supports within a row are pairwise distinct",
          not bad_distinct, "violations = %s" % (bad_distinct if bad_distinct else "none"))
    check("paper Table 1: derived p = |{I_alpha}| equals the printed p in every row",
          not bad_p, "mismatches = %s" % (bad_p if bad_p else "none"))
    # NOTE: J_alpha is BUILT as I_alpha | D_alpha, so "I_alpha subset
    # J_alpha" holds for every conceivable input and cannot falsify anything.  It is
    # recorded as a structural invariant of the reconstruction, not as a test of the
    # paper; the substantive content of that sanity gate is the disjointness check above
    # (the printed difference must not re-list an element of I_alpha).
    subset_ok = all(r["I"][k] <= r["J"][k] for r in prows for k in range(r["n_listed"]))
    sizes_ok = all(len(r["J"][k]) == len(r["I"][k]) + len(r["D"][k])
                   for r in prows for k in range(r["n_listed"]))
    check("STRUCTURAL (holds by construction of J = I u (J\\I), cannot fail): I_alpha "
          "subset J_alpha in all %d rows / %d zeros; the falsifiable content is "
          "|J_alpha| = |I_alpha| + |J_alpha \\ I_alpha|, also checked here"
          % (len(prows), sum(r["n_listed"] for r in prows)),
          subset_ok and sizes_ok, "cardinality identity holds = %s" % sizes_ok)
    return not (bad_ground or bad_disjoint or bad_p or bad_distinct)


# ---------------------------------------------------------------------------------
# DERIVED D: the paper's table against the second hand transcription of HA Table 1.
# Both were entered by hand from the same print, in different encodings and read by
# different parsers, so agreement is evidence against divergent typing and parsing slips
# and NOT evidence of agreement with the published table.
# ---------------------------------------------------------------------------------
def check_against_source(prows, hrows):
    pmap = dict((r["label"], r) for r in prows)
    hmap = dict((r["label"], r) for r in hrows)
    check("HA Table 1: 22 main components transcribed", len(hrows) == 22,
          "rows = %d" % len(hrows))
    only_ha = sorted(set(hmap) - set(pmap))
    only_paper = sorted(set(pmap) - set(hmap))
    check("the one HA component the paper does not read off its table is exactly O5",
          only_ha == ["O5"] and only_paper == [],
          "HA \\ paper = %s ; paper \\ HA = %s" % (only_ha, only_paper))
    check("22 HA components = 21 rows of the paper's table + O5",
          len(hrows) == len(prows) + 1,
          "%d = %d + 1" % (len(hrows), len(prows)))
    labels = sorted(set(pmap) & set(hmap))
    mismatch = []
    for lab in labels:
        p, h = pmap[lab], hmap[lab]
        if p["I"] != h["I"]:
            mismatch.append((lab, "I"))
        if p["D"] != h["D"]:
            mismatch.append((lab, "J\\I"))
        if p["J"] != h["J"]:
            mismatch.append((lab, "J"))
    check("all %d shared rows agree, entry for entry, with the second hand transcription "
          "of HA Table 1 (I_alpha, J_alpha \\ I_alpha and hence J_alpha) -- both tables "
          "were typed from the same print, so this catches divergent typing and parsing "
          "slips, not a misreading common to both" % len(labels),
          not mismatch, "mismatches = %s" % (mismatch if mismatch else "none"))
    # the labels 9 and 13 are the only ones the paper splits; 19 cases -> 21 rows
    base = sorted(set(lab.split(".")[0] for lab in labels), key=lambda s: int(s))
    split = sorted(lab.split(".")[0] for lab in labels if "." in lab)
    check("the paper's 21 rows come from Cases 1--19 with exactly 9 and 13 split in two",
          base == [str(i) for i in range(1, 20)] and sorted(set(split)) == ["13", "9"]
          and len(split) == 4,
          "distinct cases = %d, split labels = %s" % (len(base), sorted(set(split))))
    return not mismatch


# ---------------------------------------------------------------------------------
# DERIVED F: the relation R(beta, alpha) := [I_beta subset J_alpha], which by the
# paper's equation (2) is [u_beta^T A u_alpha = 0].
# ---------------------------------------------------------------------------------
def relation(row, extra=frozenset()):
    """The p x p boolean matrix R[b][a] = (I_{b+1} subset J_{a+1} u extra)."""
    isets, jsets = row["I"], row["J"]
    p = row["n_listed"]
    return [[isets[b] <= (jsets[a] | extra) for a in range(p)] for b in range(p)]


def asymmetric_pairs(rel):
    p = len(rel)
    return [(b + 1, a + 1) for a in range(p) for b in range(p)
            if b < a and rel[b][a] != rel[a][b]]


def offdiagonal_pairs(rel):
    """All (beta, alpha) with beta != alpha and I_beta subset J_alpha: the pairs that
    supply the hypothesis of Lemma 3(2)."""
    p = len(rel)
    return [(b + 1, a + 1) for a in range(p) for b in range(p)
            if a != b and rel[b][a]]


def rel_to_text(rel):
    return " ".join("".join("1" if v else "0" for v in row) for row in rel)


def fmt(s):
    return "{" + ",".join(str(x) for x in sorted(s)) + "}" if s else "{}"


# ---------------------------------------------------------------------------------
# DERIVED E: every row of the paper's table is decided -- either p >= n+1 = 7, or the
# printed certificate pair really does satisfy I_beta subset J_alpha with beta != alpha.
# ---------------------------------------------------------------------------------
def check_certificates(prows):
    bad, by_lemma1, by_lemma2 = [], [], []
    print("  per-row decision (Lemma 3(1) needs p >= n+1 = %d; Lemma 3(2) needs p = n = %d"
          " and one off-diagonal vanishing pair):" % (LEMMA3_BOUND, N))
    for r in prows:
        lab, p, cert = r["label"], r["p"], r["cert"]
        if cert is None:
            ok = p >= LEMMA3_BOUND
            print("    %-5s p=%d  no certificate printed -> Lemma 3(1), p >= %d : %s"
                  % (lab, p, LEMMA3_BOUND, ok))
            by_lemma1.append(lab)
        else:
            b, a = cert
            ok_idx = (b != a) and 1 <= b <= p and 1 <= a <= p
            ok_inc = ok_idx and (r["I"][b - 1] <= r["J"][a - 1])
            ok = (p == N) and ok_idx and ok_inc
            print("    %-5s p=%d  certificate I_%d subset J_%d : I_%d=%s  J_%d=%s -> %s"
                  % (lab, p, b, a, b, fmt(r["I"][b - 1]), a, fmt(r["J"][a - 1]), ok_inc))
            by_lemma2.append(lab)
        if not ok:
            bad.append(lab)
    check("every one of the %d rows of the paper's Table 1 is decided "
          "(certificate verified, or p >= %d)" % (len(prows), LEMMA3_BOUND),
          not bad, "undecided rows = %s" % (bad if bad else "none"))
    check("rows routed to Lemma 3(1) are exactly the 5 rows 14,15,16,17,18",
          sorted(by_lemma1) == ["14", "15", "16", "17", "18"],
          "%d rows: %s" % (len(by_lemma1), sorted(by_lemma1)))
    check("a certificate is printed for a row iff that row has p = 6 = n "
          "(16 rows), and omitted iff p >= 7 (5 rows)",
          all(r["p"] == N for r in prows if r["cert"] is not None)
          and all(r["p"] >= LEMMA3_BOUND for r in prows if r["cert"] is None)
          and len(by_lemma2) == 16 and len(by_lemma1) == 5,
          "p=6 rows with a certificate = %d, p>=7 rows without = %d"
          % (len(by_lemma2), len(by_lemma1)))
    ps = sorted(set(r["p"] for r in prows))
    check("the derived multiset of p values over the table is {6:16, 7:4, 8:1}",
          ps == [6, 7, 8]
          and [len([r for r in prows if r["p"] == v]) for v in (6, 7, 8)] == [16, 4, 1],
          "counts = %s" % ([(v, len([r for r in prows if r["p"] == v])) for v in ps],))
    return not bad


# ---------------------------------------------------------------------------------
# DERIVED F (gate 1): A = A^T forces u_beta^T A u_alpha = u_alpha^T A u_beta, so by (2)
# the relation I_beta subset J_alpha must be symmetric.  Any asymmetry is a defect in
# the tabulated data, not in the paper's use of it.
# ---------------------------------------------------------------------------------
def check_symmetry(prows):
    offenders = []
    print("  relation matrices R[b][a] = [I_b subset J_a], rows b = 1..p:")
    for r in prows:
        rel = relation(r)
        asym = asymmetric_pairs(rel)
        print("    %-5s R = %s   asymmetric pairs: %d"
              % (r["label"], rel_to_text(rel), len(asym)))
        if asym:
            offenders.append((r["label"], asym))
    check("symmetry gate: on all %d rows of the paper's table the relation "
          "I_beta subset J_alpha is symmetric, as A = A^T demands" % len(prows),
          not offenders, "offending rows = %s" % (offenders if offenders else "none"))
    # NOTE: the diagonal of R is I_alpha subset J_alpha, which is true by
    # construction of J (see check_data_hygiene); it is printed for completeness of the
    # relation matrix, not as a falsifiable test.
    diag_ok = all(relation(r)[k][k] for r in prows for k in range(r["n_listed"]))
    check("STRUCTURAL (cannot fail, J was built as I u (J\\I)): the diagonal of every R "
          "is all ones, i.e. I_alpha subset J_alpha, so u_alpha is a zero in every row",
          diag_ok)
    return not offenders


# ---------------------------------------------------------------------------------
# DERIVED F (the O5 row).  The paper argues O5 = B (+) 0 directly and says in a
# footnote that, under its own definition of J_alpha, the O5 row understates J_alpha for
# alpha <= 5, because A e_6 = 0 puts 6 in every J_alpha.  Derived here: as transcribed in
# DATA 2 the row breaks the symmetry gate; adding 6 to every J_alpha repairs it; and the
# repaired row then carries a Lemma 3(2) pair, which is what the paper's direct argument
# concludes.  NOTE ON STRENGTH: DATA 2 is a single hand transcription that has not been
# collated with the published print, so the failure below is a property of the row AS
# TRANSCRIBED HERE.  This program does not establish an error in the published source,
# and no such claim is made; the paper's careful footnote reading is what is reproduced.
# ---------------------------------------------------------------------------------
def check_o5(hrows):
    row = dict((r["label"], r) for r in hrows)["O5"]
    rel = relation(row)
    asym = asymmetric_pairs(rel)
    print("  O5 as transcribed in DATA 2: I = %s" % ", ".join(fmt(s) for s in row["I"]))
    print("  O5 as transcribed in DATA 2: J = %s" % ", ".join(fmt(s) for s in row["J"]))
    print("    R = %s   asymmetric pairs = %s" % (rel_to_text(rel), asym))
    check("the O5 row as transcribed here FAILS the A = A^T symmetry gate, exactly as the "
          "paper's footnote predicts for a row that understates J_alpha (no claim of an "
          "error in the published source is made: this transcription is not collated with "
          "the print, and the paper does not use the row)",
          len(asym) > 0, "asymmetric pairs = %s" % (asym,))
    e6 = frozenset([N])
    six_in_all_printed = all(e6 <= s for s in row["J"])
    check("and the reason is visible: 6 is absent from J_alpha for the alpha where "
          "the row records J_alpha \\ I_alpha = {} (so the row understates J_alpha)",
          not six_in_all_printed,
          "alphas whose transcribed J_alpha omits 6: %s"
          % [k + 1 for k in range(row["n_listed"]) if N not in row["J"][k]])
    rel6 = relation(row, extra=e6)
    asym6 = asymmetric_pairs(rel6)
    print("    after A e_6 = 0 forces 6 in every J_alpha:")
    print("    J = %s" % ", ".join(fmt(s | e6) for s in row["J"]))
    print("    R = %s   asymmetric pairs = %s" % (rel_to_text(rel6), asym6))
    check("O5 repaired by 6 in every J_alpha (justified by A e_6 = 0) satisfies the "
          "symmetry gate", not asym6, "asymmetric pairs = %d" % len(asym6))
    pairs = offdiagonal_pairs(rel6)
    check("the repaired O5 row has p = 6 = n and an off-diagonal vanishing pair, so "
          "Lemma 3(2) gives dim F_A >= %d -- the paper's direct O5 argument"
          % LEMMA3_BOUND,
          row["p"] == N and len(pairs) > 0,
          "p = %d, certificate pairs = %d, e.g. %s" % (row["p"], len(pairs), pairs[:4]))
    involving6 = [(b, a) for (b, a) in pairs if b == N or a == N]
    check("in particular the pair the paper names -- e_6 (support {6}, index 6) against "
          "another minimal zero -- is among them",
          len(involving6) > 0, "pairs involving index 6: %s" % (involving6[:6],))
    return True


# ---------------------------------------------------------------------------------
# DERIVED G: robustness sweep over the source table.  Not needed for the proof (the
# paper only needs one certificate per p = 6 row) but it shows the conclusion does not
# hinge on the particular pair printed: every component is decided with room to spare.
# ---------------------------------------------------------------------------------
def check_robustness(hrows):
    # NOTE: offdiagonal_pairs returns ORDERED pairs, and the relation is
    # symmetric, so every unordered certificate {beta,alpha} is counted twice.  A bar of
    # "at least two ordered pairs" is therefore identical to "at least one" and is not a
    # margin.  The honest quantity is the number of UNORDERED pairs; it is reported, and
    # the minimum over the p = n components is printed rather than hidden.
    weak, margins = [], []
    for r in hrows:
        rel = relation(r)
        pairs = offdiagonal_pairs(rel)
        unord = sorted(set(tuple(sorted(pr)) for pr in pairs))
        route = "p >= %d (Lemma 3(1))" % LEMMA3_BOUND if r["p"] >= LEMMA3_BOUND \
            else "%d unordered certificate pair(s) (Lemma 3(2))" % len(unord)
        print("    %-5s p=%d  ordered pairs = %-3d unordered = %-3d -> %s"
              % (r["label"], r["p"], len(pairs), len(unord), route))
        if r["p"] < LEMMA3_BOUND:
            margins.append((r["label"], len(unord)))
        if not (r["p"] >= LEMMA3_BOUND or len(unord) >= 1):
            weak.append((r["label"], r["p"], len(unord)))
    worst = min(m for _, m in margins) if margins else 0
    thin = [lab for lab, m in margins if m == worst]
    check("robustness: each of HA's %d main components has p >= %d, or at least one "
          "unordered off-diagonal pair certifying Lemma 3(2)" % (len(hrows), LEMMA3_BOUND),
          not weak and worst >= 1,
          "components with none = %s ; smallest unordered margin over the %d components "
          "with p = %d is %d, attained at %s"
          % (weak if weak else "none", len(margins), N, worst, thin))
    decided = [r["label"] for r in hrows
               if r["p"] >= LEMMA3_BOUND or offdiagonal_pairs(relation(r))]
    check("all %d main components of the source table (O5 included, using its row as "
          "transcribed, unrepaired) are decided by Lemma 3" % len(hrows),
          len(decided) == len(hrows), "decided = %d of %d" % (len(decided), len(hrows)))
    return not weak


# ---------------------------------------------------------------------------------
# DERIVED G2: Lemma 3 assumes the minimal zeros SPAN R^6.
# That hypothesis comes from Hildebrand's criterion and is cited, but it has a
# consequence that IS falsifiable on the tabulated supports: the largest rank attainable
# by any family of vectors with supports I_1,...,I_p equals the maximum matching of the
# bipartite support graph, so rank 6 forces a matching of size 6 (a system of distinct
# representatives).  If a transcription slip shrank an I_alpha, this can fail.
# ---------------------------------------------------------------------------------
def max_matching(sets):
    """Kuhn's algorithm: largest matching between the p supports and {1,...,N}."""
    match = {}                                   # coordinate -> index of the support

    def augment(a, seen):
        for k in sorted(sets[a]):
            if k in seen:
                continue
            seen.add(k)
            if k not in match or augment(match[k], seen):
                match[k] = a
                return True
        return False

    size = 0
    for a in range(len(sets)):
        if augment(a, set()):
            size += 1
    return size, match


def check_spanning_necessary(prows, hrows):
    bad = []
    for tag, rows in (("paper", prows), ("HA", hrows)):
        for r in rows:
            size, match = max_matching(r["I"])
            if size != N:
                bad.append((tag, r["label"], size))
            if tag == "paper":
                sdr = ",".join("u_%d->%d" % (match[k] + 1, k)
                               for k in sorted(match))
                print("    %-5s max matching = %d  SDR: %s" % (r["label"], size, sdr))
    check("necessary condition for Lemma 3's spanning hypothesis: the supports of every "
          "row (all %d paper rows and all %d HA rows) admit a system of distinct "
          "representatives, so rank 6 is attainable -- max matching = %d in every row"
          % (len(prows), len(hrows), N),
          not bad, "rows with a smaller matching = %s" % (bad if bad else "none"))
    return not bad


# ---------------------------------------------------------------------------------
# DERIVED H: the paper's equation (2).  Since A u_alpha >= 0 and u_beta >= 0, the sum
#   u_beta^T (A u_alpha) = sum_{k in I_beta} (u_beta)_k (A u_alpha)_k
# has no cancellation, so it vanishes iff (A u_alpha)_k = 0 for every k in I_beta, i.e.
# iff I_beta subset J_alpha.  Checked mechanically on exact nonnegative stand-ins with
# the prescribed supports and zero sets, over three different positive weightings, so
# the outcome cannot depend on the particular positive values chosen.
# ---------------------------------------------------------------------------------
def check_criterion(prows):
    weightings = [
        ("all ones", lambda k: Fraction(1)),
        ("k", lambda k: Fraction(k)),
        ("1/k", lambda k: Fraction(1, k)),
    ]
    bad = 0
    total = 0
    for r in prows:
        p = r["n_listed"]
        for _name, w in weightings:
            # u_alpha: positive exactly on I_alpha;  Au_alpha: zero exactly on J_alpha,
            # positive off it (the generic situation the table records).
            u = [[w(k) if k in r["I"][a] else Fraction(0) for k in range(1, N + 1)]
                 for a in range(p)]
            Au = [[Fraction(0) if k in r["J"][a] else w(k) for k in range(1, N + 1)]
                  for a in range(p)]
            for a in range(p):
                for b in range(p):
                    ip = sum(u[b][i] * Au[a][i] for i in range(N))
                    total += 1
                    if (ip == 0) != (r["I"][b] <= r["J"][a]):
                        bad += 1
            # the vectors must really have the supports/zero sets claimed
            for a in range(p):
                if [k for k in range(1, N + 1) if u[a][k - 1] > 0] != sorted(r["I"][a]):
                    bad += 1
                if [k for k in range(1, N + 1) if Au[a][k - 1] == 0] != sorted(r["J"][a]):
                    bad += 1
    # NOTE: this is an IDENTITY, not a verification.  u_alpha and Au_alpha
    # are defined from I_alpha and J_alpha with strictly positive weights, so
    #   u_beta . Au_alpha = sum_{k in I_beta \ J_alpha} w(k)^2,
    # which vanishes iff I_beta subset J_alpha for ANY data, right or wrong.  No matrix A
    # is constructed anywhere in this program.  What the identity illustrates is the
    # no-cancellation step of the paper's argument; the substantive input to (2) -- that
    # Au >= 0 for a zero u of a copositive A -- is a hand proof and is NOT machine-checked
    # here.  It is listed in the gap ledger printed at the end.
    check("IDENTITY (cannot fail; illustrates eq. (2) rather than testing the paper): on "
          "nonnegative stand-ins with supports I_alpha and zero sets J_alpha, "
          "u_beta^T A u_alpha = 0 <=> I_beta subset J_alpha on all %d ordered pairs "
          "(weightings %s, %d rows)"
          % (total, [nm for nm, _ in weightings], len(prows)),
          bad == 0, "violations = %d" % bad)
    return bad == 0


# ---------------------------------------------------------------------------------
# exact linear algebra over Q (no floating point, no tolerance)
# ---------------------------------------------------------------------------------
SYM_INDEX = [(i, j) for i in range(N) for j in range(i, N)]   # 21 coordinates on S^6


def sym_vec(x):
    """Vectorize the rank-one symmetric matrix x x^T in the 21 coordinates of S^6."""
    return [Fraction(x[i]) * Fraction(x[j]) for (i, j) in SYM_INDEX]


def rank(vectors):
    """Exact rank over Q by Gaussian elimination on a copy of the input."""
    rows = [list(v) for v in vectors]
    if not rows:
        return 0
    width = len(rows[0])
    r = 0
    for col in range(width):
        piv = None
        for k in range(r, len(rows)):
            if rows[k][col] != 0:
                piv = k
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][col]
        rows[r] = [v / pv for v in rows[r]]
        for k in range(len(rows)):
            if k != r and rows[k][col] != 0:
                f = rows[k][col]
                rows[k] = [rows[k][c] - f * rows[r][c] for c in range(width)]
        r += 1
        if r == len(rows):
            break
    return r


# ---------------------------------------------------------------------------------
# DERIVED I: the linear algebra of Lemma 3, exactly.  After the invertible map
# Phi(X) = U^{-1} X U^{-T} the n independent zeros become e_1,...,e_n, so the claim
# "n diagonal rank-one matrices plus one more rank-one direction are independent" is
# what has to hold.  Both branches are checked, together with the failure mode when
# the hypothesis of branch (1) is dropped.
# ---------------------------------------------------------------------------------
def check_lemma3_pairs():
    diag = [sym_vec([1 if k == i else 0 for k in range(N)]) for i in range(N)]
    r0 = rank(diag)
    check("Lemma 3: the n = %d matrices e_i e_i^T are linearly independent (rank %d)"
          % (N, N), r0 == N, "rank = %d" % r0)
    bad = []
    for i in range(N):
        for j in range(i + 1, N):
            x = [0] * N
            x[i] = 1
            x[j] = 1
            r = rank(diag + [sym_vec(x)])
            if r != N + 1:
                bad.append((i + 1, j + 1, r))
    check("Lemma 3(2): for each of the 15 pairs i<j, (e_i+e_j)(e_i+e_j)^T is independent "
          "of the diagonal family, giving rank n+1 = %d" % LEMMA3_BOUND,
          not bad, "pairs failing = %s" % (bad if bad else "none"))
    # dropping the hypothesis: a single-nonzero direction adds nothing
    stuck = []
    for i in range(N):
        for val in (1, 2, -3):
            x = [0] * N
            x[i] = val
            r = rank(diag + [sym_vec(x)])
            if r != N:
                stuck.append((i + 1, val, r))
    check("and the hypothesis is necessary: for c supported on one coordinate, c c^T is "
          "diagonal and the rank stays at %d" % N,
          not stuck, "unexpected ranks = %s" % (stuck if stuck else "none"))
    return not bad and not stuck


def check_lemma3_sweep():
    """Branch (1) of Lemma 3 for every normalized coordinate vector c over a finite
    alphabet: sum(c) = 1 (normalization) and at least two nonzero entries (distinctness)
    must force rank n+1."""
    diag = [sym_vec([1 if k == i else 0 for k in range(N)]) for i in range(N)]
    alphabet = (-1, 0, 1, 2)
    extra = [                                   # rational and larger-entry cases
        [Fraction(1, 2), Fraction(1, 2), 0, 0, 0, 0],
        [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3), 0, 0, 0],
        [Fraction(-3), Fraction(4), 0, 0, 0, 0],
        [Fraction(5), Fraction(-7), Fraction(3), 0, 0, 0],
        [Fraction(1, 6)] * 6,
        [Fraction(7, 10), Fraction(1, 10), Fraction(1, 10), Fraction(1, 10), 0, 0],
        [Fraction(-1, 2), Fraction(3, 2), 0, 0, 0, 0],
        [Fraction(2), Fraction(-1), Fraction(1), Fraction(-1), Fraction(1), Fraction(-1)],
    ]
    tested = 0
    bad = []
    singles = 0
    cands = [list(c) for c in itertools.product(alphabet, repeat=N)] + extra
    for c in cands:
        if sum(Fraction(v) for v in c) != 1:
            continue                      # not a normalized zero
        nz = [k for k in range(N) if c[k] != 0]
        if len(nz) < 2:
            singles += 1
            continue
        tested += 1
        r = rank(diag + [sym_vec(list(c))])
        if r != N + 1:
            bad.append((tuple(str(v) for v in c), r))
    check("Lemma 3(1): every normalized c over {%s}^%d with at least two nonzero entries, "
          "plus %d rational/large-entry cases (%d vectors in all), yields rank n+1 = %d"
          % (",".join(str(a) for a in alphabet), N, len(extra), tested, LEMMA3_BOUND),
          tested > 0 and not bad,
          "tested = %d, failures = %s" % (tested, bad[:4] if bad else "none"))
    check("and the excluded case is exactly c = e_k (normalization forces the single "
          "nonzero entry to be 1, i.e. u_{n+1} = u_k, contradicting distinctness)",
          singles == N, "single-support normalized c found = %d, expected %d"
          % (singles, N))
    return not bad


def check_lemma3_concrete():
    """The same two branches without appealing to Phi: an explicit basis of nonnegative
    normalized vectors of R^6, all 15 sum-pairs for branch (2), and (u_i+u_j)/2 as a
    seventh normalized zero for branch (1).  Also the normalization identity used in the
    proof: 1^T U = 1^T implies that the coordinates of a normalized vector sum to 1."""
    U = []
    for i in range(N - 1):
        U.append([Fraction(1 if k == i else 0) for k in range(N)])
    U.append([Fraction(k + 1, 21) for k in range(N)])      # sums to 21/21 = 1
    colsum_ok = all(sum(u) == 1 for u in U)
    check("explicit basis: all %d columns are nonnegative and normalized (1^T u = 1)"
          % N, colsum_ok and all(all(v >= 0 for v in u) for u in U),
          "column sums = %s" % [str(sum(u)) for u in U])
    base = [sym_vec(u) for u in U]
    r0 = rank(base)
    check("explicit basis: the %d rank-one matrices u_i u_i^T are independent (rank %d)"
          % (N, N), r0 == N, "rank = %d" % r0)
    bad2, bad1, norm_bad = [], [], []
    for i in range(N):
        for j in range(i + 1, N):
            s = [U[i][k] + U[j][k] for k in range(N)]              # u_i + u_j, a zero
            if rank(base + [sym_vec(s)]) != N + 1:
                bad2.append((i + 1, j + 1))
            h = [v / 2 for v in s]                                 # normalized: 1^T h = 1
            if sum(h) != 1:
                norm_bad.append((i + 1, j + 1, str(sum(h))))
            if any(h == U[k] for k in range(N)):                   # must be a NEW zero
                norm_bad.append((i + 1, j + 1, "coincides with some u_k"))
            if rank(base + [sym_vec(h)]) != N + 1:
                bad1.append((i + 1, j + 1))
    check("Lemma 3(2) on the explicit basis: every u_i + u_j supplies the (n+1)-st "
          "independent direction (15 pairs)", not bad2,
          "failures = %s" % (bad2 if bad2 else "none"))
    check("Lemma 3(1) on the explicit basis: (u_i+u_j)/2 is a further normalized zero, "
          "distinct from all u_k, and supplies it too (15 pairs)",
          not bad1 and not norm_bad,
          "failures = %s, normalization slips = %s"
          % (bad1 if bad1 else "none", norm_bad if norm_bad else "none"))
    ident = []
    for c in itertools.product((-1, 0, 1, 2), repeat=N):
        if sum(c) != 1:
            continue
        w = [sum(Fraction(c[t]) * U[t][k] for t in range(N)) for k in range(N)]
        if sum(w) != 1:
            ident.append(c)
    check("normalization identity 1^T U = 1^T: w = Uc has 1^T w = 1^T c, so every "
          "normalized zero has coordinates summing to 1", not ident,
          "counterexamples = %s" % (ident[:3] if ident else "none"))
    return not (bad1 or bad2 or norm_bad or ident)


# ---------------------------------------------------------------------------------
# DERIVED J: the three nonexceptional face dimensions.  In each case the face is
#   F_A = cone{ x x^T : x in R_+^6 \ {0}, x^T A x = 0 },
# so dim F_A = dim span{ x x^T : x >= 0, x^T A x = 0 }.  A lower bound comes from the
# exact rank of an explicit family of zeros, an upper bound from a linear subspace that
# visibly contains every generator; the two meet.
# ---------------------------------------------------------------------------------
def zero_family(indices):
    """xx^T for x = e_j (j in indices) and x = e_j + e_k (j<k in indices)."""
    fam = []
    idx = sorted(indices)
    for j in idx:
        x = [0] * N
        x[j] = 1
        fam.append(sym_vec(x))
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            x = [0] * N
            x[idx[a]] = 1
            x[idx[b]] = 1
            fam.append(sym_vec(x))
    return fam


def check_face_Eii():
    dims, upper, contained = [], [], True
    for i in range(N):
        rest = [k for k in range(N) if k != i]
        fam = zero_family(rest)                       # the zeros of E_ii are {x_i = 0}
        dims.append(rank(fam))
        touching = [t for t, (p, q) in enumerate(SYM_INDEX) if p == i or q == i]
        upper.append(len(SYM_INDEX) - len(touching))  # dim{M in S^6 : row/col i = 0}
        for v in fam:
            if any(v[t] != 0 for t in touching):
                contained = False
    print("    zeros of E_ii are {x : x_i = 0}; derived spans have dimensions %s"
          % dims)
    print("    every generator lies in {M : M row and column i vanish}, of dimension %s"
          % upper)
    check("nonexceptional ray E_ii: dim F = 15 for all %d choices of i "
          "(lower bound = exact rank, upper bound = containing subspace)" % N,
          dims == [15] * N and upper == [15] * N and contained,
          "ranks = %s, upper bounds = %s, containment = %s" % (dims, upper, contained))
    return dims == [15] * N


def check_face_Eij():
    tot, di, dj, inter, upper = [], [], [], [], []
    contained = True
    for i in range(N):
        for j in range(i + 1, N):
            fi = zero_family([k for k in range(N) if k != i])   # hyperplane x_i = 0
            fj = zero_family([k for k in range(N) if k != j])   # hyperplane x_j = 0
            ri, rj, rt = rank(fi), rank(fj), rank(fi + fj)
            di.append(ri)
            dj.append(rj)
            tot.append(rt)
            # the two 15-dimensional subspaces meet in the matrices vanishing on both
            # rows/columns i and j; its dimension is derived by coordinate count
            both = [t for t, (p, q) in enumerate(SYM_INDEX)
                    if p in (i, j) or q in (i, j)]
            inter.append(len(SYM_INDEX) - len(both))
            # upper bound: every generator has (i,j) entry x_i x_j = 0
            pos = SYM_INDEX.index((i, j))
            upper.append(len(SYM_INDEX) - 1)
            for v in fi + fj:
                if v[pos] != 0:
                    contained = False
    ok = (set(di) == {15} and set(dj) == {15} and set(inter) == {10}
          and set(tot) == {20} and set(upper) == {20} and contained)
    print("    zeros of E_ij are {x_i = 0} u {x_j = 0}; over the 15 pairs the derived")
    print("    dimensions are: each hyperplane span %s, their intersection %s, the sum %s"
          % (sorted(set(di)), sorted(set(inter)), sorted(set(tot))))
    check("nonexceptional ray E_ij: dim F = 15 + 15 - 10 = 20 for all 15 pairs "
          "(rank of the union, matched by the upper bound {M : M_ij = 0})",
          ok and all(a + b - c == d for a, b, c, d in zip(di, dj, inter, tot)),
          "spans %s, intersections %s, totals %s, containment %s"
          % (sorted(set(di)), sorted(set(inter)), sorted(set(tot)), contained))
    return ok


def constraint_rows_Ma(a):
    """The 6 x 21 matrix of the linear map S^6 -> R^6, M |-> M a, in the coordinates
    SYM_INDEX.  Its rank r gives dim{M in S^6 : M a = 0} = 21 - r."""
    rows = [[Fraction(0)] * len(SYM_INDEX) for _ in range(N)]
    for t, (p, q) in enumerate(SYM_INDEX):
        rows[p][t] += Fraction(a[q])
        if p != q:
            rows[q][t] += Fraction(a[p])
    return rows


def check_face_aaT():
    ok_all = True
    for a in ([1, 1, 1, 1, 1, -5], [2, -1, 3, -4, 1, -1]):
        mixed = any(v > 0 for v in a) and any(v < 0 for v in a)
        # zeros of a a^T are {x >= 0 : a^T x = 0}: for a_i > 0 > a_j take (-a_j)e_i + a_i e_j
        gens = []
        for i in range(N):
            for j in range(N):
                if a[i] > 0 > a[j]:
                    x = [0] * N
                    x[i] = -a[j]
                    x[j] = a[i]
                    gens.append(x)
        gens = gens + [[u[k] + v[k] for k in range(N)]
                       for u, v in itertools.combinations(gens, 2)]
        nonneg = all(all(v >= 0 for v in x) for x in gens)
        onperp = all(sum(x[k] * a[k] for k in range(N)) == 0 for x in gens)
        fam = [sym_vec(x) for x in gens]
        r = rank(fam)
        cr = rank(constraint_rows_Ma(a))
        upper = len(SYM_INDEX) - cr
        # every generator satisfies (xx^T) a = x (x^T a) = 0, hence lies in that subspace
        inside = all(sum(rowv[t] * f[t] for t in range(len(SYM_INDEX))) == 0
                     for f in fam for rowv in constraint_rows_Ma(a))
        print("    a = %s : %d zero generators, all nonnegative = %s, all in a^perp = %s"
              % (a, len(gens), nonneg, onperp))
        print("      derived rank of span{xx^T} = %d ; dim{M in S^6 : M a = 0} = 21 - %d "
              "= %d ; generators inside it = %s" % (r, cr, upper, inside))
        good = (mixed and nonneg and onperp and inside and r == 15 and upper == 15)
        check("nonexceptional ray a a^T with a = %s (both signs present): dim F = 15"
              % (a,), good, "rank = %d, upper bound = %d" % (r, upper))
        ok_all = ok_all and good
    return ok_all


def aaT_dim(a):
    """dim span{xx^T : x >= 0, a^T x = 0} as an exact rank, together with the upper bound
    21 - rank(M -> Ma).  Returns (rank, upper, all_generators_nonnegative_and_in_a_perp)."""
    gens = []
    for i in range(N):
        for j in range(N):
            if a[i] > 0 > a[j]:
                x = [0] * N
                x[i] = -a[j]
                x[j] = a[i]
                gens.append(x)
    gens = gens + [[u[k] + v[k] for k in range(N)]
                   for u, v in itertools.combinations(gens, 2)]
    good = (all(all(v >= 0 for v in x) for x in gens)
            and all(sum(x[k] * a[k] for k in range(N)) == 0 for x in gens))
    return rank([sym_vec(x) for x in gens]), len(SYM_INDEX) - rank(constraint_rows_Ma(a)), good


def check_face_aaT_sweep():
    """The paper claims dim F_{aa^T} = 15 for EVERY a with both signs; the two hardcoded
    a's above are a sample of size two.  Sweep every mixed-sign pattern in {-1,1}^6."""
    bad, tested = [], 0
    for signs in itertools.product((-1, 1), repeat=N):
        if not (any(s > 0 for s in signs) and any(s < 0 for s in signs)):
            continue
        tested += 1
        r, upper, good = aaT_dim(list(signs))
        if not (good and r == 15 and upper == 15):
            bad.append((signs, r, upper, good))
    check("nonexceptional ray a a^T: dim F = 15 for all %d mixed-sign patterns in "
          "{-1,1}^%d (exact rank of the zero family, matched by 21 - rank(M -> Ma))"
          % (tested, N),
          tested == 2 ** N - 2 and not bad,
          "patterns tested = %d, failures = %s" % (tested, bad[:3] if bad else "none"))
    return not bad


# ---------------------------------------------------------------------------------
# DERIVED K: the arithmetic bookkeeping of the theorem.  Exact integers throughout.
# ---------------------------------------------------------------------------------
def check_bounds():
    # NOTE: everything in this function is bookkeeping on constants fixed by
    # n = 6.  The three BOOKKEEPING checks below compare hardcoded quantities with
    # themselves and cannot be false for any input; they are labelled so that no reader
    # mistakes them for tests of the paper.  PROVENANCE: the closed form (n^2-5n+8)/2 is
    # NOT printed in the paper (which states only 6 <= low(6) <= 7 from [HZ, Lem. 4.7,
    # Thms 4.8-4.9]); it is the general-n form attributed to [HZ, Thm 4.9].  The only
    # content of the first check is that evaluating that closed form at n = 6 reproduces
    # the integer 7 the paper cites from HZ; it is arithmetic on the formula, and it would
    # fail if the formula were mistyped.
    hz_upper_formula = Fraction(N * N - 5 * N + 8, 2)
    check("the closed form (n^2 - 5n + 8)/2 attributed to [HZ, Thm 4.9] (the paper prints "
          "only the integer bound, not this formula) evaluates at n = 6 to the face "
          "dimension %d the paper cites from HZ" % HZ_FACE_DIM,
          hz_upper_formula == HZ_FACE_DIM and hz_upper_formula.denominator == 1,
          "(36 - 30 + 8)/2 = %s" % hz_upper_formula)
    kt_lo, kt_hi = KT_RANGE
    check("BOOKKEEPING (constants; cannot fail): Kostyukova-Tchemisova's n <= low(n) <= n+3 "
          "reads 6 <= low(6) <= 9 at n = 6 and so does not improve on 6 <= low(6) <= 7",
          (kt_lo, kt_hi) == (6, 9) and kt_hi > HZ_FACE_DIM,
          "KT interval = [%d, %d], HZ interval = [%d, %d]"
          % (kt_lo, kt_hi, N, HZ_FACE_DIM))
    interval = list(range(N, HZ_FACE_DIM + 1))
    check("BOOKKEEPING (constants; cannot fail): the interval the theorem has to close is "
          "{6, 7}", interval == [6, 7],
          "candidates before the new lower bound = %s" % interval)
    survivors = [v for v in interval if v >= LEMMA3_BOUND]
    check("BOOKKEEPING (constants; cannot fail): the lower bound dim F >= n+1 = %d removes "
          "6 from that interval and leaves low(6) = 7" % LEMMA3_BOUND,
          survivors == [HZ_FACE_DIM] and LEMMA3_BOUND == HZ_FACE_DIM,
          "survivors = %s" % survivors)
    return True


# ---------------------------------------------------------------------------------
# assembly: put the branches of the case analysis together.  The trichotomy itself
# (nonexceptional / exceptional of order 5 padded by a zero row and column / Cases
# 1--19) is AHD's published classification and is cited, not reproved here; what is
# derived is the dimension bound on each branch and the minimum over them.
# ---------------------------------------------------------------------------------
def check_assembly(prows, hrows):
    rest0 = [k for k in range(N) if k != 0]
    d_ii = rank(zero_family(rest0))
    d_ij = rank(zero_family([k for k in range(N) if k != 0])
                + zero_family([k for k in range(N) if k != 1]))
    gens = []
    a = [1, 1, 1, 1, 1, -5]
    for i in range(N):
        for j in range(N):
            if a[i] > 0 > a[j]:
                x = [0] * N
                x[i] = -a[j]
                x[j] = a[i]
                gens.append(x)
    gens = gens + [[u[k] + v[k] for k in range(N)]
                   for u, v in itertools.combinations(gens, 2)]
    d_aa = rank([sym_vec(x) for x in gens])
    o5 = dict((r["label"], r) for r in hrows)["O5"]
    o5_ok = (o5["p"] == N
             and len(offdiagonal_pairs(relation(o5, extra=frozenset([N])))) > 0)
    rows_ok = all((r["p"] >= LEMMA3_BOUND)
                  or (r["cert"] is not None
                      and r["I"][r["cert"][0] - 1] <= r["J"][r["cert"][1] - 1])
                  for r in prows)
    branches = [("E_ii", d_ii), ("E_ij", d_ij), ("a a^T", d_aa),
                ("O5 = B (+) 0, via Lemma 3(2)", LEMMA3_BOUND if o5_ok else 0),
                ("Cases 1--19 (%d rows), via Lemma 3" % len(prows),
                 LEMMA3_BOUND if rows_ok else 0)]
    for name, val in branches:
        print("    branch %-40s dim F_A >= %s" % (name, val))
    worst = min(v for _, v in branches)
    check("every branch of the case analysis yields dim F_A >= %d, and the minimum over "
          "the branches is exactly %d" % (LEMMA3_BOUND, LEMMA3_BOUND),
          worst == LEMMA3_BOUND and o5_ok and rows_ok,
          "branch minimum = %s" % worst)
    check("therefore low(6) >= 7; with the 7-dimensional maximal face of Holmgren-Zhang "
          "this gives low(6) = 7, the paper's Theorem 1",
          worst == HZ_FACE_DIM, "lower bound %s, upper bound %s" % (worst, HZ_FACE_DIM))
    return True


def print_gap_ledger():
    """Everything between the checks above and Theorem 1 that this program does NOT verify.
    Printed, not checked: a PASS verdict below means the checked facts hold, not that the
    theorem is machine-proved."""
    print("\n" + "-" * 86)
    print("GAP LEDGER -- premises used by the paper that this program does NOT verify")
    print("-" * 86)
    for i, txt in enumerate([
        "Lemma 2, the extreme exposer (every maximal proper face of a proper cone K is "
        "K cap y^perp for some y generating an extreme ray of K^*): hand proof, no finite "
        "computation.",
        "AHD's classification of the extreme rays of COP^6 into nonexceptional / "
        "zero-padded exceptional of order 5 / Cases 1-19 [AHD, Thm 5.1]: cited, not "
        "reproved here (38 published pages of case analysis).  If a further class existed, "
        "every check above would still pass.",
        "Hildebrand's criterion that the minimal zeros of an exceptional extreme matrix "
        "span R^6 [Hildebrand, Thm 4.5]: cited.  It is the SPANNING hypothesis of Lemma 3; "
        "only its necessary consequence (a system of distinct representatives for the "
        "supports) is checked above.",
        "The analytic step that the tabulated component of each piece is the generic one, "
        "that {I_alpha} and p are constant on a piece, and that J_alpha can only GROW off "
        "the generic locus (HA Sec. 2 / AHD Sec. 5).  This is what extends one tabulated "
        "certificate to every matrix of the piece; it is not checkable from support data.",
        "Equation (2) itself: that Au >= 0 for a zero u of a copositive A, hence no "
        "cancellation in u_beta^T A u_alpha.  Two lines by hand; the check above is an "
        "identity on stand-ins and no matrix A is ever constructed.",
        "The O5 branch as the paper argues it (from A e_6 = 0 alone, using no table row). "
        "The O5 checks above instead test the tabulated O5 row and its repair, which is a "
        "consistency check on different data, neither necessary nor sufficient for that "
        "paragraph.  That row is a single hand transcription which has not been collated "
        "with the published print, so its failure of the symmetry gate is a fact about the "
        "transcription and is NOT offered as an erratum in [HA].",
        "The ENTIRE upper bound low(6) <= 7: the integer 7 is taken from Holmgren-Zhang "
        "and no 7-dimensional maximal face is constructed here (the paper exhibits none "
        "either).  If HZ were wrong, the checks above would still pass and only the lower "
        "bound low(6) >= 7 would survive.",
        "Universality of the two finite sweeps: Lemma 3(1) is checked on c over "
        "{-1,0,1,2}^6 plus 8 rationals, and dim F_{aa^T} = 15 on all 62 sign patterns of "
        "{-1,1}^6 plus 2 further a's -- samples of infinite families, not proofs of them.",
        "Fidelity of the two tables to the printed source.  PAPER_TABLE and HA_TABLE were "
        "both entered BY HAND from the same print of [HA, Table 1] (which the paper's "
        "Table 1 reproduces), in different encodings and read by different parsers.  Their "
        "agreement, checked above, excludes divergent typing and parser-specific slips; it "
        "does NOT exclude a misreading common to both, and no automated extraction from "
        "the published table was performed.",
    ]):
        print(" %d. %s" % (i + 1, txt))


def main():
    print("=" * 86)
    print("verify.py -- low(6) = 7 for the 6x6 completely positive cone")
    print("n = %d, Lemma 3 bound n+1 = %d, exact arithmetic only (int / Fraction)"
          % (N, LEMMA3_BOUND))
    print("=" * 86)

    prows = build_rows(PAPER_TABLE, "paper")
    hrows = build_rows(HA_TABLE, "ha")

    print("\n--- 1. the paper's Table 1: parse, reconstruct J_alpha, derive p ---")
    check_data_hygiene(prows)

    print("\n--- 2. the paper's Table 1 against Table 1 of Hildebrand-Afonin ---")
    check_against_source(prows, hrows)

    print("\n--- 3. equation (2): the support criterion for u_beta^T A u_alpha = 0 ---")
    check_criterion(prows)

    print("\n--- 4. the 21 rows: certificates and the routing to Lemma 3(1)/3(2) ---")
    check_certificates(prows)

    print("\n--- 5. sanity gate: A = A^T forces the relation to be symmetric ---")
    check_symmetry(prows)

    print("\n--- 6. the O5 row of the source table (the paper argues this branch "
          "directly) ---")
    check_o5(hrows)

    print("\n--- 7. robustness over all 22 components of the source table ---")
    check_robustness(hrows)

    print("\n--- 7b. necessary condition for Lemma 3's spanning hypothesis (matchings) ---")
    check_spanning_necessary(prows, hrows)

    print("\n--- 8. Lemma 3: exact linear independence of the rank-one directions ---")
    check_lemma3_pairs()
    check_lemma3_sweep()
    check_lemma3_concrete()

    print("\n--- 9. the three nonexceptional face dimensions (15, 20, 15) ---")
    check_face_Eii()
    check_face_Eij()
    check_face_aaT()
    check_face_aaT_sweep()

    print("\n--- 10. bookkeeping of the two bounds ---")
    check_bounds()

    print("\n--- 11. assembly of the case analysis ---")
    check_assembly(prows, hrows)

    print_gap_ledger()

    total = len(RESULTS)
    failed = [lab for lab, ok in RESULTS if not ok]
    struct = len([lab for lab, _ in RESULTS
                  if lab.startswith("STRUCTURAL") or lab.startswith("IDENTITY")
                  or lab.startswith("BOOKKEEPING")])
    print("\n" + "=" * 86)
    print("%d of the %d checks are labelled STRUCTURAL / IDENTITY / BOOKKEEPING: they hold "
          "for any input\nand carry no falsification power; the remaining %d are the tests "
          "of the paper." % (struct, total, total - struct))
    if failed:
        for lab in failed:
            print("FAILED CHECK: " + lab)
        print("VERDICT: %d OF %d CHECKS FAILED" % (len(failed), total))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % total)
    return 0


if __name__ == "__main__":
    sys.exit(main())
