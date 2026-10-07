#!/usr/bin/env python3
"""Exact enumeration behind Lemma "Four weight candidates" (Appendix on
statistical privacy).

For fixed share maps of a three-state scheme with weights p = (p1, p2, p3),
every probability difference in a privacy comparison is u . p with
u in {-1,0,1}^3.  The best leakage is the linear program

    minimize s  subject to  u . p <= s (all rows u),  p1 + p2 + p3 = 1,  p >= 0.

At an optimal vertex with no zero weight, three independent rows are tight:
A p = s 1 together with p1 + p2 + p3 = 1.  This script enumerates all
C(26,3) = 2600 triples of distinct nonzero rows, solves the 4x4 system in
exact rational arithmetic, and keeps the solutions with p > 0 and
0 < s < 1/5.  The paper states that exactly four weight vectors remain (up to
permuting the states):

    (2,3,4)/9 with s = 1/9,   (1,2,4)/7 with s = 1/7,
    (2,2,3)/7 with s = 1/7,   (1,2,3)/6 with s = 1/6.

It also checks the claim used in the proof that A is invertible at every
such vertex, so that p = s A^{-1} 1 and s = 1 / (1^T A^{-1} 1).
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product

F = Fraction
BOUND = F(1, 5)


def solve(M: list[list[Fraction]], rhs: list[Fraction]) -> list[Fraction] | None:
    """Gauss-Jordan elimination over Q; None if M is singular."""
    n = len(M)
    A = [row[:] + [b] for row, b in zip(M, rhs)]
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col] != 0), None)
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        inv = 1 / A[col][col]
        A[col] = [x * inv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [A[r][n] for r in range(n)]


def det3(A: list[tuple[int, int, int]]) -> int:
    (a, b, c), (d, e, f), (g, h, i) = A
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def main() -> None:
    rows = [u for u in product((-1, 0, 1), repeat=3) if any(u)]
    assert len(rows) == 26
    triples = list(combinations(rows, 3))
    assert len(triples) == 2600

    found: dict[tuple[tuple[Fraction, ...], Fraction], int] = {}
    for A in triples:
        # Unknowns p1, p2, p3, s:  u . p - s = 0 for the three rows, sum p = 1.
        M = [[F(x) for x in u] + [F(-1)] for u in A] + [[F(1), F(1), F(1), F(0)]]
        sol = solve(M, [F(0), F(0), F(0), F(1)])
        if sol is None:
            continue
        p, s = sol[:3], sol[3]
        if all(x > 0 for x in p) and 0 < s < BOUND:
            assert det3(list(A)) != 0, ("singular A at a vertex", A)
            key = (tuple(sorted(p)), s)
            found[key] = found.get(key, 0) + 1

    expected = {
        ((F(2, 9), F(3, 9), F(4, 9)), F(1, 9)),
        ((F(1, 7), F(2, 7), F(4, 7)), F(1, 7)),
        ((F(2, 7), F(2, 7), F(3, 7)), F(1, 7)),
        ((F(1, 6), F(2, 6), F(3, 6)), F(1, 6)),
    }
    print("nonzero rows:", len(rows), " row triples:", len(triples))
    for (p, s), n in sorted(found.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        print(f"weights ({', '.join(map(str, p))})  leakage {s}  from {n} triples")
    assert set(found) == expected, sorted(found)
    print("Lemma catalogue (four weight candidates below 1/5): PASS")


if __name__ == "__main__":
    main()
