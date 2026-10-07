#!/usr/bin/env python3
"""Independent exact audit for the finite K4 classification data.

This script does not search for arbitrary three-state schemes.  It verifies:
  * the 896 odd affine tetrahedra;
  * their six Gamma_4 orbits and orbit sizes;
  * the twelve-row extension table and the two candidate orbits T6 and C5;
  * the explicit uniform-Z3 scheme on C5 and its restriction to T6;
  * the universal four-state K4 construction.
All arithmetic is exact integer arithmetic.
"""
from __future__ import annotations
from collections import Counter
from itertools import combinations, permutations, product
from typing import Iterable, Sequence

Point = tuple[int, int, int, int]


def xor(x: Point, y: Point) -> Point:
    return tuple(a ^ b for a, b in zip(x, y))  # type: ignore[return-value]


def dot(x: Point, y: Point) -> int:
    return sum(a * b for a, b in zip(x, y)) & 1


def wt(x: Point) -> int:
    return sum(x)


def dist(x: Point, y: Point) -> int:
    return wt(xor(x, y))


def point(s: str) -> Point:
    assert len(s) == 4 and set(s) <= {"0", "1"}
    return tuple(int(c) for c in s)  # type: ignore[return-value]


def mask_of(S: Iterable[Point]) -> int:
    ans = 0
    for x in S:
        i = int("".join(map(str, x)), 2)
        ans |= 1 << i
    return ans


def set_of(mask: int) -> tuple[Point, ...]:
    return tuple(point(f"{i:04b}") for i in range(16) if (mask >> i) & 1)


def spectrum(S: Sequence[Point]) -> tuple[int, ...]:
    return tuple(sorted(dist(x, y) for x, y in combinations(S, 2)))


def rank_f2(vectors: Sequence[Point]) -> int:
    rows = [sum(bit << (3 - i) for i, bit in enumerate(v)) for v in vectors]
    rank = 0
    for col in range(3, -1, -1):
        pivot = next((j for j in range(rank, len(rows)) if (rows[j] >> col) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for j in range(len(rows)):
            if j != rank and ((rows[j] >> col) & 1):
                rows[j] ^= rows[rank]
        rank += 1
    return rank


def relation_space(S: Sequence[Point]) -> tuple[Point, ...]:
    x0 = S[0]
    diffs = [xor(x, x0) for x in S]
    cube = [point(f"{i:04b}") for i in range(16)]
    return tuple(a for a in cube if all(dot(a, d) == 0 for d in diffs))


def odd_affine_tetrahedron(S: Sequence[Point]) -> bool:
    if len(S) != 4:
        return False
    x0 = S[0]
    diffs = [xor(x, x0) for x in S[1:]]
    if rank_f2(diffs) != 3:
        return False
    rel = [a for a in relation_space(S) if wt(a) > 0]
    return len(rel) == 1 and wt(rel[0]) % 2 == 1


def gamma_images(mask: int) -> set[int]:
    S = set_of(mask)
    out: set[int] = set()
    for pi in permutations(range(4)):
        for b_int in range(16):
            b = point(f"{b_int:04b}")
            img = []
            for x in S:
                y = tuple(x[pi[i]] ^ b[i] for i in range(4))
                img.append(y)  # type: ignore[arg-type]
            out.add(mask_of(img))
    return out


def verify_scheme(
    D: Sequence[Point], states: Sequence[tuple[int, ...]], probs: Sequence[int]
) -> None:
    # states has one 6-tuple per (x,r), ordered x-major.
    assert len(states) == len(D) * len(probs)
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for v in range(4):
        incident = [i for i, e in enumerate(edges) if v in e]
        laws: dict[int, Counter[tuple[int, ...]]] = {}
        for b in (0, 1):
            common = None
            for xi, x in enumerate(D):
                if x[v] != b:
                    continue
                law: Counter[tuple[int, ...]] = Counter()
                for r, p in enumerate(probs):
                    view = tuple(states[xi * len(probs) + r][j] for j in incident)
                    law[view] += p
                if common is None:
                    common = law
                else:
                    assert law == common, (v, b, x, law, common)
            if common is not None:
                laws[b] = common
        if 0 in laws and 1 in laws:
            assert set(laws[0]).isdisjoint(laws[1]), (v, laws)


def c5_shares(x: Point, r: int) -> tuple[int, ...]:
    x1, x2, x3, _x4 = x
    return (
        r % 3,
        0,
        (r - x1) % 3,
        (r - x2) % 3,
        0,
        (r - x2 - x3) % 3,
    )


def universal_k4_shares(x: Point, r: int, s: int) -> tuple[int, ...]:
    x1, x2, x3, x4 = x
    return (
        r,
        r ^ x1,
        s,
        s ^ x4 ^ x2,
        s ^ x4,
        r ^ x1 ^ x3,
    )


def main() -> None:
    cube = [point(f"{i:04b}") for i in range(16)]
    tetra = [tuple(S) for S in combinations(cube, 4) if odd_affine_tetrahedron(S)]
    assert len(tetra) == 896

    remaining = {mask_of(S) for S in tetra}
    orbits: list[set[int]] = []
    while remaining:
        m = min(remaining)
        orb = gamma_images(m) & remaining
        orbits.append(orb)
        remaining -= orb
    orbits.sort(key=lambda O: (spectrum(set_of(next(iter(O)))), min(O)))
    sizes = sorted(len(O) for O in orbits)
    assert sizes == [64, 64, 192, 192, 192, 192]
    assert len(orbits) == 6
    spectra = sorted(spectrum(set_of(next(iter(O)))) for O in orbits)
    expected_spectra = sorted(
        [
            (1, 1, 1, 2, 2, 2),
            (1, 1, 1, 2, 2, 3),
            (1, 1, 2, 2, 2, 3),
            (1, 2, 2, 2, 3, 3),
            (1, 2, 2, 3, 3, 3),
            (2, 2, 2, 3, 3, 3),
        ]
    )
    assert spectra == expected_spectra

    reps = {
        "O1": tuple(map(point, ["0000", "0001", "0010", "0100"])),
        "O2": tuple(map(point, ["0000", "0001", "0011", "0100"])),
        "O3": tuple(map(point, ["0000", "0010", "0011", "0101"])),
        "O4": tuple(map(point, ["0000", "0011", "0101", "1000"])),
        "O5": tuple(map(point, ["0000", "0101", "1000", "1011"])),
        "O6": tuple(map(point, ["0000", "0111", "1001", "1010"])),
    }
    spec_to_name = {spectrum(S): name for name, S in reps.items()}
    assert len(spec_to_name) == 6
    orbit_sizes = {}
    for name, S in reps.items():
        imgs = gamma_images(mask_of(S))
        orbit_sizes[name] = len(imgs)
        assert all(odd_affine_tetrahedron(set_of(m)) for m in imgs)
    assert orbit_sizes == {"O1": 64, "O2": 192, "O3": 192, "O4": 192, "O5": 192, "O6": 64}

    O_masks = {name: gamma_images(mask_of(S)) for name, S in reps.items()}

    # Verify the twelve-row finite extension table used in the paper.
    T6 = reps["O6"]
    extension_witnesses = {
        "0001": ("1001", "O5"),
        "0010": ("1010", "O5"),
        "0011": ("1010", "O4"),
        "0100": ("1010", "O5"),
        "0101": ("1010", "O4"),
        "0110": ("1001", "O4"),
        "1000": ("0111", "O1"),
        "1011": ("1010", "O5"),
        "1101": ("1010", "O5"),
        "1110": ("1001", "O5"),
        "1111": ("0000", "O4"),
    }
    for added, (omitted, orbit_name) in extension_witnesses.items():
        four_set = (set(T6) - {point(omitted)}) | {point(added)}
        assert mask_of(four_set) in O_masks[orbit_name]
    safe_extension = set(T6) | {point("1100")}
    safe_fours = {mask_of(S) for S in combinations(safe_extension, 4)}
    assert not any(safe_fours & O_masks[name] for name in ["O1", "O2", "O3", "O4", "O5"])
    assert len(safe_fours & O_masks["O6"]) == 4

    candidates = []
    for mask in range(1, 1 << 16):
        pts = set_of(mask)
        fours = {mask_of(S) for S in combinations(pts, 4)}
        if any(fours & O_masks[name] for name in ["O1", "O2", "O3", "O4", "O5"]):
            continue
        if fours & O_masks["O6"]:
            candidates.append(mask)
    assert len(candidates) == 80
    counts = Counter(len(set_of(m)) for m in candidates)
    assert counts == Counter({4: 64, 5: 16})
    candidate_orbits: list[set[int]] = []
    left = set(candidates)
    while left:
        m = min(left)
        orb = gamma_images(m) & left
        candidate_orbits.append(orb)
        left -= orb
    assert sorted(map(len, candidate_orbits)) == [16, 64]

    C5 = tuple(map(point, ["0000", "0111", "1001", "1010", "1100"]))
    assert set(T6) < set(C5)
    assert all((x[0] - x[1] - x[2] - x[3]) % 3 == 0 for x in C5)
    assert {x for x in cube if (x[0] - x[1] - x[2] - x[3]) % 3 == 0} == set(C5)
    states = [c5_shares(x, r) for x in C5 for r in range(3)]
    verify_scheme(C5, states, [1, 1, 1])
    assert relation_space(T6) == (point("0000"), point("1011"))
    assert relation_space(C5) == (point("0000"),)

    full_cube = tuple(cube)
    states4 = [universal_k4_shares(x, r, s) for x in full_cube for r in (0, 1) for s in (0, 1)]
    verify_scheme(full_cube, states4, [1, 1, 1, 1])

    print("odd affine tetrahedra:", len(tetra))
    print("orbit sizes:", orbit_sizes)
    print("candidate labelled domains:", len(candidates), dict(sorted(counts.items())))
    print("candidate orbit sizes:", sorted(map(len, candidate_orbits)))
    print("T6 extension table: PASS")
    print("C5/T6 explicit Z3 scheme: PASS")
    print("universal K4 four-state scheme: PASS")


if __name__ == "__main__":
    main()
