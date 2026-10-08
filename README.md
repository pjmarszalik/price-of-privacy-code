# The Price of Privacy: code and logs

Code for the paper "The Price of Privacy: Randomness Complexity of
Graph-Based Multi-Secret Sharing" by Piotr Marszalik (arXiv: http://arxiv.org/abs/2610.10002,
IACR ePrint: TODO).

In the short (IEEE) version of the paper, the search is Table IV and
Appendix E ("Exact Computation for the K4 Classification"); the programs
used there are `privacy/normal_form_search.cpp`, `run_exact_search.py` and
`k4_geometry.py`. `privacy/lp_candidates.py` belongs to the statistical-privacy
section of the longer version.

This artifact contains the programs behind the computer-assisted steps of the
paper and their outputs. All arithmetic is exact (integers only).

## Contents

| File | What it checks | Where in the paper |
|---|---|---|
| `privacy/normal_form_search.cpp` | Exhaustive search over the normal forms of three-state schemes on $K_4$ (and $K_3$). For each vertex it enumerates the local labellings of the three source states, forms the forced edge classes by union-find, and accepts exactly when every vertex separates its two bit classes. It does not assume linear shares, binary share alphabets or XOR reconstruction. | Normal form theorem and Table "Exhaustive three-state search on the six orbits" (Appendix on three states on $K_4$) |
| `run_exact_search.py` | Compiles the search, runs it, saves the log and compares every line with that table. Also runs a calibration on four published three-party domains (PARITY, F11, F13, F14): only PARITY admits three states, in strata H and Q. | same |
| `privacy/lp_candidates.py` | Enumerates all C(26,3) = 2600 triples of nonzero rows in {-1,0,1}^3, solves the vertex equations of the leakage linear program in exact rational arithmetic, and keeps p > 0, 0 < s < 1/5. Exactly four weight vectors remain: (2,3,4)/9, (1,2,4)/7, (2,2,3)/7, (1,2,3)/6. | Lemma "Four weight candidates" (Appendix on statistical privacy) |
| `k4_geometry.py` | The 896 odd affine tetrahedra of $\mathbb{F}_2^4$ and their six orbits (sizes 64, 64, 4 x 192); the extension table for $T_6$; a scan of all 65,535 nonempty domains (exactly 80 contain an $\mathcal{O}_6$ tetrahedron and no other obstruction: the orbits of $T_6$ and $C_5$, sizes 64 and 16); the explicit $\mathbb{Z}_3$ scheme on $C_5$; the universal four-state scheme on $K_4$. | Main theorem, extension table, Appendix on data |
| `logs/` | Outputs of the runs below. | |

Strata of source weights used by the search: U = (1,1,1) and Q = (2,1,1)
are the two cases left by the weight reduction; G, E, H are used only for the
calibration.

## Requirements

A C++17 compiler (g++ or clang++) and Python 3.9 or newer. No other libraries.

## How to run

```
python3 run_exact_search.py     # about 5 s; ends with "Table tab:search reproduced: PASS"
python3 k4_geometry.py          # about 30 s; ends with three lines "... PASS"
python3 privacy/lp_candidates.py   # under 1 s; ends with "... PASS"
```

Equivalent manual compilation:

```
c++ -O2 -std=c++17 privacy/normal_form_search.cpp -o privacy/normal_form_search
./privacy/normal_form_search --paper
```

Other modes: `--proof` (only the logically required cases O2-O5, strata U
and Q), `--k4` (all six orbits), `--calibration` (the $K_3$ domains), or a
single instance such as `./privacy/normal_form_search O2 U`.

## Expected output

`logs/normal_form_search.log`: orbits O1-O5 are UNSAT in U and Q, O6 is SAT in
U (the scheme on $T_6$) and UNSAT in Q; the counts after `counts=` are the
numbers of local configurations at vertices 1-4, as in the table.
`logs/k4_geometry.log`: 896 tetrahedra, orbit sizes, 80 candidates, three PASS
lines.
`logs/lp_candidates.log`: the four weight vectors with leakage 1/9, 1/7, 1/7, 1/6.

The programs are not formally verified; their correct execution is part of
the proof of the computer-assisted steps.
