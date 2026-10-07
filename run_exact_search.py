#!/usr/bin/env python3
"""Driver for privacy/normal_form_search.cpp.

Compiles the exact normal-form search, runs it in --paper mode, saves the
output to logs/normal_form_search.log and checks every line against
Table "Exhaustive three-state search on the six orbits" of the paper
(local-configuration counts and SAT/UNSAT for strata U and Q), and against
the four published K3 calibration domains.

Usage:
    python3 run_exact_search.py            # compile, run, check
    python3 run_exact_search.py --no-build # reuse an existing binary

The C++ compiler is taken from $CXX, otherwise c++, g++ or clang++.
"""
from __future__ import annotations

import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "privacy" / "normal_form_search.cpp"
BIN = HERE / "privacy" / ("normal_form_search.exe" if os.name == "nt" else "normal_form_search")
LOGS = HERE / "logs"

# Table tab:search of the paper: orbit -> (Q counts, Q result, U counts, U result).
TABLE = {
    "O1": ("33,85,85,85", "UNSAT", "298,320,320,320", "UNSAT"),
    "O2": ("33,85,85,81", "UNSAT", "298,320,320,256", "UNSAT"),
    "O3": ("33,85,81,81", "UNSAT", "298,320,256,256", "UNSAT"),
    "O4": ("85,85,85,81", "UNSAT", "320,320,320,256", "UNSAT"),
    "O5": ("81,85,85,81", "UNSAT", "256,320,320,256", "UNSAT"),
    "O6": ("81,85,81,81", "UNSAT", "256,320,256,256", "SAT"),
}
# K3 calibration: only PARITY admits three states, in strata H and Q.
K3_NAMES = ("PARITY", "F11", "F13", "F14")
K3_STRATA = ("G", "E", "U", "H", "Q")

LINE = re.compile(r"^(\w+) ([A-Z]) (SAT|UNSAT) counts=([\d,]+) ")


def find_compiler() -> str:
    for cand in (os.environ.get("CXX"), "c++", "g++", "clang++"):
        if cand and shutil.which(cand):
            return cand
    sys.exit("No C++ compiler found. Install g++ or clang++ (macOS: xcode-select --install).")


def build() -> str:
    cxx = find_compiler()
    cmd = [cxx, "-O2", "-std=c++17", str(SRC), "-o", str(BIN)]
    print("compiling:", " ".join([cxx, "-O2", "-std=c++17",
                                  "privacy/normal_form_search.cpp", "-o",
                                  "privacy/" + BIN.name]))
    subprocess.run(cmd, check=True)
    ver = subprocess.run([cxx, "--version"], capture_output=True, text=True).stdout
    return ver.splitlines()[0] if ver else cxx


def check(output: str) -> list[str]:
    errors: list[str] = []
    seen: dict[tuple[str, str], tuple[str, str]] = {}
    for line in output.splitlines():
        m = LINE.match(line)
        if m:
            seen[(m.group(1), m.group(2))] = (m.group(4), m.group(3))
    for orbit, (qc, qr, uc, ur) in TABLE.items():
        for stratum, counts, result in (("Q", qc, qr), ("U", uc, ur)):
            got = seen.get((orbit, stratum))
            if got != (counts, result):
                errors.append(f"{orbit} {stratum}: expected {counts} {result}, got {got}")
    for name in K3_NAMES:
        for stratum in K3_STRATA:
            want = "SAT" if name == "PARITY" and stratum in ("H", "Q") else "UNSAT"
            got = seen.get((name, stratum))
            if got is None or got[1] != want:
                errors.append(f"{name} {stratum}: expected {want}, got {got}")
    if "ALL CHECKS PASS" not in output:
        errors.append("program did not print ALL CHECKS PASS")
    return errors


def main() -> None:
    LOGS.mkdir(exist_ok=True)
    compiler = "prebuilt binary"
    if "--no-build" not in sys.argv[1:]:
        compiler = build()
    print("running: privacy/" + BIN.name + " --paper")
    proc = subprocess.run([str(BIN), "--paper"], capture_output=True, text=True)
    out = proc.stdout + proc.stderr
    header = (f"# compiler: {compiler}\n"
              f"# python: {platform.python_version()}, system: {platform.system()} {platform.machine()}\n"
              f"# command: privacy/{BIN.name} --paper\n")
    (LOGS / "normal_form_search.log").write_text(header + out)
    print(out, end="")
    errors = check(out) if proc.returncode == 0 else [f"exit code {proc.returncode}"]
    if errors:
        print("FAIL")
        for e in errors:
            print("  " + e)
        sys.exit(1)
    print("Table tab:search reproduced: PASS")
    print("log written to logs/normal_form_search.log")


if __name__ == "__main__":
    main()
