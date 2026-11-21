#!/usr/bin/env python3
"""Regression checks for the bounded, reproducible C++ DBLP generator."""

import argparse
import json
from pathlib import Path
import subprocess
import tempfile


parser = argparse.ArgumentParser()
parser.add_argument("--bin", type=Path, required=True)
parser.add_argument("--generator", type=Path, required=True)
parser.add_argument("--cplex-bin", type=Path)
args = parser.parse_args()
generator_bin = args.generator.resolve(strict=True)


def payload(stdout: str) -> dict:
    values = [json.loads(line[7:]) for line in stdout.splitlines()
              if line.startswith("RESULT ")]
    if len(values) != 1:
        raise AssertionError(f"expected exactly one RESULT record, got {len(values)}")
    return values[0]


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    target = root / "generated.ddt"
    command = [
        str(generator_bin), "--output", str(target),
        "--n1", "3", "--n2", "3",
        "--m-probability", "5", "--m-value", "4",
        "--seed", "20260911", "--overwrite",
    ]
    first = subprocess.run(command, check=True, capture_output=True, text=True)
    if not any(line.startswith("GENERATED ") for line in first.stdout.splitlines()):
        raise AssertionError("generator did not emit metadata")
    original = target.read_bytes()
    subprocess.run(command, check=True, capture_output=True, text=True)
    if target.read_bytes() != original:
        raise AssertionError("same generator seed is not byte-for-byte reproducible")

    text = target.read_text(encoding="utf-8")
    if text.count("\n3 ") != 9:
        raise AssertionError("generated general-row count does not match mP+mV")
    subprocess.run([str(args.bin), str(target), "--validate-only"],
                   check=True, capture_output=True, text=True)
    highs = payload(subprocess.run(
        [str(args.bin), str(target), "--single", "--method=odba"],
        check=True, capture_output=True, text=True).stdout)
    if highs.get("status") != "optimal":
        raise AssertionError("generated instance did not solve with ODBA")

    if args.cplex_bin:
        cplex = payload(subprocess.run(
            [str(args.cplex_bin), str(target), "--single", "--cplex",
             "--cplex-threads=1"],
            check=True, capture_output=True, text=True).stdout)
        if cplex.get("status") != "optimal":
            raise AssertionError("generated instance did not solve with CPLEX")
        error = max(abs(highs[key] - cplex[key]) for key in ("min", "max"))
        if error > 5e-5:
            raise AssertionError(f"generated-instance HiGHS/CPLEX disagreement: {error}")

    too_large = subprocess.run([
        str(generator_bin), "--output", str(root / "too_large.ddt"),
        "--n1", "31", "--n2", "30", "--m-probability", "0",
        "--m-value", "0"], capture_output=True, text=True)
    if too_large.returncode == 0 or "max total dimension" not in too_large.stderr:
        raise AssertionError(
            "generator did not reject 61 dimensions for the expected reason: "
            f"binary={generator_bin} rc={too_large.returncode} "
            f"stdout={too_large.stdout!r} stderr={too_large.stderr!r}")

    extended = root / "extended.ddt"
    subprocess.run([
        str(generator_bin), "--output", str(extended),
        "--n1", "31", "--n2", "30", "--m-probability", "1",
        "--m-value", "1", "--max-total-dimension", "100",
        "--seed", "20260912"], check=True, capture_output=True, text=True)
    subprocess.run([str(args.bin), str(extended), "--validate-only"],
                   check=True, capture_output=True, text=True)

    over_budget = subprocess.run([
        str(generator_bin), "--output", str(root / "over_budget.ddt"),
        "--n1", "10", "--n2", "10", "--m-probability", "10",
        "--m-value", "10", "--max-rlt-rows", "10"],
        capture_output=True, text=True)
    if over_budget.returncode == 0 or "RLT exceeds" not in over_budget.stderr:
        raise AssertionError("generator ignored the RLT-size gate")

print("bounded C++ DBLP generator passed")
