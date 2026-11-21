#!/usr/bin/env python3
"""Exercise compact Group-B conversion and the new sparse-row DDT grammar."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


parser = argparse.ArgumentParser()
parser.add_argument("--bin", type=Path, required=True)
parser.add_argument("--converter", type=Path, required=True)
parser.add_argument("--cplex-bin", type=Path)
args = parser.parse_args()

fixture = """2 2
2 2
1: [0.1 0.9] [0 1]
2: [0.1 0.9] [0 1]
3: [0.1 0.9] [0 1]
4: [0.1 0.9] [0 1]
1: P1+P2-P3<1
2: P3+P4-P1-P2<0
1: V1+V2-V3<2
2: V3+V4-V1-V2<2
"""

with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    source = root / "sample.m"
    output = root / "out"
    source.write_text(fixture, encoding="utf-8")
    subprocess.run([sys.executable, str(args.converter), str(source), "--output-root", str(output)],
                   check=True, capture_output=True, text=True)
    target = output / "sample.ddt"
    text = target.read_text(encoding="utf-8")
    if text.count("\n3 ") != 4:
        raise AssertionError("expected four general sparse DDT rows")
    subprocess.run([sys.executable, str(args.converter), str(source), "--output-root", str(output), "--check"],
                   check=True, capture_output=True, text=True)
    validation = subprocess.run([str(args.bin), str(target), "--validate-only"],
                                check=True, capture_output=True, text=True)
    validation_payloads = [json.loads(line[6:]) for line in validation.stdout.splitlines()
                           if line.startswith("VALID ")]
    if len(validation_payloads) != 1 or validation_payloads[0].get("validation_scope") != "syntax_only":
        raise AssertionError("general sparse rows failed syntax validation")
    # Syntax is not a feasibility proof. Explicitly construct and validate
    # the ORIGINAL-model focus before declaring contraction admissible.
    focal_validation = subprocess.run([str(args.bin),str(target),
        "--sequential-directions",f"--prepare-focal={output / 'focus.txt'}"],
        capture_output=True,text=True)
    if focal_validation.returncode:
        raise AssertionError((focal_validation.stdout,focal_validation.stderr))
    process = subprocess.run([str(args.bin), str(target), "--single", "--method=odba"],
                             capture_output=True, text=True)
    if process.returncode:
        trace_env = dict(os.environ, DECIDEIT_TRACE="1")
        trace = subprocess.run([str(args.bin), str(target), "--single", "--method=odba"],
                               capture_output=True, text=True, env=trace_env)
        highs = subprocess.run([str(args.bin), str(target), "--single", "--method=odba",
                                "--lp-backend=highs"], capture_output=True, text=True,
                               env=trace_env)
        raise AssertionError(
            f"converted sparse-row instance exited {process.returncode}: "
            f"stdout={process.stdout!r} stderr={process.stderr!r}; "
            f"trace_stdout={trace.stdout!r} trace_stderr={trace.stderr!r}; "
            f"highs_returncode={highs.returncode} highs_stdout={highs.stdout!r} "
            f"highs_stderr={highs.stderr!r}")
    payloads = [json.loads(line[7:]) for line in process.stdout.splitlines()
                if line.startswith("RESULT ")]
    if len(payloads) != 1 or payloads[0].get("status") != "optimal":
        raise AssertionError("converted sparse-row instance did not solve")
    path_process = subprocess.run(
        [str(args.bin), str(target), "5", "0", "--method=odba"],
        check=True, capture_output=True, text=True)
    path_payloads = [json.loads(line[7:]) for line in path_process.stdout.splitlines()
                     if line.startswith("RESULT ")]
    if (len(path_payloads) != 6 or
            sum(not level.get("direct_point") for level in path_payloads) != 5 or
            not all(level.get("status") == "optimal" for level in path_payloads)):
        raise AssertionError("general sparse-row contraction path did not solve ten DBLPs")
    inclusive_payloads = path_payloads
    inclusive_solved = [level for level in inclusive_payloads
                        if not level.get("direct_point")]
    interior_process = subprocess.run(
        [str(args.bin), str(target), "5", "0", "--method=odba",
         "--contraction-path=interior"],
        check=True, capture_output=True, text=True)
    interior_payloads = [json.loads(line[7:])
                         for line in interior_process.stdout.splitlines()
                         if line.startswith("RESULT ")]
    if (len(interior_payloads) != 5 or
            sum(not level.get("direct_point") for level in interior_payloads) != 4):
        raise AssertionError("interior contraction path did not solve eight DBLPs")
    if (inclusive_solved[0].get("level") != 0 or
            not inclusive_solved[0].get("original_problem")):
        raise AssertionError("inclusive level zero is not marked as the original problem")
    original_error = max(abs(payloads[0][key] - inclusive_solved[0][key])
                         for key in ("min", "max"))
    if original_error > 1e-7:
        raise AssertionError(
            f"inclusive level zero differs from standalone original: {original_error}")
    if args.cplex_bin:
        reference = subprocess.run(
            [str(args.cplex_bin), str(target), "--single", "--cplex", "--cplex-threads=1"],
            check=True, capture_output=True, text=True)
        cplex_payloads = [json.loads(line[7:]) for line in reference.stdout.splitlines()
                          if line.startswith("RESULT ")]
        if len(cplex_payloads) != 1 or cplex_payloads[0].get("status") != "optimal":
            raise AssertionError("CPLEX did not solve the converted sparse-row instance")
        error = max(abs(payloads[0][key] - cplex_payloads[0][key])
                    for key in ("min", "max"))
        if error > 5e-5:
            raise AssertionError(f"HiGHS/CPLEX sparse-row disagreement: {error}")
        cplex_path = subprocess.run(
            [str(args.cplex_bin), str(target), "5", "0", "--cplex",
             "--cplex-threads=1"],
            check=True, capture_output=True, text=True)
        cplex_path_payloads = [json.loads(line[7:])
                               for line in cplex_path.stdout.splitlines()
                               if line.startswith("RESULT ")]
        if len(cplex_path_payloads) != len(path_payloads):
            raise AssertionError("CPLEX sparse-row contraction path has wrong length")
        path_error = max(
            abs(odba[key] - cplex[key])
            for odba, cplex in zip(path_payloads, cplex_path_payloads)
            for key in ("min", "max"))
        if path_error > 5e-5:
            raise AssertionError(
                f"HiGHS/CPLEX sparse-row contraction disagreement: {path_error}")
        cplex_inclusive = subprocess.run(
            [str(args.cplex_bin), str(target), "5", "0", "--cplex",
             "--cplex-threads=1", "--contraction-path=inclusive"],
            check=True, capture_output=True, text=True)
        cplex_inclusive_payloads = [json.loads(line[7:])
                                    for line in cplex_inclusive.stdout.splitlines()
                                    if line.startswith("RESULT ")]
        if len(cplex_inclusive_payloads) != len(inclusive_payloads):
            raise AssertionError("CPLEX inclusive contraction path has wrong length")
        inclusive_error = max(
            abs(odba[key] - cplex[key])
            for odba, cplex in zip(inclusive_payloads, cplex_inclusive_payloads)
            for key in ("min", "max"))
        if inclusive_error > 5e-5:
            raise AssertionError(
                f"HiGHS/CPLEX inclusive contraction disagreement: {inclusive_error}")

    # More than 20 free probability coordinates takes the linear projected-
    # midpoint focal path instead of exponential vertex integration.
    large_source = root / "large_focal.m"
    large_lines = ["23 2", "0 0"]
    large_lines.extend(
        f"{index}: [0 1] [0 1]" for index in range(1, 26))
    large_source.write_text("\n".join(large_lines) + "\n", encoding="utf-8")
    subprocess.run([sys.executable, str(args.converter), str(large_source),
                    "--output-root", str(output)],
                   check=True, capture_output=True, text=True)
    focal = subprocess.run(
        [str(args.bin), str(output / "large_focal.ddt"), "1", "0",
         "--contraction-path=interior"],
        check=True, capture_output=True, text=True)
    focal_payloads = [json.loads(line[7:]) for line in focal.stdout.splitlines()
                      if line.startswith("RESULT ")]
    if (len(focal_payloads) != 1 or
            not focal_payloads[0].get("direct_point")):
        raise AssertionError("large focal-point projection fallback failed")

print("Group-B converter and sparse DDT rows passed")
