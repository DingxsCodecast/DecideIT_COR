"""Check the article's full-precision example, witnesses, and RLT gap."""

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "examples/illustrative_example.ddt"


def model():
    """Read only this example's interval and sparse-row DDT statements."""
    lines = CASE.read_text().splitlines()
    assert lines[4] == "2 11 11"
    blocks = []
    index = 7
    for nrows in (30, 36):
        assert int(lines[index]) == 22 + nrows
        constraints = []
        for line in lines[index + 1:index + 1 + 22 + nrows]:
            kind, *fields = line.split()
            coefficients = [Fraction(0)] * 22
            if kind == "1":
                alternative, variable, sign, lo, hi = fields
                assert sign == "1"
                j = 11 * (int(alternative) - 1) + int(variable) - 1
                coefficients[j] = Fraction(1)
                constraints.append((coefficients, Fraction(hi)))
                constraints.append(([-x for x in coefficients], -Fraction(lo)))
            else:
                assert kind == "3" and (len(fields) - 2) == 3 * int(fields[0])
                for a, j, v in zip(fields[1:-1:3], fields[2:-1:3], fields[3:-1:3]):
                    coefficients[11 * (int(a) - 1) + int(j) - 1] += Fraction(v)
                constraints.append((coefficients, Fraction(fields[-1])))
        blocks.append(constraints)
        index += 1 + 22 + nrows + 1  # inter-block DDT metadata
    p, v = blocks
    for a in range(2):
        row = [Fraction(i // 11 == a) for i in range(22)]
        p.extend([(row, Fraction(1)), ([-x for x in row], Fraction(-1))])
    return p, v


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def rational_vertex(constraints, point):
    """Reconstruct a nearby vertex, then test every inequality exactly."""
    numeric = [Fraction(str(x)) for x in point]
    active = [(a, b) for a, b in constraints if abs(dot(a, numeric) - b) <= Fraction("1e-9")]
    pivots = {}
    for a, b in active:
        row = list(a) + [b]
        for j in range(len(numeric)):
            step = row[j]
            if not step:
                continue
            if j in pivots:
                row = [x - step * y for x, y in zip(row, pivots[j])]
            else:
                pivots[j] = [x / step for x in row]
                break
        else:
            assert row[-1] == 0, "Inconsistent active rows"
    assert len(pivots) == len(numeric), "Witness active rows do not determine a vertex"
    exact = [Fraction(0)] * len(numeric)
    for j in reversed(range(len(numeric))):
        row = pivots[j]
        exact[j] = row[-1] - dot(row[j + 1:-1], exact[j + 1:])
    assert all(dot(a, exact) <= b for a, b in constraints), "Exact witness violates an original row"
    assert max(abs(x - y) for x, y in zip(exact, numeric)) < Fraction("1e-8")
    return exact


def packet(stdout, prefix):
    return [json.loads(line[len(prefix):]) for line in stdout.splitlines()
            if line.startswith(prefix)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, default=ROOT / "bin/DELTA_solver")
    parser.add_argument("--cplex", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as directory:
        diagnostics = Path(directory) / "witnesses.jsonl"
        run = subprocess.run([str(args.solver.resolve(strict=True)), str(CASE),
                              "--single", "--sequential-directions", "--method=odba",
                              "--trigger=3", "--trace", "--lp-backend=highs",
                              "--rlt-backend=highs", "--climb=guarded",
                              "--row-filter=cheap", "--rlt-start=staged",
                              f"--diagnostics={diagnostics}"],
                             text=True, capture_output=True, timeout=120)
        if run.returncode:
            raise RuntimeError(run.stdout[-1800:] + run.stderr[-1800:])
        traces = packet(run.stdout, "TRACE ")
        result, = packet(run.stdout, "RESULT ")
        witnesses = [event for event in map(json.loads, diagnostics.read_text().splitlines())
                     if event["event"] == "direction"]
    assert result["status"] == "optimal" and len(witnesses) == 2
    assert math.isclose(result["min"], -0.3872680777, rel_tol=0, abs_tol=1e-7)
    assert math.isclose(result["max"], 0.6871355880, rel_tol=0, abs_tol=1e-7)
    assert result["min"] < 0 < result["max"]
    p_model, v_model = model()
    findings = {}
    for witness in witnesses:
        sense = witness["sense"]
        p = rational_vertex(p_model, witness["p"])
        v = rational_vertex(v_model, witness["v"])
        objective = dot(p[:11], v[:11]) - dot(p[11:], v[11:])
        assert abs(float(objective) - result[sense]) < 1e-8
        trace_sense = sense if sense == "min" else "negated_max"
        active = [x for x in traces if x.get("sense") == trace_sense
                  and x.get("rlt_active") and "lower_bound" in x]
        last = active[-1]
        observations = [x for x in traces if x.get("sense") == trace_sense
                        and "current_search_value" in x]
        sequence = [x["incumbent"] for x in observations]
        assert sorted({x["cuts"] for x in observations}) == [0, 1, 2, 3]
        assert all(new <= old + 1e-12 for old, new in zip(sequence, sequence[1:]))
        assert last["cuts"] == 3
        sign = 1 if sense == "min" else -1
        assert abs(last["incumbent"] - sign * float(objective)) < 1e-8
        assert last["current_search_value"] >= last["incumbent"] - 1e-12
        assert last["lower_bound"] <= last["current_search_value"] + 1e-8
        assert last["gap"] == max(0, last["incumbent"] - last["lower_bound"]) == 0
        findings[sense] = {"exact_feasible_witness_objective": float(objective),
                           "internal_direction": trace_sense,
                           "trace": last, "observations": observations}
    if args.cplex:
        run = subprocess.run([str(args.cplex.resolve(strict=True)), str(CASE),
                              "--single", "--cplex", "--cplex-threads=1",
                              "--cplex-time-limit=120"],
                             text=True, capture_output=True, check=True, timeout=130)
        reference, = packet(run.stdout, "RESULT ")
        assert reference["status"] == "optimal"
        assert max(abs(reference[s] - result[s]) for s in ("min", "max")) <= 5e-5
        findings["cplex_reference"] = {s: reference[s] for s in ("min", "max")}
    print(json.dumps(findings, indent=2))


if __name__ == "__main__":
    main()
