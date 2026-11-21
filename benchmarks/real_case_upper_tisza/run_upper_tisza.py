"""Audit all nine published-data pairwise comparisons on the 10-query path.

Use the chosen DELTA_solver and optionally a locally built CPLEX comparator.
"""

import argparse
import hashlib
import json
import subprocess
import tempfile
from decimal import Decimal
from pathlib import Path

from build_upper_tisza import (
    D, EVENT_PROBABILITIES, NET_INCOME, SCALE, canonical_intervals,
    verify_source_transcription,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOLERANCE_NORMALIZED = 5e-5  # paper's declared external agreement tolerance


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def greedy_probability_extreme(coefficients, maximize=False, bounds=EVENT_PROBABILITIES):
    """Exact LP solution for boxed common p with sum(p)=1."""
    probabilities = [lo for lo, _ in bounds]
    residual = D(1) - sum(probabilities)
    for i in sorted(range(5), key=lambda i: coefficients[i], reverse=maximize):
        increment = min(residual, bounds[i][1] - probabilities[i])
        probabilities[i] += increment
        residual -= increment
    assert abs(residual) < D("1e-13")
    return sum(p * c for p, c in zip(probabilities, coefficients))


def independent_baseline(perspective, first, second):
    a = canonical_intervals(perspective, first)
    b = canonical_intervals(perspective, second)
    min_coef = [a[i][0] - b[i][1] for i in range(5)]
    max_coef = [a[i][1] - b[i][0] for i in range(5)]
    return (greedy_probability_extreme(min_coef),
            greedy_probability_extreme(max_coef, maximize=True))


def run(binary, case, cplex=False):
    command = [str(binary), str(case), "5", "0", "--method=bcm",
               "--contraction-model=fixed-statements", "--sequential-directions"]
    if cplex:
        command += ["--cplex", "--cplex-threads=1"]
    with tempfile.TemporaryDirectory() as tmp:
        log = Path(tmp) / "diagnostics.jsonl"
        completed = subprocess.run(command + [f"--diagnostics={log}"],
                                   text=True, capture_output=True, timeout=120, check=False)
        diagnostics = [json.loads(line) for line in log.read_text().splitlines()] if log.exists() else []
    if completed.returncode:
        raise RuntimeError(f"{case.name}: {binary.name} exit={completed.returncode}\n"
                           f"stdout={completed.stdout[-1500:]}\n"
                           f"stderr={completed.stderr[-1500:]}")
    results = [json.loads(line[7:]) for line in completed.stdout.splitlines()
               if line.startswith("RESULT ")]
    levels = {record["level"]: record for record in results}
    if set(levels) != set(range(6)):
        raise AssertionError(f"{case.name}: expected levels 0..5, got {sorted(levels)}")
    if any(levels[k]["status"] != "optimal" for k in levels):
        raise AssertionError(f"{case.name}: non-optimal result")
    if any(x["min"] > x["max"] + 1e-10 for x in results):
        raise AssertionError(f"{case.name}: reversed endpoints")
    for previous, current in zip(results, results[1:]):
        if current["min"] < previous["min"] - 1e-10 or current["max"] > previous["max"] + 1e-10:
            raise AssertionError(f"{case.name}: contraction intervals are not nested")
    witnesses = [x for x in diagnostics if x["event"] == "direction"]
    if len(witnesses) != 10 or any(x["witness_violation"] > 1e-8 or x["objective_residual"] > 1e-6
                                     for x in witnesses):
        raise AssertionError(f"{case.name}: invalid or missing witnesses")
    if not levels[5].get("direct_point") or abs(levels[5]["max"] - levels[5]["min"]) > 1e-10:
        raise AssertionError(f"{case.name}: 100% contraction was not a direct point")
    completed_event = next(x for x in diagnostics if x.get("event") == "chain_complete")
    if completed_event["dblp_queries"] != 10:
        raise AssertionError(f"{case.name}: not exactly 10 DBLP queries")
    focal = next(x for x in diagnostics if x.get("event") == "focal")
    input_event = next(x for x in diagnostics if x.get("event") == "input")
    if input_event["contraction_model"] != "fixed-statements":
        raise AssertionError(f"{case.name}: wrong contraction model")
    if focal["original_model_violation"] > 1e-8:
        raise AssertionError(f"{case.name}: focal point infeasible")
    return {
        "command": command,
        "levels": results,
        "chain": completed_event,
        "focal_sha256": focal["focal_sha256"],
        "focal_original_model_violation": focal["original_model_violation"],
        "input_sha256": input_event["input_sha256"],
        "focus": {"p": focal["p"], "v": focal["v"]},
    }


def render_figure(records, target):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.0), sharex=True)
    labels = {"government": "Government", "insurers": "Insurers",
              "property_owners": "Property owners"}
    for ax, perspective in zip(axes, labels):
        record = next(r for r in records if r["perspective"] == perspective
                      and r["first"] == 1 and r["second"] == 2)
        levels = record["bcm"]["levels"]
        x = [20 * item["level"] for item in levels]
        lo = [item["min"] * float(SCALE) for item in levels]
        hi = [item["max"] * float(SCALE) for item in levels]
        ax.fill_between(x, lo, hi, color="#2a78b4", alpha=0.22)
        ax.plot(x, lo, color="#144b77", marker="o", markersize=2.8,
                linewidth=1.3, label="Lower")
        ax.plot(x, hi, color="#bc6731", marker="s", markersize=2.8,
                linewidth=1.3, label="Upper")
        ax.axhline(0, color="0.45", linestyle="--", linewidth=0.8)
        ax.set_title(labels[perspective])
        ax.set_xticks([0, 20, 40, 60, 80, 100])
        ax.grid(alpha=0.18)
    axes[0].set_ylabel("Scenario 1 minus 2 (million HUF)")
    for ax in axes:
        ax.set_xlabel("Fixed-statement contraction (%)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(target, bbox_inches="tight", metadata={"Creator": "DecideIT case audit"})
    plt.close(fig)
    return str(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, default=ROOT / "bin/DELTA_solver")
    parser.add_argument("--cplex", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    solver = args.solver.resolve(strict=True)
    reference = args.cplex.resolve(strict=True) if args.cplex else None
    verify_source_transcription()
    records = []
    for perspective in NET_INCOME:
        for first, second in ((1, 2), (1, 3), (2, 3)):
            case = HERE / "cases" / f"upper_tisza_{perspective}_s{first}_s{second}.ddt"
            bcm = run(solver, case)
            cplex = run(reference, case, True) if reference else None
            maximum_error = None
            if cplex:
                if bcm["focal_sha256"] != cplex["focal_sha256"]:
                    raise AssertionError(f"{case.name}: unequal focal point")
                maximum_error = max(abs(a[k] - b[k])
                                    for a, b in zip(bcm["levels"], cplex["levels"])
                                    for k in ("min", "max"))
                if maximum_error > TOLERANCE_NORMALIZED:
                    raise AssertionError(f"{case.name}: CPLEX gap {maximum_error}")
            independent = independent_baseline(perspective, first, second)
            independent_norm = [float(x / SCALE) for x in independent]
            independent_error = max(abs(bcm["levels"][0][k] - x)
                                    for k, x in zip(("min", "max"), independent_norm))
            if independent_error > 1e-9:
                raise AssertionError(f"{case.name}: baseline oracle error {independent_error}")
            all_level_errors = []
            for level in bcm["levels"]:
                lam = D(level["level"]) / 5
                p = [D(str(x)) for x in bcm["focus"]["p"][:5]]
                v = [D(str(x)) * SCALE - D("2450") for x in bcm["focus"]["v"]]
                p_box = [(lo*(1-lam)+lam*x, hi*(1-lam)+lam*x)
                         for (lo, hi), x in zip(EVENT_PROBABILITIES, p)]
                # At the final focal point use the saved values directly. A
                # finite decimal representation can sum to 1 only approximately.
                a_box, b_box = [
                    [(lo*(1-lam)+lam*x, hi*(1-lam)+lam*x)
                     for (lo, hi), x in zip(canonical_intervals(perspective, scenario), v[offset:offset+5])]
                    for scenario, offset in [(first, 0), (second, 5)]]
                lower = [a[0]-b[1] for a, b in zip(a_box, b_box)]
                upper = [a[1]-b[0] for a, b in zip(a_box, b_box)]
                expected = [greedy_probability_extreme(lower, bounds=p_box),
                            greedy_probability_extreme(upper, maximize=True, bounds=p_box)]
                error = max(abs(level[k]-float(x/SCALE)) for k, x in zip(("min", "max"), expected))
                if error > 1e-10:
                    raise AssertionError(f"{case.name} level {level['level']}: independent boxed-LP error {error}")
                all_level_errors.append(error)
            records.append({
                "perspective": perspective, "first": first, "second": second,
                "case": str(case), "case_sha256": sha256(case),
                "baseline_exact_mhuf": [str(x) for x in independent],
                "baseline_oracle_error_normalized": independent_error,
                "maximum_cplex_error_normalized": maximum_error,
                "maximum_all_level_lp_error_normalized": max(all_level_errors),
                "bcm": bcm, "cplex": cplex,
            })
            print(f"PASS {perspective} {first}v{second}: "
                  f"raw=[{independent[0]}, {independent[1]}] MHUF; "
                  f"max boxed-LP error={max(all_level_errors):.3g}", flush=True)
    figure = render_figure(records, args.output / "upper_tisza_s1_s2_contraction.pdf") if args.plot else None
    summary = {
        "source_url": "https://cdn.aaai.org/FLAIRS/2003/Flairs03-096.pdf",
        "source_tables": "1--4", "objective_unit": "million HUF",
        "normalization": {"lower": str(D("-2450")), "upper": str(D("2.96")),
                          "difference_multiplier": str(SCALE)},
        "max_abs_bcm_cplex_error_normalized": max(
            (r["maximum_cplex_error_normalized"] for r in records if reference), default=None),
        "max_abs_all_level_lp_error_normalized": max(r["maximum_all_level_lp_error_normalized"] for r in records),
        "max_abs_bcm_baseline_oracle_error_normalized": max(
            r["baseline_oracle_error_normalized"] for r in records),
        "cases": 9, "nonpoint_queries": 90, "direct_points": 9,
        "bcm_binary_sha256": sha256(solver), "cplex_binary_sha256": sha256(reference) if reference else None,
        "generator_sha256": sha256(HERE / "build_upper_tisza.py"),
        "figure": figure, "records": records,
    }
    (args.output / "upper_tisza_results.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("AUDIT COMPLETE: 9 comparisons, 90 DBLP queries, 9 direct points")


if __name__ == "__main__":
    main()
