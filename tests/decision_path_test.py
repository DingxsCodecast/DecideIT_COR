"""Feasible-focus, shared-path and failure-record regression (no timing claims)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def diagnostics(process):
    return [json.loads(x[11:]) for x in process.stdout.splitlines()
            if x.startswith("DIAGNOSTIC ")]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, required=True)
    parser.add_argument("--cplex-bin", type=Path)
    args = parser.parse_args()
    env = {k: v for k, v in os.environ.items() if not k.startswith("DECIDEIT_")}
    with tempfile.TemporaryDirectory(prefix="decideit_focal_") as tmp:
        tmp = Path(tmp)
        text = (ROOT / "tests/smoke.ddt").read_text()
        # Both candidates violate ORIGINAL general rows. V's corrected focus
        # cannot be the midpoint, exercising asymmetric bound contraction.
        text = text.replace("alternative_B\n4\n", "alternative_B\n5\n")
        text = text.replace("values\n4\n", "3 2 1 1 -1 2 1 -1 -1.6\nvalues\n5\n")
        text += "3 2 1 1 -1 1 2 -1 -1.5\n"
        source = tmp / "general.ddt"
        source.write_text(text)
        focal = tmp / "focus.txt"
        def run(binary, extra):
            return subprocess.run([str(binary), str(source), "5", "0",
                "--contraction-model=fixed-statements", "--sequential-directions", *extra],
                capture_output=True, text=True, timeout=60, env=env)
        prepared = run(args.bin, [f"--prepare-focal={focal}"])
        assert prepared.returncode == 0, prepared.stderr
        records = diagnostics(prepared)
        focus = next(x for x in records if x["event"] == "focal")
        assert focus["projection_lps"] == 2 and focus["original_model_violation"] <= 1e-8
        assert next(x for x in records if x["event"] == "input")["input_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
        payload = "\n".join(focal.read_text().splitlines()[:7]) + "\n"
        assert hashlib.sha256(payload.encode()).hexdigest() == focus["focal_sha256"]
        results = []
        for binary in [args.bin] + ([args.cplex_bin] if args.cplex_bin else []):
            log = tmp / (binary.name + ".jsonl")
            extra = [f"--focal-file={focal}", f"--diagnostics={log}"]
            if binary == args.cplex_bin:
                extra += ["--cplex", "--cplex-threads=1"]
            process = run(binary, extra)
            assert process.returncode == 0, (process.stdout, process.stderr)
            diag = diagnostics(process)
            assert diag == [json.loads(x) for x in log.read_text().splitlines()]
            model = [x for x in diag if x["event"] == "level_model"]
            assert len(model) == 5 and len({x["rows_sha256"] for x in model}) == 1
            for item in model:
                fraction = item["fraction"]
                for block in ("p", "v"):
                    for j, f in enumerate(focus[block]):
                        assert abs(item[block + "_lower"][j] - fraction*f) <= 1e-12
                        assert abs(item[block + "_upper"][j] - (1-fraction+fraction*f)) <= 1e-12
            directions = [x for x in diag if x["event"] == "direction"]
            assert len(directions) == 10
            assert all(x["witness_violation"] <= 1e-8 and x["objective_residual"] <= 1e-6 for x in directions)
            assert next(x for x in diag if x["event"] == "chain_complete")["dblp_queries"] == 10
            results.append([json.loads(x[7:]) for x in process.stdout.splitlines() if x.startswith("RESULT ")])
            assert run(binary, extra).returncode != 0  # no log overwrite
        if len(results) == 2:
            assert all(abs(a[s] - b[s]) <= 1e-6 for a, b in zip(*results) for s in ("min", "max"))
        focal.write_text(focal.read_text().replace(focus["focal_sha256"], "0"*64))
        invalid = run(args.bin, [f"--focal-file={focal}"])
        assert invalid.returncode != 0
        assert diagnostics(invalid)[-1]["event"] == "error"
        assert "RESULT " not in invalid.stdout
        failed = subprocess.run([str(args.bin), str(ROOT/"tests/internal_basis.ddt"),
            "--single", "--method=cp", "--sequential-directions"],
            capture_output=True,text=True,timeout=30,env={**env,"DECIDEIT_REVERSE_SEARCH":"tnp"})
        assert failed.returncode != 0
        directions = [x for x in diagnostics(failed) if x["event"]=="direction"]
        assert len(directions)==2 and directions[0]["attempted"] and not directions[1]["attempted"]
        assert directions[0]["termination_flag"]==-9
    print("Decision path: shared feasible focus, unchanged rows, asymmetric contraction, ten queries, witnesses and artifact guards passed")


if __name__ == "__main__":
    main()
