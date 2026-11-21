# LP-only dominance queries for interval decision analysis

Research code accompanying *An LP-Only Solver for Repeated Dominance
Queries in Interval Decision Analysis*. The solver evaluates both
endpoints of an expected-value difference over separate linear
probability and value bases. It implements pure cutting planes (CP),
immediate RLT bounding (BCM), and delayed RLT bounding (ODBA).

This release uses HiGHS for every LP. The paper's production timing
profile uses a separately held Pivot-Copy backend for eligible small LPs
and HiGHS for RLT. Its timings are tied to that build. The public
HiGHS configurations use the same decision-query and certificate logic
and have distinct benchmark profile names. The optional CPLEX adapter
independently solves the original nonconvex formulation using a local
CPLEX Studio installation.

## Build

A POSIX environment (Linux or macOS), a C++11 compiler, GNU Make, HiGHS
with its C++ headers and library, and GNU Scientific Library with its
headers and libraries are required. The local verification used HiGHS
1.12.0 on macOS. Python 3.10 or later runs the checks and experiment drivers.
Matplotlib is optional for the application figure.

```sh
make -j4
make test PYTHON=python3
```

On macOS the build uses the Homebrew prefix if available. Set
`HIGHS_ROOT` and `GSL_ROOT` to the installed prefixes when they differ
from the default. `HIGHS_LIBDIR` and `GSL_LIBDIR` can override the
library directories independently. On a Mac with Homebrew HiGHS and GSL
under `/usr/local`, pass `GSL_ROOT=/usr/local` to both commands.

For the separate comparator and the RLT replay tool:

```sh
make cplex CPLEX_ROOT=/path/to/CPLEX_Studio
make rlt-replay
```

CPLEX 22.2 is the documented comparison version. Set `CPLEX_ARCH` or
`CPLEX_LIBDIR` if the installer uses another library layout.

## Verify the paper's examples

```sh
python3 tests/check_illustrative.py
python3 benchmarks/real_case_upper_tisza/run_upper_tisza.py --output results/tisza
```

Both commands accept `--cplex bin/DELTA_solver_cplex` after the optional
build. The second command has a `--plot` option and requires a new output
directory. It checks all nine pairwise comparisons across three perspectives, ten DBLP
queries per pair, their shrinking intervals, and the final direct focal
points against a separate boxed-LP calculation. Its input provenance is
documented in `build_upper_tisza.py`, based on Tables 1 through 4 of the
[original Upper Tisza study](https://cdn.aaai.org/FLAIRS/2003/Flairs03-096.pdf).

The first command checks the full-precision example, evaluates exact
rational witnesses in the original constraints, and checks the ODBA
incumbent and RLT gap. For the lower endpoint at the first RLT test,
the current value is about −0.315555438, the best feasible value is
−0.3872680777, and the remaining-region bound is about −0.331499856.
The residual gap is zero. The upper direction is minimized with its
objective negated internally and is converted back to a positive upper
endpoint on output. Default CPLEX and ODBA values in the paper differ
by about `1.30e-6` at that upper endpoint. Verified output is saved in
`examples/illustrative_audit.json` and
`benchmarks/real_case_upper_tisza/upper_tisza_results.json`.

To inspect the raw trace:

```sh
bin/DELTA_solver examples/illustrative_example.ddt --single \
  --method=odba --trigger=3 --climb=guarded --row-filter=cheap \
  --lp-backend=highs --rlt-backend=highs --rlt-start=staged \
  --sequential-directions --trace
```

## Complete comparison paths

```sh
python3 benchmarks/run_revision_paths.py \
  --manifest benchmarks/revision_20260927/local_smoke_manifest.json \
  --profiles highs_bcm --workers 1 --output results/smoke
```

Use `independent_manifest.json` (210 fixed inputs) or
`stress_manifest.json` (654 screened inputs) to run the corresponding
cohort. Use profiles `cp`, `odba_t1`, `odba_t3`, `odba_t5`,
and `highs_bcm` for the specialized policies. Add `cplex_global`
after building the optional comparator. The runner prepares one shared
feasible focus per case/repetition, holds the original statements fixed,
solves min/max at 0/20/40/60/80 percent and evaluates the 100 percent
point directly. It preserves timeouts and uncertified paths in its
records. Use `--workers 1` and `--repetitions 3` for serial latency
comparisons. The HiGHS-only smoke run on macOS certified 20 of the 21
smoke paths and recorded the one expected unresolved path.

`benchmarks/run_mechanism_completion.py` runs the warm-state factorial
and a separately instrumented RLT study with `make cplex rlt-replay`.
It supports `--smoke` for a short check. `bin/rlt_replay` reports
the linked HiGHS version with no arguments. When passed an MPS model,
HiGHS options file and basis, it replays three paired cold/warm solves.
The archive uses absolute primal/dual LP feasibility checks and
floating-point optimality tests. Inspect the path status and witnesses
before using bounds for close decisions.

## Source and release rights

`src_cpp/` and `include/` contain the solver. `tools/` contains the
DDT converter, instance generator and replay utility. The frozen data and
experimental protocols are in `benchmarks/`. The layout excludes the
copyright-restricted Pivot-Copy implementation and headers, licensed
solver binaries, raw machine logs and previously built executables.

The authors should record the applicable source license and confirm
redistribution rights before making this directory public. The build
uses separately installed HiGHS and GSL, and optionally CPLEX. Existing
copyright and dependency terms apply. Cite the accompanying article
using its final bibliographic details once available.
