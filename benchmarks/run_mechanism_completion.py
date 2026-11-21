#!/usr/bin/env python3
"""Fixed factorial and separately instrumented RLT study; no outcome filtering."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))
from run_revision_paths import write_json, parse_lines, clean_environment

PROFILES = ['highs_bcm', 'highs_bcm_focal', 'highs_bcm_rlt_cold',
            'highs_bcm_cold', 'cplex_global']
CHILD = None


def interrupt(signum, frame):
    if CHILD is not None and CHILD.poll() is None:
        CHILD.terminate()  # Path coordinator propagates TERM to its process groups.
    raise KeyboardInterrupt


def run(cmd, log, timeout=None):
    global CHILD
    with log.open('w') as stream:
        CHILD = subprocess.Popen([str(x) for x in cmd], stdout=stream,
                                 stderr=subprocess.STDOUT, env=clean_environment())
        try:
            rc = CHILD.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt):
            CHILD.terminate()
            try:
                CHILD.wait(timeout=10)
            except subprocess.TimeoutExpired:
                CHILD.kill(); CHILD.wait()
            raise
        finally:
            CHILD = None
    return rc


def median(values):
    return statistics.median(values) if values else None


def factorial_summary(folder):
    records = json.loads((folder / 'records.json').read_text())
    keys = sorted({r['case_id'] for r in records})
    paired = []
    for case in keys:
        times = {}
        for p in PROFILES:
            rs = [r for r in records if r['case_id'] == case and r['profile'] == p]
            if rs and all(r['status'] == 'certified' for r in rs):
                times[p] = median([r['interactive_path_wall_s'] for r in rs])
        if len(times) == 5:
            a, b, c, d = (times[p] for p in PROFILES[:4])
            paired.append(dict(case_id=case, times_s=times,
                transport_benefit_with_rlt_warm=b/a,
                transport_benefit_with_rlt_cold=d/c,
                rlt_warm_benefit_with_transport=c/a,
                rlt_warm_benefit_with_focal=d/b,
                multiplicative_interaction=d*a/(b*c),
                baseline_over_cplex=a/times['cplex_global']))
    fields = ['transport_benefit_with_rlt_warm', 'transport_benefit_with_rlt_cold',
              'rlt_warm_benefit_with_transport', 'rlt_warm_benefit_with_focal',
              'multiplicative_interaction', 'baseline_over_cplex']
    return dict(attempted_cases=len(keys), common_all_repetitions_certified_cases=len(paired),
                excluded_from_conditional_ratios=len(keys)-len(paired),
                ratio_interpretation='benefit >1 means enabling that component was faster; interaction is descriptive',
                timing_unit='per-case median complete ten-query path, including shared preparation',
                conditional_paired_medians={f: median([r[f] for r in paired]) for f in fields},
                paired_cases=paired,
                all_attempts_summary=json.loads((folder/'summary.json').read_text()))


def diagnostics_summary(folder, replay_bin, out, timeout):
    calls, attempts, fallbacks = [], [], []
    for log in sorted(folder.glob('runs/*/rep_*/highs_bcm_diagnostic/stderr.log')):
        identity = str(log.parent.relative_to(folder))
        record = json.loads((log.parent/'record.json').read_text())
        for prefix, target in [('RLT_AUDIT ', calls), ('RLT_ATTEMPT ', attempts),
                               ('RLT_FALLBACK ', fallbacks)]:
            for i, row in enumerate(parse_lines(log, prefix)):
                target.append(dict(path=identity, case_id=record['case_id'], sequence=i+1, **row))
    write_json(out/'rlt_calls.json', calls)
    write_json(out/'rlt_attempts.json', attempts)
    write_json(out/'rlt_fallbacks.json', fallbacks)
    replays = []
    for directory in sorted(folder.glob('runs/*/rep_*/highs_bcm_diagnostic/snapshots')):
        candidates = sorted(directory.glob('*.basis'), key=lambda p: int(p.stem.split('_')[-1]))
        # Predeclared selection, independent of speed: first/middle/last available basis.
        chosen = sorted({0, len(candidates)//2, len(candidates)-1}) if candidates else []
        for idx in chosen:
            base = candidates[idx].with_suffix('')
            logfile = base.with_suffix('.replay.log')
            if logfile.exists():
                raise ValueError(f'Existing replay log: {logfile}; use a new output directory')
            rc = run([replay_bin, base.with_suffix('.mps'), base.with_suffix('.options'),
                      base.with_suffix('.basis')], logfile, timeout)
            rows = parse_lines(logfile, 'REPLAY ')
            error = max((abs(r['objective']-rows[0]['objective']) for r in rows), default=None)
            tolerance = 1e-7 * max(1, abs(rows[0]['objective'])) if rows else None
            valid = rc == 0 and len(rows) == 6 and error <= tolerance
            cold = [r['wall_s'] for r in rows if not r['warm']]
            warm = [r['wall_s'] for r in rows if r['warm']]
            artifacts = {suffix: hashlib.sha256(base.with_suffix(suffix).read_bytes()).hexdigest()
                         for suffix in ('.mps','.options','.basis')}
            replays.append(dict(snapshot=str(base.relative_to(folder)), returncode=rc,
                basis_source='before solving this exact LP, not its subsequently optimal basis',
                artifact_sha256=artifacts,
                validated=valid, objective_error=error, objective_tolerance=tolerance,
                cold_over_warm=median(cold)/median(warm) if valid else None, records=rows))
    write_json(out/'replay_records.json', replays)
    return dict(attempted_lp_calls=len(attempts), completed_lp_calls=len(calls),
        cold_fallback_attempts=len(fallbacks),
        supplied_valid_basis_calls=sum(r['basis_valid_before'] for r in attempts),
        supplied_primal_calls=sum(r.get('primal_valid_before',False) for r in attempts),
        basic_variable_count=sum(r['basic_variables'] for r in calls),
        basic_at_bound_counts={str(t): sum(r[f'basic_at_bound_1e{t}'] for r in calls) for t in (8,7,6)},
        lp_calls_with_basic_at_bound={str(t): sum(r[f'basic_at_bound_1e{t}']>0 for r in calls) for t in (8,7,6)},
        simplex_iteration_total=sum(r['iterations'] for r in calls),
        max_primal_infeasibility=max((r['max_primal_infeasibility'] for r in calls),default=None),
        max_dual_infeasibility=max((r['max_dual_infeasibility'] for r in calls),default=None),
        replay_snapshots=len(replays), replay_invalid=sum(not r['validated'] for r in replays),
        descriptive_median_snapshot_cold_over_warm=median([r['cold_over_warm'] for r in replays if r['validated']]),
        limitations=['Returned-basis degeneracy is representation- and tolerance-dependent; not zero-step pivot counts.',
                    'Snapshots are capped at the first 32 RLT calls per path; replay uses first/middle/last available basis.',
                    'Diagnostic execution includes serialization overhead and is excluded from whole-path timing comparisons.',
                    'Fresh cold model versus supplied-basis warm model; identical LP and saved solver options.',
                    'Snapshot observations from one path are correlated; descriptive medians are not case-level inference.'])


def main():
    global CHILD
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bin-dir', type=Path, default=ROOT/'bin')
    parser.add_argument('--timeout', type=float, default=1800)
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error('timeout must be positive')
    out, bins = args.output.resolve(), args.bin_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out/'completion_summary.json').exists():
        raise ValueError('Completed output exists; use a new run ID')
    manifest = ROOT/'benchmarks/revision_20260927/mechanism_manifest.json'
    diagnostic_manifest = manifest.parent/'local_gate_manifest.json'
    if args.smoke:
        original = json.loads((manifest.parent/'local_smoke_manifest.json').read_text())
        original['instances'] = [x for x in original['instances'] if x['dimension'] <= 18]
        original['selection'] = 'Pre-existing local-smoke strata, dimension <=18, before measuring outcomes'
        manifest = out/'smoke_manifest.json'; write_json(manifest, original)
        historical = json.loads(diagnostic_manifest.read_text())
        chosen = {'group_b/1/36/20260422_112354568_0123_00p.ddt',
                  'group_b/2/443/20260424_173332982_2656_00p.ddt',
                  'group_b/3/412/20260425_114121370_8820_00p.ddt'}
        original['instances'] = original['instances'] + [x for x in historical['instances'] if x['case_id'] in chosen]
        original['selection'] = 'Nine small independent strata plus three frozen historical cut-iteration controls; diagnostic cohort only'
        diagnostic_manifest = out/'diagnostic_smoke_manifest.json'; write_json(diagnostic_manifest, original)
    for stage, profiles, reps in [('factorial', PROFILES, 1 if args.smoke else 3),
                                  ('diagnostic', ['highs_bcm_diagnostic'], 1)]:
        write_json(out/'progress.json', dict(state='RUNNING', stage=stage, workers=1))
        (out/'CURRENT_STEP').write_text(stage+'\n')
        cmd = [sys.executable, ROOT/'benchmarks/run_revision_paths.py', '--manifest',
               diagnostic_manifest if stage == 'diagnostic' else manifest,
               '--bin', bins/'DELTA_solver', '--cplex-bin', bins/'DELTA_solver_cplex',
               '--output', out/stage, '--profiles', *profiles, '--workers', '1',
               '--repetitions', str(reps), '--timeout', str(args.timeout)]
        if run(cmd, out/f'{stage}.log'):
            raise ValueError(f'{stage} failed; raw logs retained')
        summary = json.loads((out/stage/'summary.json').read_text())
        if not summary['validation_passed']:
            raise ValueError(f'{stage} has an invalid result')
        # Timeouts and explicit non-certification remain censored observations.
        # They do not block a correctly completed experiment or become successes.
    write_json(out/'progress.json', dict(state='RUNNING', stage='fixed_lp_replay', workers=1))
    (out/'CURRENT_STEP').write_text('fixed_lp_replay\n')
    diag = diagnostics_summary(out/'diagnostic', bins/'rlt_replay', out, args.timeout)
    result = dict(schema='decideit-mechanism-completion-v1', smoke=args.smoke,
        factorial=factorial_summary(out/'factorial'), rlt_diagnostics=diag,
        diagnostic_path_outcomes=json.loads((out/'diagnostic/summary.json').read_text()),
        all_stages_completed=True, replay_validation_passed=diag['replay_invalid']==0 and diag['replay_snapshots']>0)
    write_json(out/'completion_summary.json', result)
    if diag['replay_invalid'] or diag['replay_snapshots']==0:
        raise ValueError('Replay validation failed or zero replay samples; inspect replay_records.json')
    write_json(out/'progress.json', dict(state='COMPLETE', stage='suite_complete', workers=1))
    print(json.dumps({k: v for k,v in result.items() if k not in ('factorial','diagnostic_path_outcomes')},indent=2))


if __name__ == '__main__':
    signal.signal(signal.SIGTERM, interrupt); signal.signal(signal.SIGINT, interrupt)
    try:
        main()
    except KeyboardInterrupt:
        raise SystemExit(130)
    except (ValueError, OSError, subprocess.TimeoutExpired) as exc:
        print(f'mechanism completion error: {exc}', file=sys.stderr)
        raise SystemExit(2)
