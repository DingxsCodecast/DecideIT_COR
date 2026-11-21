#!/usr/bin/env python3
"""Paired DELTA interaction paths; one process = ten DBLP queries + direct focus.

Standard-library only. All jobs use this activated interpreter. Raw stdout,
stderr and flushed numerical diagnostics survive timeouts/interruption.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import statistics
import subprocess
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "benchmarks/revision_20260927"
ALLOWED = ("highs_bcm", "highs_bcm_cold", "highs_bcm_focal", "highs_bcm_rlt_cold", "highs_bcm_diagnostic", "cp", "odba_t1", "odba_t3", "odba_t5", "cplex_global")
STOP = threading.Event()
ACTIVE = set()
ACTIVE_LOCK = threading.Lock()


def interrupt(signum, frame):
    STOP.set()
    with ACTIVE_LOCK:
        for pid in ACTIVE:
            try:
                os.killpg(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_digest():
    paths = [ROOT / "makefile"]
    for folder in ("include", "src_cpp", "src"):
        paths.extend(p for p in (ROOT/folder).rglob("*") if p.suffix in (".h", ".cpp", ".c"))
    content = "".join(f"{p.relative_to(ROOT)} {digest(p)}\n" for p in sorted(paths))
    return hashlib.sha256(content.encode()).hexdigest()


def write_json(path, obj):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
    temp.replace(path)


def clean_environment():
    env = {k: v for k, v in os.environ.items() if not k.startswith("DECIDEIT_")}
    env.update({x: "1" for x in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS")})
    return env


def process(command, directory, timeout):
    if STOP.is_set():
        raise InterruptedError("Experiment interrupted; partial logs retained")
    started = time.perf_counter()
    with (directory/"stdout.log").open("w") as out, (directory/"stderr.log").open("w") as err:
        child = subprocess.Popen(command, stdout=out, stderr=err, env=clean_environment(), start_new_session=True)
        with ACTIVE_LOCK:
            ACTIVE.add(child.pid)
            if STOP.is_set():
                os.killpg(child.pid,signal.SIGTERM)
        timed_out = False
        try:
            rc = child.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt):
            timed_out = True
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                rc = child.wait(timeout=2)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                rc = child.wait()
        finally:
            with ACTIVE_LOCK:
                ACTIVE.discard(child.pid)
    return rc, timed_out, time.perf_counter()-started


def parse_lines(path, prefix=""):
    if not path.exists():
        return []
    result = []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(prefix):
            try:
                result.append(json.loads(line[len(prefix):],parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f"Nonfinite JSON: {x}"))))
            except json.JSONDecodeError:
                raise ValueError(f"Malformed record in {path}")
    return result


def command(profile, source, focal, log, args):
    binary = args.cplex_bin if profile == "cplex_global" else args.bin
    cmd = [str(binary.resolve()), str(source.resolve()), "5", "0", "--sequential-directions",
           "--contraction-path=inclusive", "--contraction-order=inward",
           "--contraction-model=fixed-statements", f"--focal-file={focal}", f"--diagnostics={log}"]
    if profile == "cplex_global":
        return cmd + ["--cplex", "--cplex-threads=1", f"--cplex-time-limit={args.timeout}"]
    method = "cp" if profile == "cp" else "odba" if profile.startswith("odba_") else "bcm"
    trigger = int(profile[-1]) if profile.startswith("odba_") else 3
    cmd += [f"--method={method}", f"--trigger={trigger}", "--climb=guarded", "--lp-backend=highs",
            "--rlt-backend=highs", "--row-filter=cheap", "--contraction-start=transport",
            "--contraction-reduction=refresh", "--contraction-precheck=solve", "--rlt-start=staged",
            "--rlt-simplex-scale-strategy=auto", "--rlt-scale4-max-columns=2500", "--cut-limit=5000"]
    if profile in ("highs_bcm_cold","highs_bcm_rlt_cold"):
        cmd += ["--cold-rlt"]
    if profile in ("highs_bcm_cold","highs_bcm_focal"):
        cmd += ["--contraction-start=focal"]
    if profile=="highs_bcm_diagnostic":
        snapshots=log.parent/"snapshots";snapshots.mkdir(exist_ok=True)
        cmd += ["--rlt-diagnostics", "--rlt-snapshot-dir="+str(snapshots)]
    return cmd


def check_chain(diag, results, focus, original_hash):
    problems = []
    if any(not all(isinstance(x.get(s),(int,float)) and math.isfinite(x[s]) for s in ("min","max","wall_s")) for x in results):
        return ["nonfinite_or_missing_result"]
    models = [x for x in diag if x["event"] == "level_model"]
    directions = [x for x in diag if x["event"] == "direction"]
    if len(results) != 6 or {x["level"] for x in results} != set(range(6)) or len(models) != 5 or len(directions) != 10:
        problems.append("not_exactly_ten_dblp_queries_plus_direct_point")
    for m in models:
        if m["rows_sha256"] != original_hash or m["focal_sha256"] != focus["focal_sha256"]:
            problems.append("changed_statements_or_focal")
        f = m["fraction"]
        for block in ("p", "v"):
            base = models[0]
            for j, centre in enumerate(focus[block]):
                for side in ("lower", "upper"):
                    expected = (1-f)*base[block+"_"+side][j]+f*centre
                    if abs(m[block+"_"+side][j]-expected)>1e-12:
                        problems.append("incorrect_bound_contraction")
    for d in directions:
        if d["witness_violation"] is None or d["witness_violation"]>1e-8 or d["objective_residual"] is None or d["objective_residual"]>1e-6:
            problems.append("invalid_witness_or_objective")
    ordered = sorted(results, key=lambda x:x["level"])
    for previous, current in zip(ordered, ordered[1:]):
        if current["min"] < previous["min"]-5e-5 or current["max"]>previous["max"]+5e-5:
            problems.append("nonnested_endpoint_alert")
    end = [x for x in diag if x["event"] == "chain_complete"]
    if len(end)!=1 or end[0]["dblp_queries"]!=10 or not results[-1].get("direct_point"):
        problems.append("missing_direct_or_chain_marker")
    return sorted(set(problems))


def job(index, item, repetition, args):
    if STOP.is_set():
        raise InterruptedError("Experiment interrupted")
    identity = hashlib.sha256(item["case_id"].encode()).hexdigest()[:16]
    directory = args.output / "runs" / identity / f"rep_{repetition:02d}"
    directory.mkdir(parents=True, exist_ok=True)
    source = ROOT / item["path"]
    if digest(source) != item["sha256"]:
        raise ValueError(f"Input checksum changed: {source}")
    focal = directory / "focus.txt"
    preparation = directory / "preparation.json"
    if not preparation.exists():
        prep_dir = directory / "prepare"
        prep_dir.mkdir()  # partial writes are never silently overwritten
        cmd = [str(args.bin.resolve()),str(source.resolve()),"--sequential-directions",
               "--contraction-model=fixed-statements",f"--prepare-focal={focal}"]
        rc, timeout, wall = process(cmd, prep_dir, args.timeout)
        ds = parse_lines(prep_dir/"stdout.log", "DIAGNOSTIC ")
        if rc or timeout:
            raise ValueError(f"Feasible focus preparation failed: {directory}")
        focus = next(x for x in ds if x["event"]=="focal")
        inp = next(x for x in ds if x["event"]=="input")
        write_json(preparation, {"focus":focus,"input":inp,"process_wall_s":wall,"artifact_sha256":digest(focal)})
    prepared = json.loads(preparation.read_text())
    if digest(focal) != prepared["artifact_sha256"]:
        raise ValueError("Shared focus artifact checksum changed")
    focus = prepared["focus"]
    original_hash = prepared["input"]["original_rows_sha256"]
    profiles = args.profiles
    rotation = (index+repetition-1)%len(profiles)
    order = profiles[rotation:]+profiles[:rotation]
    records = []
    for profile in order:
        if STOP.is_set():
            raise InterruptedError("Experiment interrupted")
        work = directory / profile
        record_path = work / "record.json"
        if record_path.exists():
            records.append(json.loads(record_path.read_text()))
            continue
        work.mkdir()  # refuse ambiguous partial run; use a new run ID
        cmd = command(profile,source,focal,work/"diagnostics.jsonl",args)
        rc, timeout, wall = process(cmd,work,args.timeout)
        diag = parse_lines(work/"diagnostics.jsonl")
        results = parse_lines(work/"stdout.log", "RESULT ")
        errors = [x for x in diag if x["event"]=="error"]
        record = {"case_id":item["case_id"],"input_sha256":item["sha256"],"dimension_per_block":item["dimension"],
            "category":item.get("category"),"profile":profile,"repetition":repetition,"command":cmd,
            "returncode":rc,"process_timeout":timeout,"path_process_wall_s":wall,
            "shared_preparation_process_s":prepared["process_wall_s"],
            "interactive_path_wall_s":wall+prepared["process_wall_s"],
            "focal_preparation_s":focus["preparation_s"],"focal_sha256":focus["focal_sha256"],
            "rows_sha256":original_hash,"diagnostics":str((work/"diagnostics.jsonl").relative_to(args.output)),
            "results":results,"last_error":errors[-1] if errors else None}
        record["reuse_counters"]={name:sum(x.get(name,0) for x in results) for name in
            ("warm_solves","cold_solves","cross_level_basis_uses","cross_level_primal_uses")}
        record["transport_seed_layers"]=sum(x.get("previous_available",False) and
            x.get("contraction_start")=="transport" and x.get("start_min_feasible",False) and
            x.get("start_max_feasible",False) for x in results)
        if timeout:
            record["status"]="timeout"
            record["validation_problems"]=[]
        elif rc:
            controlled = errors and "not certified" in errors[-1]["message"].lower()
            record["status"]="uncertified" if controlled else "execution_error"
            record["validation_problems"]=[]
        else:
            problems = check_chain(diag,results,focus,original_hash)
            record["validation_problems"]=problems
            record["status"]="certified" if not problems else "invalid_certificate"
        write_json(record_path,record)
        records.append(record)
        print(json.dumps({"case_id":item["case_id"],"profile":profile,"rep":repetition,"status":record["status"],
                          "interactive_path_s":record["interactive_path_wall_s"]}),flush=True)
    pairs = []
    reference = next((r for r in records if r["profile"]=="cplex_global"),None)
    if reference:
        for r in records:
            if r is reference or r["status"]!="certified" or reference["status"]!="certified":
                continue
            error = max(abs(a[s]-b[s]) for a,b in zip(r["results"],reference["results"]) for s in ("min","max"))
            pairs.append({"case_id":item["case_id"],"repetition":repetition,"profile":r["profile"],
                          "endpoint_max_error":error,"agreement":error<=5e-5,"precision_alert":error>1e-6,
                          "ratio":r["interactive_path_wall_s"]/reference["interactive_path_wall_s"]})
    write_json(directory/"pairs.json",pairs)
    return records,pairs


def main():
    signal.signal(signal.SIGTERM,interrupt)
    signal.signal(signal.SIGINT,interrupt)
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest",type=Path,default=FROZEN/"local_gate_manifest.json")
    parser.add_argument("--bin",type=Path,default=ROOT/"bin/DELTA_solver")
    parser.add_argument("--cplex-bin",type=Path,default=ROOT/"bin/DELTA_solver_cplex")
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--profiles",nargs="+",choices=ALLOWED,default=["highs_bcm"])
    parser.add_argument("--repetitions",type=int,default=1)
    parser.add_argument("--timeout",type=float,default=1800)
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--local-gate",action="store_true")
    args = parser.parse_args()
    if args.repetitions<1 or args.workers<1 or args.workers>10 or args.timeout<=0 or len(set(args.profiles))!=len(args.profiles):
        parser.error("invalid protocol settings")
    args.output=args.output.resolve()
    manifest=json.loads(args.manifest.read_text())
    items=manifest["instances"]
    if args.local_gate and (len(items)!=59 or args.repetitions!=1 or args.profiles!=["highs_bcm","cplex_global"]):
        parser.error("local gate requires frozen59, one paired repetition and prescribed profiles")
    config={"schema":"decideit-path-execution-v1","manifest_sha256":digest(args.manifest),"solver_source_sha256":source_digest(),
        "bin_sha256":digest(args.bin),
        "cplex_bin_sha256":digest(args.cplex_bin) if "cplex_global" in args.profiles else None,
        "profiles":args.profiles,
        "repetitions":args.repetitions,"timeout":args.timeout,"workers":args.workers,"contraction_model":"fixed-statements",
        "primary_unit":"ten-query interactive path; shared preparation charged once per profile",
        "python":sys.executable,"runner_sha256":digest(Path(__file__))}
    args.output.mkdir(parents=True,exist_ok=True)
    config_path=args.output/"config.json"
    if config_path.exists():
        if json.loads(config_path.read_text())!=config:
            raise ValueError("Resume configuration or source/binary changed; use new run ID")
    else:
        write_json(config_path,config)
    tasks=[(index,item,rep,args) for rep in range(1,args.repetitions+1) for index,item in enumerate(items)]
    records=[];pairs=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for rs,ps in pool.map(lambda task:job(*task),tasks):
            records.extend(rs);pairs.extend(ps)
    statuses={name:sum(r["status"]==name for r in records) for name in sorted({r["status"] for r in records})}
    disagreements=[p for p in pairs if not p["agreement"]]
    # Explicitly unresolved/timeout results are retained but do not become
    # false successes. An execution/numerical error or paired disagreement
    # blocks release; historical failure cases need not all become certified.
    unexpected=[]
    if args.local_gate:
        allowed=set()
        for filename in ("independent_unresolved_manifest.json","tail_disagreement_manifest.json"):
            allowed.update(x["case_id"] for x in json.loads((FROZEN/filename).read_text())["instances"])
        unexpected=[r["case_id"] for r in records if r["status"]=="uncertified" and
                    (r["profile"]!="highs_bcm" or r["case_id"] not in allowed)]
    passed=not disagreements and not unexpected and not any(r["status"] in ("execution_error","invalid_certificate") for r in records)
    pending=sum(r["status"]=="timeout" for r in records)
    summary={"instances":len(items),"records":len(records),"statuses":statuses,"paired_complete":len(pairs),
        "disagreements":disagreements,"precision_alerts":sum(p["precision_alert"] for p in pairs),
        "unexpected_uncertified":unexpected,
        "unresolved_retained":sum(r["status"]=="uncertified" for r in records),"timeouts_pending_review":pending,
        "median_path_ratio":statistics.median(p["ratio"] for p in pairs) if pairs else None,
        "validation_passed":passed,"release_gate_passed":bool(args.local_gate and passed and pending==0 and len(pairs)>0),
        "solver_source_sha256":config["solver_source_sha256"],"manifest_sha256":config["manifest_sha256"]}
    def quantile(values,p):
        if not values:return None
        values=sorted(values);pos=p*(len(values)-1);lo=int(pos);hi=min(lo+1,len(values)-1)
        return values[lo]+(pos-lo)*(values[hi]-values[lo])
    def path_stats(rs):
        success=[r for r in rs if r["status"]=="certified"]
        latency=[r["interactive_path_wall_s"] for r in success]
        return {"attempted_paths":len(rs),"certified_paths":len(success),
            "conditional_success_latency_s":{str(p):quantile(latency,p) for p in (0.5,0.9,0.95,0.99,1.0)},
            "certified_within_deadline_rate_all_attempts":{str(t):sum(v<=t for v in latency)/len(rs) for t in (1,2,5,10)},
            "uncertified":sum(r["status"]=="uncertified" for r in rs),"timeouts":sum(r["status"]=="timeout" for r in rs)}
    summary["interactive_paths_by_profile"]={profile:path_stats([r for r in records if r["profile"]==profile]) for profile in args.profiles}
    summary["interactive_paths_by_dimension"]={str(d):{profile:path_stats([r for r in records if r["profile"]==profile and r["dimension_per_block"]==d])
        for profile in args.profiles} for d in sorted({r["dimension_per_block"] for r in records})}
    write_json(args.output/"summary.json",summary)
    write_json(args.output/"records.json",records)
    write_json(args.output/"pairs.json",pairs)
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
    if not passed or (args.local_gate and not summary["release_gate_passed"]):
        return 2
    return 0


if __name__=="__main__":
    try:
        raise SystemExit(main())
    except (ValueError,OSError,StopIteration) as exc:
        print(f"revision path error: {exc}",file=sys.stderr)
        raise SystemExit(2)
