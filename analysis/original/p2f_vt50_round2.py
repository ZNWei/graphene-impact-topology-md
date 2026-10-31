#!/usr/bin/env python3
"""Generate, run, QC, and pool the preregistered 30-run P2-F round-2 stage.

The first two P2-F stages are read-only snapshots.  This program never
rewrites them and never creates any new R1L0 trajectory.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
Q1 = ROOT.parent
sys.path.insert(0, str(Q1 / "p2f_vt50_resolution_pipeline" / "scripts"))
sys.path.insert(0, str(Q1 / "scripts"))
import p2f_vt50_resolution as common  # noqa: E402
from generate_campaign import TRACK_COLUMNS, base_input, impact_input, make_spec  # noqa: E402
from analyze_campaign import parse_track  # noqa: E402


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def protocol(path: Path) -> dict[str, Any]:
    cfg = read_json(path)
    cells = cfg["factor_cells"]
    expected = {(0, 0, 20500), (0, 1, 21500), (1, 1, 21500)}
    observed = {(int(c["cross_rebo"]), int(c["cross_lj"]), int(c["resolution_velocity_mps"])) for c in cells}
    if observed != expected or len(cells) != 3:
        raise SystemExit("round 2 must be exactly R0L0@20500, R0L1@21500, R1L1@21500")
    if int(cfg["replicates_per_cell"]) != 10 or len(cfg["impact_offsets_angstrom"]) != 10:
        raise SystemExit("round 2 requires exactly ten new independent replicates per included cell")
    if len({tuple(map(float, x)) for x in cfg["impact_offsets_angstrom"]}) != 10:
        raise SystemExit("round-2 impact offsets must be unique")
    if cfg["excluded_estimable_cell"].get("cross_rebo") != 1 or cfg["excluded_estimable_cell"].get("cross_lj") != 0:
        raise SystemExit("R1L0 must remain excluded from round 2")
    b = cfg["baseline_overrides"]
    if b.get("box_zhalf_angstrom") != 2000.0 or b.get("recovery_ps") != 5.0 or not b.get("trajectory_enabled") or int(b.get("trajectory_every_steps", 0)) != 1000:
        raise SystemExit("round 2 must preserve P2-F box, recovery, and trajectory cadence")
    ep = cfg["primary_endpoint"]
    if ep["projectile_type"] != 4 or float(ep["control_plane_angstrom"]) != -500.0 or float(ep["event_threshold_fraction"]) != 0.5:
        raise SystemExit("round 2 must use the unchanged P2-F endpoint")
    return cfg


def directories() -> None:
    for name in ["generated/base", "generated/base_inputs", "generated/inputs", "generated/metadata", "runs", "results"]:
        (ROOT / name).mkdir(parents=True, exist_ok=True)


def campaign(cfg: dict[str, Any]) -> dict[str, Any]:
    c = read_json(Q1 / "config/campaign.json")
    c["campaign_id"] = cfg["campaign_id"]
    c["potentials"] = {key: "../" + value for key, value in c["potentials"].items()}
    c["baseline"].update(cfg["baseline_overrides"])
    return c


def specs(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    c = campaign(cfg)
    result = []
    for cell in cfg["factor_cells"]:
        for index, offset in enumerate(cfg["impact_offsets_angstrom"]):
            spec = make_spec(c, "p2f_vt50_round2", int(cell["cross_rebo"]), int(cell["cross_lj"]), float(cell["resolution_velocity_mps"]), int(cfg["seed_start"]) + index, list(offset))
            spec["potentials"]["lj_table"] = "../" + spec["potentials"]["lj_table"]
            result.append(spec)
    return result


def generate(args: argparse.Namespace) -> None:
    cfg = protocol(args.protocol.resolve()); directories(); rows = specs(cfg)
    if len(rows) != 30 or len({x["run_id"] for x in rows}) != 30:
        raise SystemExit("round-2 grid must be exactly 3 cells x 10 repeats = 30 unique impacts")
    bases = {x["base_id"]: x for x in rows}
    for base_id, spec in bases.items():
        (ROOT / "generated/base_inputs" / f"{base_id}.in").write_text(base_input(spec, f"generated/base/{base_id}.data"), encoding="utf-8")
    for spec in rows:
        input_path = ROOT / "generated/inputs" / f"{spec['run_id']}.in"
        input_path.write_text(impact_input(spec), encoding="utf-8")
        spec["input_file"] = str(input_path.relative_to(ROOT)); spec["input_sha256"] = common.sha256(input_path)
        spec["base_input_file"] = str((ROOT / "generated/base_inputs" / f"{spec['base_id']}.in").relative_to(ROOT)); spec["track_columns"] = TRACK_COLUMNS
        common.write_json(ROOT / "generated/metadata" / f"{spec['run_id']}.json", spec); (ROOT / "runs" / spec["run_id"]).mkdir(exist_ok=True)
    manifest = {"protocol": cfg, "protocol_sha256": common.sha256(args.protocol.resolve()), "potential_sha256": {name: common.sha256(ROOT / value) for name, value in campaign(cfg)["potentials"].items()}, "runs": rows, "base_inputs": sorted(bases), "design": "R0L0@20500 (n=10), R0L1@21500 (n=10), R1L1@21500 (n=10); R1L0 excluded = 30 new impacts"}
    common.write_json(ROOT / "generated/p2f_vt50_round2_manifest.json", manifest)
    print(json.dumps({"generated_base_states": len(bases), "generated_impact_runs": len(rows), "status": "complete"}, indent=2))


def parent_data(cfg: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    failures: list[str] = []; merged: list[dict[str, Any]] = []
    p80 = cfg["parent_snapshots"]["p2f_parent_80_results"]; r1 = cfg["parent_snapshots"]["p2f_round1_40_results"]
    p80_path = ROOT / p80["path"]; r1_path = ROOT / r1["path"]
    if not p80_path.exists() or common.sha256(p80_path) != p80["sha256"]: failures.append("parent 80-run snapshot hash mismatch")
    if not r1_path.exists() or common.sha256(r1_path) != r1["sha256"]: failures.append("round-1 40-run snapshot hash mismatch")
    if failures: return merged, failures
    p80_data = read_json(p80_path); r1_data = read_json(r1_path)
    if len(p80_data.get("valid_results", [])) != int(p80["expected_valid_runs"]) or p80_data.get("invalid") or p80_data.get("incomplete"):
        failures.append("parent 80-run snapshot quality mismatch")
    if len(r1_data.get("new_valid_results", [])) != int(r1["expected_valid_runs"]) or r1_data.get("new_invalid") or r1_data.get("new_incomplete"):
        failures.append("round-1 40-run snapshot quality mismatch")
    if not failures: merged = list(p80_data["valid_results"]) + list(r1_data["new_valid_results"])
    return merged, failures


def validate_inputs(args: argparse.Namespace) -> None:
    cfg = protocol(args.protocol.resolve()); manifest_path = ROOT / "generated/p2f_vt50_round2_manifest.json"
    if not manifest_path.exists(): raise SystemExit("round-2 manifest absent: run generate first")
    manifest = read_json(manifest_path); failures: list[str] = []
    expected = {(int(c["cross_rebo"]), int(c["cross_lj"]), int(c["resolution_velocity_mps"])) for c in cfg["factor_cells"]}
    counts = Counter((int(x["cross_rebo"]), int(x["cross_lj"]), int(x["velocity_mps"])) for x in manifest["runs"])
    if set(counts) != expected or any(x != 10 for x in counts.values()) or len(manifest["runs"]) != 30: failures.append("manifest is not the exact registered 3 x 1 x 10 design")
    if len(manifest["base_inputs"]) != 10: failures.append("round 2 must have ten shared base states")
    for run in manifest["runs"]:
        path = ROOT / run["input_file"]; text = path.read_text(encoding="utf-8") if path.exists() else ""
        if not path.exists() or common.sha256(path) != run["input_sha256"]: failures.append(f"input hash mismatch: {run['run_id']}")
        if "dump traj all custom 1000" not in text or "write_data runs/" not in text or "RUN_COMPLETE" not in text: failures.append(f"required output missing: {run['run_id']}")
    parent, parent_failures = parent_data(cfg); failures.extend(parent_failures)
    if parent and len(parent) != 120: failures.append("pooled parent snapshots do not contain exactly 120 accepted trajectories")
    report = {"expected_new_impact_runs": 30, "generated_new_impact_runs": len(manifest["runs"]), "per_cell_velocity_counts": {str(k):v for k,v in sorted(counts.items())}, "parent_pooled_runs": len(parent), "status": "pass" if not failures else "fail", "failures": failures}
    common.write_json(ROOT / "results/p2f_vt50_round2_input_validation.json", report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


def lmp(command: str, input_file: Path, log_file: Path) -> None:
    with log_file.open("w", encoding="utf-8") as handle:
        result = subprocess.run(shlex.split(command) + ["-in", str(input_file.relative_to(ROOT))], cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if result.returncode: raise RuntimeError(f"LAMMPS failed: {log_file}")


def run(args: argparse.Namespace) -> None:
    manifest = read_json(ROOT / "generated/p2f_vt50_round2_manifest.json"); bases = {x["base_id"]:x for x in manifest["runs"]}
    for base_id, spec in bases.items():
        data = ROOT / "generated/base" / f"{base_id}.data"
        if not data.exists(): lmp(args.lmp, ROOT / spec["base_input_file"], ROOT / "generated/base" / f"{base_id}.screen.log")
        if not data.exists(): raise RuntimeError(f"base state absent after LAMMPS execution: {base_id}")
    pending = [x for x in manifest["runs"] if args.force or not common.completed(ROOT / "runs" / x["run_id"])]
    def one(spec: dict[str, Any]) -> str:
        folder = ROOT / "runs" / spec["run_id"]; lmp(args.lmp, ROOT / spec["input_file"], folder / "screen.log")
        if not common.completed(folder): raise RuntimeError(f"RUN_COMPLETE absent: {spec['run_id']} (inspect {folder / 'screen.log'})")
        return spec["run_id"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in concurrent.futures.as_completed([pool.submit(one, x) for x in pending]): print(future.result(), flush=True)
    print(json.dumps({"requested":len(manifest["runs"]),"ran":len(pending),"status":"complete"}, indent=2))


def grouped_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows: groups[int(row["velocity_mps"])].append(row)
    out = []
    for velocity in sorted(groups):
        g = groups[velocity]; n = len(g); events = sum(bool(x["primary_event"]) for x in g)
        out.append({"velocity_mps":velocity,"n":n,"events":events,"event_probability":events/n,"wilson_95_ci":common.wilson(events,n)})
    return out


def logistic_fit(summary: list[dict[str, Any]], cfg: dict[str, Any]) -> dict[str, Any]:
    """Finite ordinary grouped-binomial logistic fit for three or more speed groups."""
    model = cfg["analysis_model"]; center=float(model["velocity_center_mps"]); scale=float(model["velocity_scale_mps"])
    if len(summary) < 3 or any(int(x["n"]) != 10 for x in summary): return {"status":"not_estimable_incomplete_design","reason":"at least three complete ten-replicate velocity groups are required"}
    probs=[float(x["event_probability"]) for x in summary]
    if any(b<a for a,b in zip(probs,probs[1:])): return {"status":"not_estimable_nonmonotonic","reason":"observed event fractions decrease with velocity"}
    if not (probs[0]<0.5<probs[-1]): return {"status":"not_estimable_unbracketed","reason":"0.5 event probability is not bracketed"}
    mixed=[x for x in summary if 0<int(x["events"])<int(x["n"])]
    if len(mixed)<2: return {"status":"not_estimable_monotone_likelihood","reason":"fewer than two distinct non-extreme velocity groups leave the ordinary grouped-logistic likelihood unbounded"}
    x=[(float(v["velocity_mps"])-center)/scale for v in summary]; y=[int(v["events"]) for v in summary]; n=[int(v["n"]) for v in summary]
    b0,b1=0.0,1.0
    for _ in range(int(model["maximum_iterations"])):
        eta=[b0+b1*xi for xi in x]; p=[1/(1+math.exp(-max(-35,min(35,e)))) for e in eta]; w=[ni*pi*(1-pi) for ni,pi in zip(n,p)]
        g0=sum(yi-ni*pi for yi,ni,pi in zip(y,n,p)); g1=sum(xi*(yi-ni*pi) for xi,yi,ni,pi in zip(x,y,n,p)); h00=sum(w); h01=sum(wi*xi for wi,xi in zip(w,x)); h11=sum(wi*xi*xi for wi,xi in zip(w,x)); det=h00*h11-h01*h01
        if det<=1e-14: return {"status":"not_estimable_singular_information","reason":"logistic information matrix is singular"}
        d0=(h11*g0-h01*g1)/det; d1=(-h01*g0+h00*g1)/det; b0+=d0; b1+=d1
        if max(abs(d0),abs(d1))<float(model["convergence_tolerance"]): break
    else: return {"status":"not_estimable_nonconvergent","reason":"finite grouped-logistic fit did not converge"}
    if b1<=0: return {"status":"not_estimable_nonpositive_slope","reason":"fitted slope is not positive"}
    eta=[b0+b1*xi for xi in x]; p=[1/(1+math.exp(-max(-35,min(35,e)))) for e in eta]; w=[ni*pi*(1-pi) for ni,pi in zip(n,p)]; h00=sum(w); h01=sum(wi*xi for wi,xi in zip(w,x)); h11=sum(wi*xi*xi for wi,xi in zip(w,x)); det=h00*h11-h01*h01
    if det<=1e-14: return {"status":"not_estimable_singular_information","reason":"final information matrix is singular"}
    vt=center-scale*b0/b1; low=min(x["velocity_mps"] for x in summary); high=max(x["velocity_mps"] for x in summary)
    if not low<=vt<=high: return {"status":"not_estimable_threshold_outside_range","reason":"fitted V_T,50 lies outside tested range"}
    cov00,cov01,cov11=h11/det,-h01/det,h00/det; g0=-scale/b1; g1=scale*b0/(b1*b1); se=math.sqrt(max(0,g0*g0*cov00+2*g0*g1*cov01+g1*g1*cov11))
    return {"status":"estimated","estimator":"pooled grouped binomial logistic regression","slope_per_1000_mps":b1,"v_t_50_mps":vt,"wald_95_ci_mps":[vt-1.959963984540054*se,vt+1.959963984540054*se],"warning":"This is a fragmentation-mediated projectile-origin mass-transmission threshold, not conventional ballistic V50."}


def analyse(args: argparse.Namespace) -> None:
    cfg=protocol(args.protocol.resolve()); manifest=read_json(ROOT / "generated/p2f_vt50_round2_manifest.json"); ep=cfg["primary_endpoint"]; planes=[float(ep["control_plane_angstrom"])]+[float(x) for x in cfg["secondary_control_planes_angstrom"]]; expected_p=int(cfg["expected_atom_counts"]["projectile_type_4"]); expected_t=int(cfg["expected_atom_counts"]["target_types_1_to_3"])
    valid=[]; invalid=[]; incomplete=[]
    for spec in manifest["runs"]:
        folder=ROOT / "runs" / spec["run_id"]; needed=[folder/x for x in ["final_state.txt","final.data","track.dat","trajectory.lammpstrj"]]
        if not common.completed(folder) or any(not x.exists() for x in needed): incomplete.append({"run_id":spec["run_id"],"reason":"missing_RUN_COMPLETE_or_required_output"}); continue
        if common.lost_atoms(folder): invalid.append({"run_id":spec["run_id"],"reason":"Lost atoms in screen.log"}); continue
        try:
            np,nt,fractions=common.terminal_fractions(common.read_data_atoms(folder/"final.data"),int(ep["projectile_type"]),planes)
            if (np,nt)!=(expected_p,expected_t): raise ValueError(f"terminal atom-count mismatch: projectile={np}, target={nt}")
            frames=common.trajectory_check(folder/"trajectory.lammpstrj",int(ep["projectile_type"]),expected_p); final_track=parse_track(folder/"track.dat")[-1]
        except (ValueError,OSError) as exc: invalid.append({"run_id":spec["run_id"],"reason":str(exc)}); continue
        f=fractions[str(float(ep["control_plane_angstrom"]))]
        valid.append({"run_id":spec["run_id"],"cross_rebo":spec["cross_rebo"],"cross_lj":spec["cross_lj"],"velocity_mps":spec["velocity_mps"],"seed":spec["seed"],"impact_offset_angstrom":spec["impact_offset_angstrom"],"terminal_snapshot_definition":ep["definition"],"projectile_origin_atom_count":np,"target_atom_count":nt,"terminal_projectile_origin_fraction_by_plane":fractions,"primary_fraction":f,"primary_event":f>=float(ep["event_threshold_fraction"]),"trajectory_frame_count":frames,"final_track":final_track,"legacy_intact_body_outcome":"not_evaluated_by_design"})
    parent,parent_failures=parent_data(cfg); pooled=parent+valid; by=defaultdict(list)
    for row in pooled: by[(int(row["cross_rebo"]),int(row["cross_lj"]))].append(row)
    estimates={}
    for a,b in [(0,0),(0,1),(1,0),(1,1)]:
        s=grouped_summary(by[(a,b)]); fit=logistic_fit(s,cfg); fit["event_proportions"]=s; estimates[f"R{a}_L{b}"]=fit
    counts=Counter((int(x["cross_rebo"]),int(x["cross_lj"]),int(x["velocity_mps"])) for x in valid); expected={(int(c["cross_rebo"]),int(c["cross_lj"]),int(c["resolution_velocity_mps"])) for c in cfg["factor_cells"]}; new_ok=not parent_failures and len(valid)==30 and not invalid and not incomplete and set(counts)==expected and all(n==10 for n in counts.values())
    all_estimated=all(x["status"]=="estimated" for x in estimates.values()); gate={"status":"pass" if new_ok and all_estimated else "not_ready","scope":"P2-F pooled V_T,50 estimability gate only; not a chemical-mechanism, conventional-V50, DOI, or Q1-submission gate.","conditions":{"locked_parent_120_run_snapshots_pass":not parent_failures,"all_30_round2_runs_valid_and_complete":new_ok,"all_factor_cells_have_estimable_v_t_50":all_estimated}}
    payload={"protocol_sha256":common.sha256(args.protocol.resolve()),"new_valid_results":valid,"new_invalid":invalid,"new_incomplete":incomplete,"pooled_run_count":len(pooled),"threshold_estimates":estimates,"p2f_vt50_round2_gate":gate}; common.write_json(ROOT/"results/p2f_vt50_round2_results.json",payload); common.write_json(ROOT/"results/p2f_vt50_round2_gate.json",gate)
    with (ROOT/"results/p2f_vt50_round2_summary.csv").open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h);w.writerow(["factor_cell","velocity_mps","n","events","event_probability","wilson_95_low","wilson_95_high","estimator_status","v_t_50_mps","v_t_50_wald_95_low_mps","v_t_50_wald_95_high_mps"])
        for label,fit in estimates.items():
            ci=fit.get("wald_95_ci_mps",[None,None])
            for x in fit["event_proportions"]: w.writerow([label,x["velocity_mps"],x["n"],x["events"],x["event_probability"],x["wilson_95_ci"][0],x["wilson_95_ci"][1],fit["status"],fit.get("v_t_50_mps"),ci[0],ci[1]])
    print(json.dumps({"new_valid":len(valid),"new_invalid":len(invalid),"new_incomplete":len(incomplete),"pooled_runs":len(pooled),"gate":gate["status"]},indent=2))


def status(_args: argparse.Namespace) -> None:
    m=read_json(ROOT/"generated/p2f_vt50_round2_manifest.json"); states=Counter("complete" if common.completed(ROOT/"runs"/x["run_id"]) else "pending" for x in m["runs"]); failed=[{"run_id":x["run_id"],"screen_log":str((ROOT/"runs"/x["run_id"] / "screen.log").relative_to(ROOT))} for x in m["runs"] if (ROOT/"runs"/x["run_id"] / "screen.log").exists() and not common.completed(ROOT/"runs"/x["run_id"])]
    print(json.dumps({"total_new_runs":len(m["runs"]),"states":dict(states),"incomplete_runs_with_logs":failed},indent=2))


def main() -> None:
    parser=argparse.ArgumentParser(); sub=parser.add_subparsers(dest="command",required=True)
    for name,fn in [("generate",generate),("validate-inputs",validate_inputs),("analyse",analyse),("status",status)]:
        p=sub.add_parser(name);p.add_argument("--protocol",type=Path,default=ROOT/"config/p2f_vt50_round2_protocol.json");p.set_defaults(function=fn)
    p=sub.add_parser("run");p.add_argument("--lmp",required=True);p.add_argument("--jobs",type=int,default=1);p.add_argument("--force",action="store_true",help="Re-run completed round-2 impacts only after an explicit discarded-run decision.");p.set_defaults(function=run)
    args=parser.parse_args()
    if getattr(args,"jobs",1)<1: raise SystemExit("--jobs must be >= 1")
    args.function(args)


if __name__=="__main__": main()
