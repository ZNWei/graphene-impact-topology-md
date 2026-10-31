#!/usr/bin/env python3
"""Run and analyse pre-registered P2-F 21 km/s threshold-resolution data.

This program creates only the forty new trajectories defined in the protocol.
It never alters the completed P2-F parent data; pooled analysis reads a
SHA-locked parent-results snapshot and writes new results in this package.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PARENT = ROOT.parent
sys.path.insert(0, str(PARENT / "scripts"))
from analyze_campaign import parse_dump, parse_track  # noqa: E402
from generate_campaign import TRACK_COLUMNS, base_input, impact_input, make_spec  # noqa: E402


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def cfg(path: Path) -> dict[str, Any]:
    value = read_json(path)
    cells = value["factor_cells"]
    expected = {(0, 0), (0, 1), (1, 0), (1, 1)}
    observed = {(int(x["cross_rebo"]), int(x["cross_lj"])) for x in cells}
    if observed != expected or len(cells) != 4:
        raise SystemExit("resolution protocol must specify each 2x2 factor cell exactly once")
    if any(int(x["resolution_velocity_mps"]) != 21000 for x in cells):
        raise SystemExit("the registered round-1 resolution velocity is 21000 m/s for every cell")
    if int(value["replicates_per_cell"]) != 10 or len(value["impact_offsets_angstrom"]) != 10:
        raise SystemExit("the resolution stage must contain exactly ten independent repeats per cell")
    if len({tuple(map(float, x)) for x in value["impact_offsets_angstrom"]}) != 10:
        raise SystemExit("resolution impact offsets must be unique")
    b = value["baseline_overrides"]
    if b.get("box_zhalf_angstrom") != 2000.0 or b.get("recovery_ps") != 5.0:
        raise SystemExit("resolution stage must preserve the validated 2000 A half-box and 5 ps recovery")
    if not b.get("trajectory_enabled") or int(b.get("trajectory_every_steps", 0)) != 1000:
        raise SystemExit("resolution stage requires 0.1 ps origin-labelled trajectory output")
    endpoint = value["primary_endpoint"]
    if endpoint["projectile_type"] != 4 or float(endpoint["control_plane_angstrom"]) != -500.0 or float(endpoint["event_threshold_fraction"]) != 0.5:
        raise SystemExit("resolution endpoint must exactly match completed P2-F")
    return value


def completed(directory: Path) -> bool:
    marker = directory / "final_state.txt"
    return marker.exists() and "RUN_COMPLETE" in marker.read_text(encoding="utf-8", errors="replace")


def lost_atoms(directory: Path) -> bool:
    log = directory / "screen.log"
    return log.exists() and "Lost atoms" in log.read_text(encoding="utf-8", errors="replace")


def prepare_directories() -> None:
    for name in ["generated/base", "generated/base_inputs", "generated/inputs", "generated/metadata", "runs", "results"]:
        (ROOT / name).mkdir(parents=True, exist_ok=True)


def prepared_campaign(protocol: dict[str, Any]) -> dict[str, Any]:
    campaign = read_json(PARENT / "config/campaign.json")
    campaign["campaign_id"] = protocol["campaign_id"]
    campaign["potentials"] = {key: "../" + value for key, value in campaign["potentials"].items()}
    campaign["baseline"].update(protocol["baseline_overrides"])
    return campaign


def build_specs(protocol: dict[str, Any]) -> list[dict[str, Any]]:
    campaign = prepared_campaign(protocol)
    specs: list[dict[str, Any]] = []
    for cell in protocol["factor_cells"]:
        for index, offset in enumerate(protocol["impact_offsets_angstrom"]):
            spec = make_spec(campaign, "p2f_vt50_resolution", int(cell["cross_rebo"]), int(cell["cross_lj"]), float(cell["resolution_velocity_mps"]), int(protocol["seed_start"]) + index, list(offset))
            spec["potentials"]["lj_table"] = "../" + spec["potentials"]["lj_table"]
            specs.append(spec)
    return specs


def generate(args: argparse.Namespace) -> None:
    protocol = cfg(args.protocol.resolve())
    prepare_directories()
    specs = build_specs(protocol)
    if len(specs) != 40 or len({x["run_id"] for x in specs}) != 40:
        raise SystemExit("resolution grid must be exactly 4 cells x 10 repeats = 40 unique impacts")
    base_specs = {x["base_id"]: x for x in specs}
    for base_id, spec in base_specs.items():
        (ROOT / "generated/base_inputs" / f"{base_id}.in").write_text(base_input(spec, f"generated/base/{base_id}.data"), encoding="utf-8")
    for spec in specs:
        input_file = ROOT / "generated/inputs" / f"{spec['run_id']}.in"
        input_file.write_text(impact_input(spec), encoding="utf-8")
        spec["input_file"] = str(input_file.relative_to(ROOT))
        spec["input_sha256"] = sha256(input_file)
        spec["base_input_file"] = str((ROOT / "generated/base_inputs" / f"{spec['base_id']}.in").relative_to(ROOT))
        spec["track_columns"] = TRACK_COLUMNS
        write_json(ROOT / "generated/metadata" / f"{spec['run_id']}.json", spec)
        (ROOT / "runs" / spec["run_id"]).mkdir(exist_ok=True)
    manifest = {
        "protocol": protocol,
        "protocol_sha256": sha256(args.protocol.resolve()),
        "potential_sha256": {name: sha256(ROOT / relative) for name, relative in prepared_campaign(protocol)["potentials"].items()},
        "runs": specs,
        "base_inputs": sorted(base_specs),
        "design": "4 factor cells x 1 preregistered 21000 m/s resolution velocity x 10 independent replicates = 40 new impacts"
    }
    write_json(ROOT / "generated/p2f_vt50_resolution_manifest.json", manifest)
    print(json.dumps({"generated_base_states": len(base_specs), "generated_impact_runs": len(specs), "status": "complete"}, indent=2))


def validate_inputs(args: argparse.Namespace) -> None:
    protocol = cfg(args.protocol.resolve())
    manifest_path = ROOT / "generated/p2f_vt50_resolution_manifest.json"
    if not manifest_path.exists():
        raise SystemExit("resolution manifest absent: run generate first")
    manifest = read_json(manifest_path)
    failures: list[str] = []
    counts = Counter((int(x["cross_rebo"]), int(x["cross_lj"]), int(x["velocity_mps"])) for x in manifest["runs"])
    expected = {(int(x["cross_rebo"]), int(x["cross_lj"]), int(x["resolution_velocity_mps"])) for x in protocol["factor_cells"]}
    if set(counts) != expected or any(n != 10 for n in counts.values()) or len(manifest["runs"]) != 40:
        failures.append("manifest is not the exact registered 4 x 1 x 10 design")
    if len(manifest["base_inputs"]) != 10:
        failures.append("resolution stage must contain ten shared base states")
    for spec in manifest["runs"]:
        path = ROOT / spec["input_file"]
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        if not path.exists() or sha256(path) != spec["input_sha256"]:
            failures.append(f"input hash mismatch: {spec['run_id']}")
        if "dump traj all custom 1000" not in text or "trajectory.lammpstrj" not in text:
            failures.append(f"trajectory output absent: {spec['run_id']}")
        if "write_data runs/" not in text or "RUN_COMPLETE" not in text:
            failures.append(f"terminal-output marker absent: {spec['run_id']}")
    parent_file = ROOT / protocol["parent_snapshot"]["result_file"]
    if not parent_file.exists() or sha256(parent_file) != protocol["parent_snapshot"]["result_sha256"]:
        failures.append("parent P2-F result snapshot is absent or its hash differs from the preregistered source")
    parent_manifest = ROOT / "parent_snapshot/p2f_manifest_original.json"
    parent_protocol = ROOT / "parent_snapshot/p2f_protocol_original.json"
    if not parent_manifest.exists() or sha256(parent_manifest) != protocol["parent_snapshot"]["manifest_sha256"]:
        failures.append("parent P2-F manifest snapshot is absent or its hash differs from the preregistered source")
    if not parent_protocol.exists() or sha256(parent_protocol) != protocol["parent_snapshot"]["protocol_sha256"]:
        failures.append("parent P2-F protocol snapshot is absent or its hash differs from the preregistered source")
    parent = read_json(parent_file) if parent_file.exists() else {}
    if len(parent.get("valid_results", [])) != int(protocol["parent_snapshot"]["accepted_parent_run_count"]) or parent.get("invalid") or parent.get("incomplete"):
        failures.append("parent P2-F snapshot is not the accepted 80-valid-run record")
    report = {"expected_new_impact_runs": 40, "generated_new_impact_runs": len(manifest["runs"]), "per_cell_velocity_counts": {str(k): v for k, v in sorted(counts.items())}, "status": "pass" if not failures else "fail", "failures": failures}
    write_json(ROOT / "results/p2f_vt50_resolution_input_validation.json", report)
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


def lmp(command: str, input_file: Path, log_file: Path) -> None:
    with log_file.open("w", encoding="utf-8") as handle:
        result = subprocess.run(shlex.split(command) + ["-in", str(input_file.relative_to(ROOT))], cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise RuntimeError(f"LAMMPS failed: {log_file}")


def run(args: argparse.Namespace) -> None:
    manifest = read_json(ROOT / "generated/p2f_vt50_resolution_manifest.json")
    bases = {x["base_id"]: x for x in manifest["runs"]}
    for base_id, spec in bases.items():
        data = ROOT / "generated/base" / f"{base_id}.data"
        if not data.exists():
            lmp(args.lmp, ROOT / spec["base_input_file"], ROOT / "generated/base" / f"{base_id}.screen.log")
        if not data.exists():
            raise RuntimeError(f"base state absent after LAMMPS execution: {base_id}")
    pending = [x for x in manifest["runs"] if args.force or not completed(ROOT / "runs" / x["run_id"])]
    def one(spec: dict[str, Any]) -> str:
        folder = ROOT / "runs" / spec["run_id"]
        lmp(args.lmp, ROOT / spec["input_file"], folder / "screen.log")
        if not completed(folder):
            raise RuntimeError(f"RUN_COMPLETE absent: {spec['run_id']} (inspect {folder / 'screen.log'})")
        return spec["run_id"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in concurrent.futures.as_completed([pool.submit(one, x) for x in pending]):
            print(future.result(), flush=True)
    print(json.dumps({"requested": len(manifest["runs"]), "ran": len(pending), "status": "complete"}, indent=2))


def read_data_atoms(path: Path) -> list[tuple[int, int, float]]:
    values: list[tuple[int, int, float]] = []
    in_atoms = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("Atoms"):
            in_atoms = True; continue
        if not in_atoms:
            continue
        if line.startswith(("Velocities", "Bonds", "Angles", "Dihedrals", "Impropers")):
            break
        fields = line.split()
        if len(fields) >= 7 and fields[0].lstrip("+-").isdigit():
            values.append((int(fields[0]), int(fields[2]), float(fields[6])))
    if not values:
        raise ValueError(f"no Atoms # full records in {path}")
    return values


def terminal_fractions(atoms: Iterable[tuple[int, int, float]], projectile_type: int, planes: list[float]) -> tuple[int, int, dict[str, float]]:
    atom_rows = list(atoms)
    projectile = [z for _i, kind, z in atom_rows if kind == projectile_type]
    target = [z for _i, kind, z in atom_rows if kind in (1, 2, 3)]
    if not projectile or not target:
        raise ValueError("projectile or target atoms missing from terminal data")
    payload = {str(p): sum(z <= p for z in projectile) / len(projectile) for p in planes}
    return len(projectile), len(target), payload


def trajectory_check(path: Path, projectile_type: int, expected: int) -> int:
    nframes = 0
    for _step, columns, rows in parse_dump(path):
        positions = {name: idx for idx, name in enumerate(columns)}
        if "type" not in positions or "z" not in positions:
            raise ValueError(f"trajectory lacks type/z columns: {path}")
        if sum(int(row[positions["type"]]) == projectile_type for row in rows) != expected:
            raise ValueError(f"trajectory projectile atom count differs from expected count: {path}")
        nframes += 1
    if nframes == 0:
        raise ValueError(f"trajectory contains no readable frames: {path}")
    return nframes


def wilson(successes: int, trials: int) -> list[float]:
    z = 1.959963984540054
    p = successes / trials
    den = 1 + z*z/trials
    mid = (p + z*z/(2*trials)) / den
    half = z * math.sqrt(p*(1-p)/trials + z*z/(4*trials*trials)) / den
    return [max(0.0, mid-half), min(1.0, mid+half)]


def grouped_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[int(row["velocity_mps"])].append(row)
    result = []
    for velocity in sorted(groups):
        group = groups[velocity]
        n = len(group); events = sum(bool(x["primary_event"]) for x in group)
        result.append({"velocity_mps": velocity, "n": n, "events": events, "event_probability": events/n if n else None, "wilson_95_ci": wilson(events, n) if n else [None, None]})
    return result


def logistic_fit(summary: list[dict[str, Any]], protocol: dict[str, Any]) -> dict[str, Any]:
    """Pure-Python grouped-binomial Newton fit, with transparent gate failures."""
    model = protocol["analysis_model"]
    center = float(model["velocity_center_mps"]); scale = float(model["velocity_scale_mps"])
    if len(summary) != 3 or any(x["n"] != 10 for x in summary):
        return {"status": "not_estimable_incomplete_design", "reason": "three complete ten-replicate velocity groups are required"}
    probabilities = [float(x["event_probability"]) for x in summary]
    if any(b < a for a, b in zip(probabilities, probabilities[1:])):
        return {"status": "not_estimable_nonmonotonic", "reason": "observed event fractions decrease with velocity"}
    if not (probabilities[0] < 0.5 < probabilities[-1]):
        return {"status": "not_estimable_unbracketed", "reason": "the 0.5 event probability is not bracketed by the tested velocity span"}
    mixed = [item for item in summary if 0 < int(item["events"]) < int(item["n"])]
    # With only one mixed velocity group, a binomial logistic likelihood is
    # unbounded when all lower groups are 0/n and all upper groups are n/n.
    # Newton iterations then increase the slope indefinitely; that is
    # monotone likelihood caused by separation, not a numerical convergence
    # failure and must not be reported as a finite V_T,50.
    if len(mixed) <= 1:
        return {"status": "not_estimable_monotone_likelihood", "reason": "only one (or no) mixed event-count group remains while lower groups are all failures and upper groups are all events; the ordinary grouped-logistic likelihood is unbounded"}
    x = [(float(item["velocity_mps"]) - center) / scale for item in summary]
    y = [int(item["events"]) for item in summary]
    n = [int(item["n"]) for item in summary]
    beta0, beta1 = 0.0, 1.0
    for _ in range(int(model["maximum_iterations"])):
        eta = [beta0 + beta1 * xi for xi in x]
        p = [1.0/(1.0 + math.exp(-max(-35.0, min(35.0, e)))) for e in eta]
        w = [ni*pi*(1-pi) for ni, pi in zip(n, p)]
        g0 = sum(yi-ni*pi for yi, ni, pi in zip(y, n, p))
        g1 = sum(xi*(yi-ni*pi) for xi, yi, ni, pi in zip(x, y, n, p))
        h00 = sum(w); h01 = sum(wi*xi for wi, xi in zip(w, x)); h11 = sum(wi*xi*xi for wi, xi in zip(w, x))
        det = h00*h11 - h01*h01
        if det <= 1e-14:
            return {"status": "not_estimable_singular_information", "reason": "logistic information matrix is singular"}
        d0 = (h11*g0-h01*g1)/det; d1 = (-h01*g0+h00*g1)/det
        beta0 += d0; beta1 += d1
        if max(abs(d0), abs(d1)) < float(model["convergence_tolerance"]):
            break
    else:
        return {"status": "not_estimable_nonconvergent", "reason": "grouped logistic fit did not converge"}
    if beta1 <= 0:
        return {"status": "not_estimable_nonpositive_slope", "reason": "fitted logistic slope is not positive"}
    eta = [beta0 + beta1 * xi for xi in x]
    p = [1.0/(1.0 + math.exp(-max(-35.0, min(35.0, e)))) for e in eta]
    w = [ni*pi*(1-pi) for ni, pi in zip(n, p)]
    h00 = sum(w); h01 = sum(wi*xi for wi, xi in zip(w, x)); h11 = sum(wi*xi*xi for wi, xi in zip(w, x)); det = h00*h11-h01*h01
    if det <= 1e-14:
        return {"status": "not_estimable_singular_information", "reason": "final logistic information matrix is singular"}
    cov00, cov01, cov11 = h11/det, -h01/det, h00/det
    vt50 = center - scale*beta0/beta1
    grad0 = -scale/beta1; grad1 = scale*beta0/(beta1*beta1)
    variance = grad0*grad0*cov00 + 2*grad0*grad1*cov01 + grad1*grad1*cov11
    lo, hi = min(x["velocity_mps"] for x in summary), max(x["velocity_mps"] for x in summary)
    if not (lo <= vt50 <= hi):
        return {"status": "not_estimable_threshold_outside_range", "reason": "fitted V_T,50 lies outside the observed velocity span"}
    se = math.sqrt(max(0.0, variance))
    return {"status": "estimated", "estimator": "three-speed grouped binomial logistic regression", "slope_per_1000_mps": beta1, "v_t_50_mps": vt50, "wald_95_ci_mps": [vt50-1.959963984540054*se, vt50+1.959963984540054*se], "warning": "This is a fragmentation-mediated projectile-origin mass-transmission threshold, not conventional ballistic V50."}


def analyse(args: argparse.Namespace) -> None:
    protocol = cfg(args.protocol.resolve())
    manifest = read_json(ROOT / "generated/p2f_vt50_resolution_manifest.json")
    endpoint = protocol["primary_endpoint"]; planes = [float(endpoint["control_plane_angstrom"])] + [float(x) for x in protocol["secondary_control_planes_angstrom"]]
    expected_p = int(protocol["expected_atom_counts"]["projectile_type_4"]); expected_t = int(protocol["expected_atom_counts"]["target_types_1_to_3"])
    valid: list[dict[str, Any]] = []; invalid: list[dict[str, Any]] = []; incomplete: list[dict[str, Any]] = []
    for spec in manifest["runs"]:
        folder = ROOT / "runs" / spec["run_id"]
        required = [folder / x for x in ["final_state.txt", "final.data", "track.dat", "trajectory.lammpstrj"]]
        if not completed(folder) or any(not x.exists() for x in required):
            incomplete.append({"run_id": spec["run_id"], "reason": "missing_RUN_COMPLETE_or_required_output"}); continue
        if lost_atoms(folder):
            invalid.append({"run_id": spec["run_id"], "reason": "Lost atoms in screen.log"}); continue
        try:
            np, nt, fractions = terminal_fractions(read_data_atoms(folder / "final.data"), int(endpoint["projectile_type"]), planes)
            if (np, nt) != (expected_p, expected_t):
                raise ValueError(f"terminal atom-count mismatch: projectile={np}, target={nt}")
            frames = trajectory_check(folder / "trajectory.lammpstrj", int(endpoint["projectile_type"]), expected_p)
            final_track = parse_track(folder / "track.dat")[-1]
        except (ValueError, OSError) as exc:
            invalid.append({"run_id": spec["run_id"], "reason": str(exc)}); continue
        value = fractions[str(float(endpoint["control_plane_angstrom"]))]
        valid.append({"run_id": spec["run_id"], "cross_rebo": spec["cross_rebo"], "cross_lj": spec["cross_lj"], "velocity_mps": spec["velocity_mps"], "seed": spec["seed"], "impact_offset_angstrom": spec["impact_offset_angstrom"], "terminal_snapshot_definition": endpoint["definition"], "projectile_origin_atom_count": np, "target_atom_count": nt, "terminal_projectile_origin_fraction_by_plane": fractions, "primary_fraction": value, "primary_event": value >= float(endpoint["event_threshold_fraction"]), "trajectory_frame_count": frames, "final_track": final_track, "legacy_intact_body_outcome": "not_evaluated_by_design"})
    parent_path = ROOT / protocol["parent_snapshot"]["result_file"]
    parent = read_json(parent_path) if parent_path.exists() else {}
    parent_ok = parent_path.exists() and sha256(parent_path) == protocol["parent_snapshot"]["result_sha256"] and len(parent.get("valid_results", [])) == 80 and not parent.get("invalid") and not parent.get("incomplete")
    pooled = list(parent.get("valid_results", [])) + valid
    by_cell: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for row in pooled: by_cell[(int(row["cross_rebo"]), int(row["cross_lj"]))].append(row)
    estimates: dict[str, Any] = {}
    new_by_cell = Counter((int(x["cross_rebo"]), int(x["cross_lj"])) for x in valid)
    for rebo, lj in [(0,0),(0,1),(1,0),(1,1)]:
        label = f"R{rebo}_L{lj}"; summary = grouped_summary(by_cell[(rebo,lj)])
        new_summary = [x for x in summary if int(x["velocity_mps"]) == 21000]
        fit = logistic_fit(summary, protocol)
        if len(new_summary) != 1 or new_summary[0]["events"] in (0,10):
            fit["resolution_stage_diagnostic"] = "new_21000_mps_group_is_extreme"
            fit["resolution_stage_diagnostic_reason"] = "The new 21000 m/s group has 0/10 or 10/10 events. This is recorded as a design diagnostic; the ordinary grouped-logistic estimator status above remains the authoritative estimability result."
        fit["event_proportions"] = summary
        estimates[label] = fit
    all_new_valid = len(valid) == 40 and not invalid and not incomplete and all(new_by_cell[(a,b)] == 10 for a,b in [(0,0),(0,1),(1,0),(1,1)])
    all_estimated = all(x["status"] == "estimated" for x in estimates.values())
    gate = {"status": "pass" if parent_ok and all_new_valid and all_estimated else "not_ready", "scope": "P2-F pooled V_T,50 estimability gate only; not a chemical-mechanism, conventional-V50, DOI, or Q1-submission gate.", "conditions": {"parent_80_run_snapshot_hash_and_quality_pass": parent_ok, "all_40_new_resolution_runs_valid_and_complete": all_new_valid, "all_factor_cells_have_estimable_v_t_50": all_estimated}}
    payload = {"protocol_sha256": sha256(args.protocol.resolve()), "parent_snapshot_sha256": sha256(parent_path) if parent_path.exists() else None, "new_valid_results": valid, "new_invalid": invalid, "new_incomplete": incomplete, "pooled_parent_plus_resolution_run_count": len(pooled), "threshold_estimates": estimates, "p2f_vt50_resolution_gate": gate}
    write_json(ROOT / "results/p2f_vt50_resolution_results.json", payload); write_json(ROOT / "results/p2f_vt50_resolution_gate.json", gate)
    with (ROOT / "results/p2f_vt50_resolution_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle); writer.writerow(["factor_cell","velocity_mps","n","events","event_probability","wilson_95_low","wilson_95_high","estimator_status","v_t_50_mps","v_t_50_wald_95_low_mps","v_t_50_wald_95_high_mps"])
        for label, fit in estimates.items():
            ci = fit.get("wald_95_ci_mps", [None,None])
            for item in fit["event_proportions"]:
                writer.writerow([label,item["velocity_mps"],item["n"],item["events"],item["event_probability"],item["wilson_95_ci"][0],item["wilson_95_ci"][1],fit["status"],fit.get("v_t_50_mps"),ci[0],ci[1]])
    print(json.dumps({"new_valid":len(valid),"new_invalid":len(invalid),"new_incomplete":len(incomplete),"pooled_runs":len(pooled),"gate":gate["status"]}, indent=2))


def status(_args: argparse.Namespace) -> None:
    manifest = read_json(ROOT / "generated/p2f_vt50_resolution_manifest.json")
    states = Counter("complete" if completed(ROOT / "runs" / x["run_id"]) else "pending" for x in manifest["runs"])
    failures = [{"run_id": x["run_id"], "screen_log": str((ROOT / "runs" / x["run_id"] / "screen.log").relative_to(ROOT))} for x in manifest["runs"] if (ROOT / "runs" / x["run_id"] / "screen.log").exists() and not completed(ROOT / "runs" / x["run_id"])]
    print(json.dumps({"total_new_runs":len(manifest["runs"]),"states":dict(states),"incomplete_runs_with_logs":failures}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(); subs = parser.add_subparsers(dest="command", required=True)
    for name, func in [("generate",generate),("validate-inputs",validate_inputs),("analyse",analyse),("status",status)]:
        command = subs.add_parser(name); command.add_argument("--protocol",type=Path,default=ROOT / "config/p2f_vt50_resolution_protocol.json"); command.set_defaults(function=func)
    command = subs.add_parser("run"); command.add_argument("--lmp",required=True); command.add_argument("--jobs",type=int,default=1); command.add_argument("--force",action="store_true",help="Re-run completed impacts only after an explicit discarded-run decision."); command.set_defaults(function=run)
    args = parser.parse_args()
    if getattr(args,"jobs",1) < 1: raise SystemExit("--jobs must be >= 1")
    args.function(args)


if __name__ == "__main__":
    main()
