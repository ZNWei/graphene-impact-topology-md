#!/usr/bin/env python3
"""Generate, run, quality-control, and analyse the isolated P2-F campaign.

P2-F is a new confirmatory experiment.  Its primary endpoint is the terminal
mass fraction of origin-labelled projectile atoms downstream of a fixed plane.
It intentionally does not invoke the legacy intact-body separation classifier.
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


def protocol(path: Path) -> dict[str, Any]:
    cfg = read_json(path)
    cells = cfg["factor_cells"]
    observed = {(int(cell["cross_rebo"]), int(cell["cross_lj"])) for cell in cells}
    expected = {(0, 0), (0, 1), (1, 0), (1, 1)}
    if observed != expected or len(cells) != 4:
        raise SystemExit("P2-F requires the complete 2x2 factor design exactly once")
    n = int(cfg["replicates_per_velocity"])
    if n != 10 or len(cfg["impact_offsets_angstrom"]) != n:
        raise SystemExit("P2-F is pre-registered for exactly ten independent offsets per cell and velocity")
    if len({tuple(map(float, offset)) for offset in cfg["impact_offsets_angstrom"]}) != n:
        raise SystemExit("P2-F offsets must be unique")
    for cell in cells:
        velocities = cell["candidate_velocities_mps"]
        if len(velocities) != 2 or velocities != sorted(velocities) or velocities[0] <= 0:
            raise SystemExit("every P2-F cell requires two strictly increasing positive candidate velocities")
    endpoint = cfg["primary_endpoint"]
    if float(endpoint["control_plane_angstrom"]) >= 0 or not (0 < float(endpoint["event_threshold_fraction"]) < 1):
        raise SystemExit("invalid downstream control plane or binary-event fraction")
    overrides = cfg["baseline_overrides"]
    if overrides.get("box_zhalf_angstrom") != 2000.0 or overrides.get("recovery_ps") != 5.0:
        raise SystemExit("P2-F requires the validated 2000 A half-box and 5 ps recovery")
    if not overrides.get("trajectory_enabled") or int(overrides.get("trajectory_every_steps", 0)) != 1000:
        raise SystemExit("P2-F requires a 0.1 ps trajectory cadence at the prescribed timestep")
    return cfg


def completed(directory: Path) -> bool:
    marker = directory / "final_state.txt"
    return marker.exists() and "RUN_COMPLETE" in marker.read_text(encoding="utf-8", errors="replace")


def lost_atoms(directory: Path) -> bool:
    log = directory / "screen.log"
    return log.exists() and "Lost atoms" in log.read_text(encoding="utf-8", errors="replace")


def prepare_directories() -> None:
    for directory in [ROOT / "generated/base", ROOT / "generated/base_inputs", ROOT / "generated/inputs", ROOT / "generated/metadata", ROOT / "runs", ROOT / "results"]:
        directory.mkdir(parents=True, exist_ok=True)


def prepared_campaign(cfg: dict[str, Any]) -> dict[str, Any]:
    campaign = read_json(PARENT / "config/campaign.json")
    campaign["campaign_id"] = cfg["campaign_id"]
    campaign["potentials"] = {key: "../" + value for key, value in campaign["potentials"].items()}
    campaign["baseline"].update(cfg["baseline_overrides"])
    return campaign


def build_specs(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    campaign = prepared_campaign(cfg)
    specs: list[dict[str, Any]] = []
    for cell in cfg["factor_cells"]:
        for velocity in cell["candidate_velocities_mps"]:
            for index, offset in enumerate(cfg["impact_offsets_angstrom"]):
                spec = make_spec(campaign, "p2f_threshold", int(cell["cross_rebo"]), int(cell["cross_lj"]), float(velocity), int(cfg["seed_start"]) + index, list(offset))
                spec["potentials"]["lj_table"] = "../" + spec["potentials"]["lj_table"]
                specs.append(spec)
    return specs


def generate(args: argparse.Namespace) -> None:
    cfg = protocol(args.protocol.resolve())
    prepare_directories()
    specs = build_specs(cfg)
    expected_runs = 4 * 2 * int(cfg["replicates_per_velocity"])
    if len(specs) != expected_runs or len({spec["run_id"] for spec in specs}) != expected_runs:
        raise SystemExit("generated P2-F grid is not exactly 80 unique impact runs")
    bases = {spec["base_id"]: spec for spec in specs}
    for base_id, spec in bases.items():
        path = ROOT / "generated/base_inputs" / f"{base_id}.in"
        path.write_text(base_input(spec, f"generated/base/{base_id}.data"), encoding="utf-8")
    for spec in specs:
        path = ROOT / "generated/inputs" / f"{spec['run_id']}.in"
        path.write_text(impact_input(spec), encoding="utf-8")
        spec["input_file"] = str(path.relative_to(ROOT))
        spec["input_sha256"] = sha256(path)
        spec["base_input_file"] = str((ROOT / "generated/base_inputs" / f"{spec['base_id']}.in").relative_to(ROOT))
        spec["track_columns"] = TRACK_COLUMNS
        write_json(ROOT / "generated/metadata" / f"{spec['run_id']}.json", spec)
        (ROOT / "runs" / spec["run_id"]).mkdir(exist_ok=True)
    potential_hashes = {name: sha256(ROOT / rel) for name, rel in prepared_campaign(cfg)["potentials"].items()}
    manifest = {
        "protocol": cfg,
        "protocol_sha256": sha256(args.protocol.resolve()),
        "potential_sha256": potential_hashes,
        "runs": specs,
        "base_inputs": sorted(bases),
        "design": "4 factor cells x 2 pre-registered velocities x 10 independent seed/offset replicates = 80 new impacts"
    }
    write_json(ROOT / "generated/p2f_manifest.json", manifest)
    print(json.dumps({"generated_base_states": len(bases), "generated_impact_runs": len(specs), "status": "complete"}, indent=2))


def validate_inputs(args: argparse.Namespace) -> None:
    cfg = protocol(args.protocol.resolve())
    manifest_path = ROOT / "generated/p2f_manifest.json"
    if not manifest_path.exists():
        raise SystemExit("P2-F manifest absent: run generate first")
    manifest = read_json(manifest_path)
    expected = {(int(cell["cross_rebo"]), int(cell["cross_lj"]), int(v)) for cell in cfg["factor_cells"] for v in cell["candidate_velocities_mps"]}
    counts = Counter((int(run["cross_rebo"]), int(run["cross_lj"]), int(run["velocity_mps"])) for run in manifest["runs"])
    failures: list[str] = []
    if set(counts) != expected:
        failures.append("manifest is not the exact pre-registered factor-cell/velocity grid")
    if any(count != int(cfg["replicates_per_velocity"]) for count in counts.values()):
        failures.append("at least one factor-cell/velocity combination does not contain exactly ten replicates")
    if len(manifest["runs"]) != 80 or len({run["run_id"] for run in manifest["runs"]}) != 80:
        failures.append("manifest does not contain exactly 80 unique runs")
    if len(manifest["base_inputs"]) != 10:
        failures.append("P2-F must use ten shared base states, one per seed/offset, across all factor cells and velocities")
    for run in manifest["runs"]:
        input_path = ROOT / run["input_file"]
        text = input_path.read_text(encoding="utf-8") if input_path.exists() else ""
        if not input_path.exists() or sha256(input_path) != run["input_sha256"]:
            failures.append(f"input hash mismatch: {run['run_id']}")
        if "trajectory.lammpstrj" not in text or "dump traj all custom 1000" not in text:
            failures.append(f"required origin-labelled trajectory output is absent: {run['run_id']}")
        if "write_data runs/" not in text or "RUN_COMPLETE" not in text:
            failures.append(f"terminal-data or completion marker is absent: {run['run_id']}")
    report = {"expected_impact_runs": 80, "generated_impact_runs": len(manifest["runs"]), "base_states": len(manifest["base_inputs"]), "per_cell_velocity_counts": {str(key): value for key, value in sorted(counts.items())}, "status": "pass" if not failures else "fail", "failures": failures}
    write_json(ROOT / "results/p2f_input_validation.json", report)
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


def lmp(command: str, input_file: Path, log_file: Path) -> None:
    with log_file.open("w", encoding="utf-8") as handle:
        result = subprocess.run(shlex.split(command) + ["-in", str(input_file.relative_to(ROOT))], cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise RuntimeError(f"LAMMPS failed: {log_file}")


def run(args: argparse.Namespace) -> None:
    manifest_path = ROOT / "generated/p2f_manifest.json"
    if not manifest_path.exists():
        raise SystemExit("P2-F manifest absent: run generate and validate-inputs first")
    manifest = read_json(manifest_path)
    bases = {run["base_id"]: run for run in manifest["runs"]}
    for base_id, spec in bases.items():
        data = ROOT / "generated/base" / f"{base_id}.data"
        if not data.exists():
            lmp(args.lmp, ROOT / spec["base_input_file"], ROOT / "generated/base" / f"{base_id}.screen.log")
        if not data.exists():
            raise RuntimeError(f"base state absent after LAMMPS execution: {base_id}")
    pending = [run for run in manifest["runs"] if args.force or not completed(ROOT / "runs" / run["run_id"])]
    def one(spec: dict[str, Any]) -> str:
        directory = ROOT / "runs" / spec["run_id"]
        lmp(args.lmp, ROOT / spec["input_file"], directory / "screen.log")
        if not completed(directory):
            raise RuntimeError(f"RUN_COMPLETE absent: {spec['run_id']} (inspect {directory / 'screen.log'})")
        return spec["run_id"]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(one, spec) for spec in pending]
        for future in concurrent.futures.as_completed(futures):
            print(future.result(), flush=True)
    print(json.dumps({"requested": len(manifest["runs"]), "ran": len(pending), "status": "complete"}, indent=2))


def read_data_atoms(path: Path) -> list[tuple[int, int, float]]:
    """Read id, type, z from LAMMPS write_data Atoms # full without assumptions about other sections."""
    atoms: list[tuple[int, int, float]] = []
    in_atoms = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("Atoms"):
                in_atoms = True
                continue
            if not in_atoms:
                continue
            if line.startswith(("Velocities", "Bonds", "Angles", "Dihedrals", "Impropers")):
                break
            fields = line.split()
            if len(fields) >= 7 and fields[0].lstrip("+-").isdigit():
                atoms.append((int(fields[0]), int(fields[2]), float(fields[6])))
    if not atoms:
        raise ValueError(f"no Atoms # full records in {path}")
    return atoms


def fractions_from_atoms(atoms: Iterable[tuple[int, int, float]], projectile_type: int, planes: list[float]) -> tuple[int, int, dict[str, float]]:
    rows = list(atoms)
    projectile = [z for _ident, atom_type, z in rows if atom_type == projectile_type]
    target = [z for _ident, atom_type, z in rows if atom_type in (1, 2, 3)]
    if not projectile or not target:
        raise ValueError("projectile or target atoms missing from terminal data")
    payload = {str(plane): sum(z <= plane for z in projectile) / len(projectile) for plane in planes}
    payload.update({"target_" + str(plane): sum(z <= plane for z in target) / len(target) for plane in planes})
    return len(projectile), len(target), payload


def trajectory_series(path: Path, timestep_ps: float, projectile_type: int, planes: list[float], expected_projectile: int) -> list[dict[str, Any]]:
    series: list[dict[str, Any]] = []
    for step, columns, rows in parse_dump(path):
        positions = {name: index for index, name in enumerate(columns)}
        if "type" not in positions or "z" not in positions:
            raise ValueError(f"trajectory lacks type/z columns: {path}")
        pz = [row[positions["z"]] for row in rows if int(row[positions["type"]]) == projectile_type]
        if len(pz) != expected_projectile:
            raise ValueError(f"trajectory projectile atom count differs from expected count at step {step}: {path}")
        series.append({"step": step, "time_ps": step * timestep_ps, "projectile_origin_fraction_by_plane": {str(plane): sum(z <= plane for z in pz) / expected_projectile for plane in planes}})
    if not series:
        raise ValueError(f"trajectory contains no readable frames: {path}")
    return series


def wilson_interval(successes: int, trials: int) -> list[float | None]:
    if trials <= 0:
        return [None, None]
    z = 1.959963984540054
    p = successes / trials
    den = 1 + z*z/trials
    center = (p + z*z/(2*trials)) / den
    half = z * math.sqrt(p*(1-p)/trials + z*z/(4*trials*trials)) / den
    return [max(0.0, center-half), min(1.0, center+half)]


def threshold_estimate(group: list[dict[str, Any]]) -> dict[str, Any]:
    """Two-speed logistic estimate; rejects complete separation rather than inventing a finite threshold."""
    by_velocity: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in group:
        by_velocity[int(row["velocity_mps"])].append(row)
    velocities = sorted(by_velocity)
    summary: list[dict[str, Any]] = []
    for velocity in velocities:
        rows = by_velocity[velocity]
        successes = sum(bool(row["primary_event"]) for row in rows)
        summary.append({"velocity_mps": velocity, "n": len(rows), "events": successes, "event_probability": successes/len(rows) if rows else None, "wilson_95_ci": wilson_interval(successes, len(rows))})
    result: dict[str, Any] = {"event_proportions": summary, "estimator": "two-speed binomial logistic model with an increasing velocity slope"}
    if len(velocities) != 2 or any(item["n"] != 10 for item in summary):
        result.update({"status": "not_estimable_incomplete_design", "reason": "two completed ten-replicate speed groups are required"})
        return result
    low, high = summary
    if low["events"] in (0, 10) or high["events"] in (0, 10):
        result.update({"status": "not_estimable_complete_separation", "reason": "at least one event proportion is 0 or 1; report the bracket and add pre-registered adaptive midpoint runs rather than forcing a finite logistic V_T,50."})
        return result
    p1, p2 = float(low["event_probability"]), float(high["event_probability"])
    if p2 <= p1:
        result.update({"status": "not_estimable_nonmonotonic", "reason": "the observed event probability did not increase across the two pre-registered speeds"})
        return result
    x1, x2 = float(low["velocity_mps"]), float(high["velocity_mps"])
    logit1, logit2 = math.log(p1/(1-p1)), math.log(p2/(1-p2))
    slope = (logit2-logit1)/(x2-x1)
    intercept = logit1-slope*x1
    vt50 = -intercept/slope
    weights = [10*p1*(1-p1), 10*p2*(1-p2)]
    s0 = sum(weights); s1 = weights[0]*x1 + weights[1]*x2; s2 = weights[0]*x1*x1 + weights[1]*x2*x2
    determinant = s0*s2-s1*s1
    if determinant <= 0:
        result.update({"status": "not_estimable_singular_information", "reason": "logistic information matrix is singular"})
        return result
    cov00, cov01, cov11 = s2/determinant, -s1/determinant, s0/determinant
    grad0, grad1 = -1/slope, intercept/(slope*slope)
    variance = grad0*grad0*cov00 + 2*grad0*grad1*cov01 + grad1*grad1*cov11
    se = math.sqrt(max(0.0, variance))
    result.update({"status": "estimated", "slope_per_mps": slope, "v_t_50_mps": vt50, "wald_95_ci_mps": [vt50-1.959963984540054*se, vt50+1.959963984540054*se], "warning": "This is a fragmentation-mediated mass-transmission threshold, not a conventional ballistic V50."})
    return result


def analyse(args: argparse.Namespace) -> None:
    cfg = protocol(args.protocol.resolve())
    manifest = read_json(ROOT / "generated/p2f_manifest.json")
    endpoint = cfg["primary_endpoint"]
    planes = [float(endpoint["control_plane_angstrom"])] + [float(value) for value in cfg["secondary_control_planes_angstrom"]]
    expected_projectile = int(cfg["expected_atom_counts"]["projectile_type_4"])
    expected_target = int(cfg["expected_atom_counts"]["target_types_1_to_3"])
    valid: list[dict[str, Any]] = []
    invalid: list[dict[str, Any]] = []
    incomplete: list[dict[str, Any]] = []
    for spec in manifest["runs"]:
        directory = ROOT / "runs" / spec["run_id"]
        required = [directory / "final_state.txt", directory / "final.data", directory / "track.dat", directory / "trajectory.lammpstrj"]
        if not completed(directory) or any(not path.exists() for path in required):
            incomplete.append({"run_id": spec["run_id"], "reason": "missing_RUN_COMPLETE_or_required_output"})
            continue
        if lost_atoms(directory):
            invalid.append({"run_id": spec["run_id"], "reason": "Lost atoms in screen.log"})
            continue
        try:
            nprojectile, ntarget, fractions = fractions_from_atoms(read_data_atoms(directory / "final.data"), int(endpoint["projectile_type"]), planes)
            if nprojectile != expected_projectile or ntarget != expected_target:
                raise ValueError(f"terminal atom-count mismatch: projectile={nprojectile}, target={ntarget}")
            series = trajectory_series(directory / "trajectory.lammpstrj", float(spec["baseline"]["timestep_ps"]), int(endpoint["projectile_type"]), planes, expected_projectile)
            final_track = parse_track(directory / "track.dat")[-1]
        except (ValueError, OSError) as exc:
            invalid.append({"run_id": spec["run_id"], "reason": str(exc)})
            continue
        primary_fraction = fractions[str(float(endpoint["control_plane_angstrom"]))]
        valid.append({
            "run_id": spec["run_id"], "cross_rebo": spec["cross_rebo"], "cross_lj": spec["cross_lj"], "velocity_mps": spec["velocity_mps"], "seed": spec["seed"], "impact_offset_angstrom": spec["impact_offset_angstrom"],
            "terminal_snapshot_definition": endpoint["definition"], "projectile_origin_atom_count": nprojectile, "target_atom_count": ntarget,
            "terminal_projectile_origin_fraction_by_plane": {str(key): fractions[str(key)] for key in planes},
            "terminal_target_origin_fraction_by_plane": {str(key): fractions["target_" + str(key)] for key in planes},
            "primary_fraction": primary_fraction, "primary_event": primary_fraction >= float(endpoint["event_threshold_fraction"]),
            "final_track": final_track, "trajectory_projectile_origin_fraction_by_plane": series,
            "legacy_intact_body_outcome": "not_evaluated_by_design"
        })
    grouped: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for row in valid:
        grouped[(int(row["cross_rebo"]), int(row["cross_lj"]))].append(row)
    estimates = {f"R{rebo}_L{lj}": threshold_estimate(grouped[(rebo, lj)]) for rebo, lj in [(0, 0), (0, 1), (1, 0), (1, 1)]}
    expected_groups = {(int(cell["cross_rebo"]), int(cell["cross_lj"]), int(v)) for cell in cfg["factor_cells"] for v in cell["candidate_velocities_mps"]}
    observed_counts = Counter((int(row["cross_rebo"]), int(row["cross_lj"]), int(row["velocity_mps"])) for row in valid)
    complete_n10 = set(observed_counts) == expected_groups and all(count == 10 for count in observed_counts.values()) and not invalid and not incomplete
    gate = {"status": "pass" if complete_n10 and all(item["status"] == "estimated" for item in estimates.values()) else "not_ready", "scope": "P2-F data-completeness and estimability gate only; not a manuscript, chemical-mechanism, conventional-V50, DOI, or Q1-submission gate.", "conditions": {"all_80_runs_valid_and_complete": complete_n10, "all_factor_cells_have_estimable_v_t_50": all(item["status"] == "estimated" for item in estimates.values())}}
    payload = {"protocol_sha256": sha256(args.protocol.resolve()), "primary_endpoint": endpoint, "valid_results": valid, "invalid": invalid, "incomplete": incomplete, "replicate_counts": {str(key): value for key, value in sorted(observed_counts.items())}, "all_cells_velocity_n10": complete_n10, "threshold_estimates": estimates, "p2f_gate": gate}
    write_json(ROOT / "results/p2f_transmission_results.json", payload)
    write_json(ROOT / "results/p2f_gate.json", gate)
    with (ROOT / "results/p2f_threshold_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["factor_cell", "velocity_mps", "n", "events", "event_probability", "wilson_95_low", "wilson_95_high", "estimator_status", "v_t_50_mps", "v_t_50_wald_95_low_mps", "v_t_50_wald_95_high_mps"])
        for label, estimate in estimates.items():
            for item in estimate["event_proportions"]:
                ci = item["wilson_95_ci"]
                vtci = estimate.get("wald_95_ci_mps", [None, None])
                writer.writerow([label, item["velocity_mps"], item["n"], item["events"], item["event_probability"], ci[0], ci[1], estimate["status"], estimate.get("v_t_50_mps"), vtci[0], vtci[1]])
    print(json.dumps({"valid": len(valid), "invalid": len(invalid), "incomplete": len(incomplete), "all_cells_velocity_n10": complete_n10, "p2f_gate": gate["status"]}, indent=2))


def status(_args: argparse.Namespace) -> None:
    manifest = read_json(ROOT / "generated/p2f_manifest.json")
    states = Counter("complete" if completed(ROOT / "runs" / run["run_id"]) else "pending" for run in manifest["runs"])
    failures = [{"run_id": run["run_id"], "screen_log": str((ROOT / "runs" / run["run_id"] / "screen.log").relative_to(ROOT))} for run in manifest["runs"] if (ROOT / "runs" / run["run_id"] / "screen.log").exists() and not completed(ROOT / "runs" / run["run_id"])]
    print(json.dumps({"total": len(manifest["runs"]), "states": dict(states), "incomplete_runs_with_logs": failures}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    for name, function in [("generate", generate), ("validate-inputs", validate_inputs), ("analyse", analyse), ("status", status)]:
        command = commands.add_parser(name)
        command.add_argument("--protocol", type=Path, default=ROOT / "config/p2f_protocol.json")
        command.set_defaults(function=function)
    command = commands.add_parser("run")
    command.add_argument("--lmp", required=True)
    command.add_argument("--jobs", type=int, default=1)
    command.add_argument("--force", action="store_true", help="Re-run already completed P2-F impact directories; do not use without an explicit rerun decision.")
    command.set_defaults(function=run)
    args = parser.parse_args()
    if getattr(args, "jobs", 1) < 1:
        raise SystemExit("--jobs must be >= 1")
    args.function(args)


if __name__ == "__main__":
    main()
