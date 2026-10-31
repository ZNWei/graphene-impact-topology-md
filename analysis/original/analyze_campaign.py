#!/usr/bin/env python3
"""Analyse P1-P4 generated runs without touching legacy data or scripts."""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

TRACK_COLUMNS = [
    "step", "time_ps", "stage", "ke_all_eV", "pe_all_eV", "etotal_eV",
    "ke_target_eV", "ke_projectile_eV", "ke_projectile_cm_eV",
    "ke_projectile_nontrans_eV", "pe_target_assigned_eV",
    "pe_projectile_assigned_eV", "e_reactive_eV", "e_zbl_eV", "e_kc_eV",
    "e_lj_tail_eV", "target_zmin_A", "target_zmax_A", "projectile_zmin_A",
    "projectile_zmax_A", "projectile_zcm_A", "projectile_vz_mps"
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_track(path: Path) -> list[dict[str, float]]:
    columns: list[str] | None = None
    records: list[dict[str, float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("# step time_ps"):
            columns = line[2:].split()
        elif line.startswith("#") or not line.strip():
            continue
        else:
            if columns is None:
                columns = TRACK_COLUMNS
            values = line.split()
            if len(values) != len(columns):
                continue
            records.append({key: float(value) for key, value in zip(columns, values)})
    if not records:
        raise ValueError(f"no numerical track records in {path}")
    return records


def parse_initial(path: Path) -> dict[str, float]:
    values = path.read_text(encoding="utf-8", errors="replace").strip().split()
    if len(values) != 8 or values[0] != "INITIAL":
        raise ValueError(f"malformed initial state: {path}")
    keys = ["etotal_eV", "ke_projectile_eV", "ke_target_eV", "pe_all_eV", "ke_projectile_cm_eV", "ke_projectile_nontrans_eV"]
    return {key: float(value) for key, value in zip(keys, values[2:])}


def plateau_sd(records: list[dict[str, float]], window_ps: float) -> float:
    end = records[-1]["time_ps"]
    selected = [r["projectile_vz_mps"] for r in records if r["time_ps"] >= end-window_ps]
    return statistics.pstdev(selected) if len(selected) >= 2 else float("inf")


def classify(records: list[dict[str, float]], cfg: dict[str, Any]) -> tuple[str, dict[str, float]]:
    last = records[-1]
    cutoff = cfg["termination_cutoff_angstrom"]
    sd = plateau_sd(records, cfg["plateau_window_ps"])
    separated_bottom = last["projectile_zmax_A"] < last["target_zmin_A"] - cutoff
    separated_top = last["projectile_zmin_A"] > last["target_zmax_A"] + cutoff
    stable = sd <= cfg["plateau_vz_sd_mps"]
    if separated_bottom and last["projectile_vz_mps"] < 0 and stable:
        outcome = "perforated"
    elif separated_top and last["projectile_vz_mps"] > 0 and stable:
        outcome = "rebounded"
    else:
        outcome = "unresolved"
    details = {"plateau_vz_sd_mps": sd, "bottom_separated": float(separated_bottom), "top_separated": float(separated_top), "stable": float(stable)}
    return outcome, details


def parse_dump(path: Path) -> Iterable[tuple[int, list[str], list[list[float]]]]:
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        while True:
            line = fh.readline()
            if not line:
                return
            if line.strip() != "ITEM: TIMESTEP":
                continue
            timestep = int(fh.readline().strip())
            if fh.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError(f"malformed dump {path}")
            n = int(fh.readline().strip())
            if not fh.readline().startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"missing box bounds in {path}")
            for _ in range(3):
                fh.readline()
            atoms_header = fh.readline().strip().split()
            if atoms_header[:2] != ["ITEM:", "ATOMS"]:
                raise ValueError(f"missing atom header in {path}")
            cols = atoms_header[2:]
            rows = []
            for _ in range(n):
                rows.append([float(v) for v in fh.readline().split()])
            yield timestep, cols, rows


def structure_metrics(path: Path, metadata: dict[str, Any]) -> dict[str, Any]:
    """Distance/persistence analysis; never labels a distance contact as a bond."""
    try:
        import numpy as np
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise RuntimeError("trajectory analysis needs numpy and scipy (scipy.spatial.cKDTree)") from exc
    cfg = metadata["baseline"]
    dt_ps = cfg["timestep_ps"]
    cutoffs = cfg["bond_distance_cutoffs_angstrom"]
    contact_cutoff = cutoffs[1]
    impact_center = np.array([cfg["box_xy_angstrom"]/2 + metadata["impact_offset_angstrom"][0], cfg["box_xy_angstrom"]/2 + metadata["impact_offset_angstrom"][1]])
    local_radius = cfg["local_gap_radius_angstrom"]
    edge_width = cfg["edge_width_angstrom"]
    sxlo = (cfg["box_xy_angstrom"] - cfg["sheet_size_angstrom"])/2
    sxhi = (cfg["box_xy_angstrom"] + cfg["sheet_size_angstrom"])/2
    active_pt: dict[tuple[int, int], float] = {}
    active_ll: dict[tuple[int, int], float] = {}
    pt_lifetimes: list[float] = []
    ll_lifetimes: list[float] = []
    pt_counts: list[int] = []
    ll_counts: list[int] = []
    gaps: list[dict[str, float]] = []
    baseline_damage: float | None = None
    final_damage: float | None = None
    first_time: float | None = None
    last_time: float | None = None

    def arrays(cols: list[str], rows: list[list[float]]):
        ix = {name: i for i, name in enumerate(cols)}
        a = np.asarray(rows, dtype=float)
        return a[:, ix["id"]].astype(int), a[:, ix["type"]].astype(int), a[:, [ix["x"], ix["y"], ix["z"]]]

    def target_damage(ids, types, xyz, edge_ids: set[int]) -> int:
        damaged = 0
        for layer in (1, 2, 3):
            mask = types == layer
            layer_ids, coords = ids[mask], xyz[mask]
            tree = cKDTree(coords)
            coord = [len(v)-1 for v in tree.query_ball_point(coords, contact_cutoff)]
            damaged += sum(1 for atom_id, n in zip(layer_ids, coord) if atom_id not in edge_ids and n < 2.5)
        return damaged

    edge_ids: set[int] = set()
    previous_time: float | None = None
    for timestep, cols, rows in parse_dump(path):
        ids, types, xyz = arrays(cols, rows)
        time_ps = timestep * dt_ps
        if first_time is None:
            first_time = time_ps
            target = (types >= 1) & (types <= 3)
            for atom_id, (x, y, _z) in zip(ids[target], xyz[target]):
                if x <= sxlo+edge_width or x >= sxhi-edge_width or y <= sxlo+edge_width or y >= sxhi-edge_width:
                    edge_ids.add(int(atom_id))
        frame_delta = 0.0 if previous_time is None else max(0.0, time_ps-previous_time)
        previous_time = time_ps
        last_time = time_ps
        projectile_mask = types == 4
        target_mask = (types >= 1) & (types <= 3)
        p_ids, p_xyz = ids[projectile_mask], xyz[projectile_mask]
        t_ids, t_xyz = ids[target_mask], xyz[target_mask]
        p_tree, t_tree = cKDTree(p_xyz), cKDTree(t_xyz)
        pair_lists = p_tree.query_ball_tree(t_tree, contact_cutoff)
        pt_pairs = {(int(p_ids[i]), int(t_ids[j])) for i, neighbors in enumerate(pair_lists) for j in neighbors}
        for key in set(active_pt).difference(pt_pairs):
            pt_lifetimes.append(active_pt[key])
            del active_pt[key]
        for key in pt_pairs:
            active_pt[key] = active_pt.get(key, 0.0) + frame_delta
        pt_counts.append(len(pt_pairs))

        ll_pairs: set[tuple[int, int]] = set()
        gap_values: list[float] = []
        for lo, hi in ((1, 2), (2, 3)):
            lo_mask, hi_mask = types == lo, types == hi
            lo_ids, hi_ids = ids[lo_mask], ids[hi_mask]
            lo_xyz, hi_xyz = xyz[lo_mask], xyz[hi_mask]
            lo_tree = cKDTree(lo_xyz)
            hits = cKDTree(hi_xyz).query_ball_tree(lo_tree, contact_cutoff)
            ll_pairs.update((int(hi_ids[i]), int(lo_ids[j])) for i, js in enumerate(hits) for j in js)
            radial = np.linalg.norm(hi_xyz[:, :2]-impact_center, axis=1) <= local_radius
            if radial.any():
                xy_tree = cKDTree(lo_xyz[:, :2])
                _d, nearest = xy_tree.query(hi_xyz[radial, :2], k=1)
                z_gaps = np.abs(hi_xyz[radial, 2] - lo_xyz[np.asarray(nearest), 2])
                gap_values.extend(float(v) for v in z_gaps)
        for key in set(active_ll).difference(ll_pairs):
            ll_lifetimes.append(active_ll[key])
            del active_ll[key]
        for key in ll_pairs:
            active_ll[key] = active_ll.get(key, 0.0) + frame_delta
        ll_counts.append(len(ll_pairs))
        if gap_values:
            gaps.append({"time_ps": time_ps, "median_A": float(np.median(gap_values)), "q05_A": float(np.quantile(gap_values, .05)), "q95_A": float(np.quantile(gap_values, .95)), "n": len(gap_values)})
        current_damage = target_damage(ids, types, xyz, edge_ids)
        if baseline_damage is None:
            baseline_damage = current_damage
        final_damage = current_damage
    pt_lifetimes.extend(active_pt.values())
    ll_lifetimes.extend(active_ll.values())
    if first_time is None or final_damage is None:
        raise ValueError(f"no frames in {path}")
    return {
        "definition": "distance-threshold contacts and persistence, not bond-order evidence",
        "contact_cutoff_A": contact_cutoff,
        "frames": len(pt_counts), "duration_ps": (last_time-first_time) if last_time is not None else 0.0,
        "projectile_target_contact_count_peak": max(pt_counts, default=0),
        "projectile_target_contact_count_final": pt_counts[-1] if pt_counts else 0,
        "projectile_target_contact_lifetime_ps_max": max(pt_lifetimes, default=0.0),
        "projectile_target_contact_lifetime_ps_median": statistics.median(pt_lifetimes) if pt_lifetimes else 0.0,
        "layer_layer_contact_count_peak": max(ll_counts, default=0),
        "layer_layer_contact_count_final": ll_counts[-1] if ll_counts else 0,
        "layer_layer_contact_lifetime_ps_max": max(ll_lifetimes, default=0.0),
        "layer_layer_contact_lifetime_ps_median": statistics.median(ll_lifetimes) if ll_lifetimes else 0.0,
        "target_damage_edge_excluded_baseline": baseline_damage,
        "target_damage_edge_excluded_final": final_damage,
        "target_damage_edge_excluded_net": final_damage-baseline_damage,
        "local_interlayer_gap_time_series": gaps,
    }


def run_result(root: Path, metadata: dict[str, Any]) -> dict[str, Any] | None:
    run_dir = root / "runs" / metadata["run_id"]
    track_path, initial_path = run_dir / "track.dat", run_dir / "initial_state.txt"
    if not (track_path.exists() and initial_path.exists() and (run_dir / "final_state.txt").exists()):
        return None
    track = parse_track(track_path)
    initial = parse_initial(initial_path)
    outcome, terminal = classify(track, metadata["baseline"])
    final = track[-1]
    drift = (final["etotal_eV"]-initial["etotal_eV"])/max(1.0, abs(initial["etotal_eV"]))
    result = {
        "run_id": metadata["run_id"], "stage": metadata["stage"], "cross_rebo": metadata["cross_rebo"], "cross_lj": metadata["cross_lj"],
        "velocity_mps": metadata["velocity_mps"], "seed": metadata["seed"], "impact_offset_angstrom": metadata["impact_offset_angstrom"],
        "reactive_model": metadata["reactive_model"], "variation": metadata["variation"], "base_id": metadata["base_id"],
        "outcome": outcome, "terminal": terminal, "initial": initial,
        "final": final, "energy_drift_fraction": drift,
        "projectile_ke_loss_eV": initial["ke_projectile_eV"]-final["ke_projectile_eV"],
        "projectile_com_ke_loss_eV": initial["ke_projectile_cm_eV"]-final["ke_projectile_cm_eV"],
        "projectile_nontrans_ke_change_eV": final["ke_projectile_nontrans_eV"]-initial["ke_projectile_nontrans_eV"],
        "energy_closure_note": "NVE during impact/recovery; no reflecting z wall; frozen edge velocity is explicitly zero."
    }
    trajectory = run_dir / "trajectory.lammpstrj"
    if trajectory.exists():
        result["structure"] = structure_metrics(trajectory, metadata)
    return result


def mean_ci(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"n": 0, "mean": None, "sd": None, "sem": None}
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return {"n": len(values), "mean": statistics.mean(values), "sd": sd, "sem": sd/math.sqrt(len(values))}


def threshold_plan(results: list[dict[str, Any]], campaign: dict[str, Any]) -> dict[str, Any]:
    cells = []
    screen = [r for r in results if r["stage"] == "screen"]
    for rebo, lj in ((0, 0), (0, 1), (1, 0), (1, 1)):
        rows = sorted((r for r in screen if r["cross_rebo"] == rebo and r["cross_lj"] == lj), key=lambda x: x["velocity_mps"])
        rebound = [r["velocity_mps"] for r in rows if r["outcome"] == "rebounded"]
        perforated = [r["velocity_mps"] for r in rows if r["outcome"] == "perforated"]
        if rebound and perforated:
            low, high = max(rebound), min(perforated)
            if low >= high:
                status, velocities = "ambiguous", []
            else:
                half, step = campaign["threshold"]["refine_half_width_mps"], campaign["threshold"]["refine_step_mps"]
                centre = (low+high)/2
                start, end = int(math.floor((centre-half)/step)*step), int(math.ceil((centre+half)/step)*step)
                velocities = list(range(start, end+1, step))
                status = "ready"
        else:
            status, velocities = "screen_incomplete", []
        cells.append({"cross_rebo": rebo, "cross_lj": lj, "status": status, "screen_rebound_mps": rebound, "screen_perforated_mps": perforated, "velocities_mps": velocities})
    return {"campaign_id": campaign["campaign_id"], "definition": "Threshold repeats are generated only from resolved screen brackets; unresolved cases are not interpolated.", "cells": cells}


def build_submission_gate(results: list[dict[str, Any]], campaign: dict[str, Any], root: Path) -> dict[str, Any]:
    modes = {(r["cross_rebo"], r["cross_lj"]) for r in results if r["stage"] == "mechanism"}
    mechanism = [r for r in results if r["stage"] == "mechanism"]
    structure_ok = all("structure" in r for r in mechanism) and len(mechanism) >= 20
    threshold = [r for r in results if r["stage"] == "threshold"]
    repeats_by_cell_velocity: dict[tuple, int] = defaultdict(int)
    for r in threshold:
        repeats_by_cell_velocity[(r["cross_rebo"], r["cross_lj"], r["velocity_mps"])] += 1
    p2 = bool(repeats_by_cell_velocity) and all(n >= 10 for n in repeats_by_cell_velocity.values())
    convergence = [r for r in results if r["stage"] == "convergence"]
    expected_vars = {"sheet_size_angstrom", "edge_width_angstrom", "timestep_ps", "recovery_ps", "epsilon_eV", "sigma_angstrom", "reactive_model"}
    done_vars = {r["variation"]["name"] for r in convergence}
    doi_file = root / "results" / "publication.json"
    doi_ok = False
    if doi_file.exists():
        doi_ok = bool(read_json(doi_file).get("doi"))
    checks = {
        "P1_strict_2x2_and_shared_initial_state": {"pass": modes == {(0,0),(0,1),(1,0),(1,1)}, "n_mechanism_runs": len(mechanism)},
        "P2_threshold_repeats": {"pass": p2, "replicate_counts": {str(k): v for k, v in repeats_by_cell_velocity.items()}},
        "P3_structure_and_energy_closure": {"pass": structure_ok and all(abs(r["energy_drift_fraction"]) <= 0.01 for r in mechanism), "structure_runs": sum("structure" in r for r in mechanism), "energy_drift_limit": 0.01},
        "P4_convergence_and_provenance": {"pass": expected_vars.issubset(done_vars), "completed_variations": sorted(done_vars), "expected_variations": sorted(expected_vars)},
        "P4_public_dataset_DOI": {"pass": doi_ok, "external_action_required": not doi_ok}
    }
    return {"campaign_id": campaign["campaign_id"], "status": "pass" if all(v["pass"] for v in checks.values()) else "not_ready", "checks": checks, "claim_prohibitions_until_pass": ["Do not call distance contacts bond order, covalent crosslinks, or chemical gating.", "Do not report V50 without resolved repeated binary outcomes and a fitted uncertainty interval.", "Do not call projectile KE loss target absorbed energy."]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", required=True, type=Path)
    parser.add_argument("--stage", default="all", choices=["smoke", "screen", "threshold", "mechanism", "convergence", "all"])
    parser.add_argument("--write-threshold-plan", action="store_true")
    args = parser.parse_args()
    root = args.campaign.resolve().parent.parent
    manifest = read_json(root / "generated" / "run_manifest.json")
    selected_stages = {args.stage} if args.stage != "all" else {"smoke", "screen", "threshold", "mechanism", "convergence"}
    results = []
    missing = []
    for metadata in manifest["runs"]:
        if metadata["stage"] not in selected_stages:
            continue
        result = run_result(root, metadata)
        if result is None:
            missing.append(metadata["run_id"])
        else:
            results.append(result)
    write_json(root / "results" / f"run_results_{args.stage}.json", {"results": results, "missing": missing})
    if args.write_threshold_plan:
        plan = threshold_plan(results, manifest["campaign"])
        write_json(root / "results" / "threshold_plan.json", plan)
    if args.stage == "all":
        all_results = []
        for metadata in manifest["runs"]:
            result = run_result(root, metadata)
            if result is not None:
                all_results.append(result)
        gate = build_submission_gate(all_results, manifest["campaign"], root)
        write_json(root / "results" / "q1_submission_gate.json", gate)
    print(json.dumps({"analysed": len(results), "missing": len(missing), "stage": args.stage}, indent=2))


if __name__ == "__main__":
    main()
