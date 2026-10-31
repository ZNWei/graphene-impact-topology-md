#!/usr/bin/env python3
"""Generate a provenance-locked 2x2 MD campaign for P1-P4.

This script only writes new inputs/manifests in q1_mechanism_pipeline. It never
opens, edits, or executes the legacy Experiment 1-4 scripts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any


FACTOR_LEVELS = (
    (0, 0, "R0_L0"),
    (0, 1, "R0_L1"),
    (1, 0, "R1_L0"),
    (1, 1, "R1_L1"),
)
TRACK_COLUMNS = [
    "step", "time_ps", "stage", "ke_all_eV", "pe_all_eV", "etotal_eV",
    "ke_target_eV", "ke_projectile_eV", "ke_projectile_cm_eV",
    "ke_projectile_nontrans_eV", "pe_target_assigned_eV",
    "pe_projectile_assigned_eV", "e_reactive_eV", "e_zbl_eV", "e_kc_eV",
    "e_lj_tail_eV", "target_zmin_A", "target_zmax_A", "projectile_zmin_A",
    "projectile_zmax_A", "projectile_zcm_A", "projectile_vz_mps"
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def smoothstep5(x: float) -> tuple[float, float]:
    """C2 quintic step and derivative with respect to x, 0 <= x <= 1."""
    x = min(1.0, max(0.0, x))
    return 6*x**5 - 15*x**4 + 10*x**3, 30*x**4 - 60*x**3 + 30*x**2


def generate_lj_table(cfg: dict[str, Any], output: Path) -> dict[str, Any]:
    """Create a short-range-switched LJ tail; force is -dU/dr in eV/A."""
    eps = float(cfg["epsilon_eV"])
    sig = float(cfg["sigma_angstrom"])
    rcut = float(cfg["cutoff_angstrom"])
    r_off = float(cfg["inner_off_angstrom"])
    r_on = float(cfg["inner_on_angstrom"])
    n = int(cfg["grid_points"])
    if not (0.1 < r_off < r_on < rcut):
        raise ValueError("LJ switch requires 0 < inner_off < inner_on < cutoff")
    outer_start = max(r_on + 0.5, rcut - 1.0)
    if outer_start >= rcut:
        raise ValueError("LJ outer switch is degenerate")
    rmin = 0.25
    dr = (rcut - rmin) / (n - 1)

    rows: list[tuple[int, float, float, float]] = []
    max_boundary_force = 0.0
    for i in range(n):
        r = rmin + i * dr
        sr6 = (sig / r) ** 6
        u_lj = 4.0 * eps * (sr6 * sr6 - sr6)
        du_lj = 4.0 * eps * (-12.0 * sr6 * sr6 / r + 6.0 * sr6 / r)
        if r <= r_off:
            s_in, ds_in = 0.0, 0.0
        elif r < r_on:
            s_in, deriv = smoothstep5((r-r_off)/(r_on-r_off))
            ds_in = deriv/(r_on-r_off)
        else:
            s_in, ds_in = 1.0, 0.0
        if r <= outer_start:
            s_out, ds_out = 1.0, 0.0
        else:
            s, deriv = smoothstep5((r-outer_start)/(rcut-outer_start))
            s_out, ds_out = 1.0-s, -deriv/(rcut-outer_start)
        s_total = s_in*s_out
        ds_total = ds_in*s_out + s_in*ds_out
        u = s_total*u_lj
        force = -(ds_total*u_lj + s_total*du_lj)
        if i in (0, n-1):
            max_boundary_force = max(max_boundary_force, abs(force), abs(u))
        rows.append((i+1, r, u, force))

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        fh.write("# Smoothly switched C-C LJ tail for cross-interface 2x2 design\n")
        fh.write(f"# epsilon={eps} eV sigma={sig} A inner=[{r_off},{r_on}] A cutoff={rcut} A\n\n")
        fh.write("CC_LJ_TAIL\n")
        fh.write(f"N {n} R {rmin:.12g} {rcut:.12g}\n\n")
        for idx, r, u, force in rows:
            fh.write(f"{idx} {r:.12g} {u:.16g} {force:.16g}\n")
    audit = {
        "model": "U(r)=S_inner(r)*S_outer(r)*4epsilon[(sigma/r)^12-(sigma/r)^6]",
        "epsilon_eV": eps, "sigma_angstrom": sig, "inner_off_angstrom": r_off,
        "inner_on_angstrom": r_on, "outer_start_angstrom": outer_start,
        "cutoff_angstrom": rcut, "grid_points": n,
        "boundary_energy_or_force_abs_max": max_boundary_force,
        "table_sha256": sha256(output)
    }
    return audit


def pair_block(cross_rebo: int, cross_lj: int, model: str, potentials: dict[str, str]) -> str:
    """Return static LAMMPS pair_style and pair_coeff commands for one factor cell."""
    rebo_file, airebo_file, kc_file = potentials["rebo"], potentials["airebo"], potentials["kc"]
    lj_table = potentials.get("lj_table", "potentials/switched_cc_lj.table")
    if model not in {"rebo", "airebo_rebo_only"}:
        raise ValueError(f"Unknown reactive model: {model}")
    reactive_style = "rebo" if model == "rebo" else "airebo 3.0 0 0"
    reactive_file = rebo_file if model == "rebo" else airebo_file
    if cross_rebo:
        styles = [reactive_style, "zbl 0.95 1.00", "kolmogorov/crespi/z 16.0"]
        if cross_lj:
            styles.append("table linear 20000")
        lines = ["pair_style hybrid/overlay " + " ".join(styles)]
        lines.append(f"pair_coeff * * {reactive_style.split()[0]} {reactive_file} C C C C")
        lines.append("pair_coeff * * zbl 6.0 6.0")
        lines += [
            f"pair_coeff 1 2 kolmogorov/crespi/z {kc_file} C C C NULL",
            f"pair_coeff 2 3 kolmogorov/crespi/z {kc_file} C C C NULL",
            f"pair_coeff 1 3 kolmogorov/crespi/z {kc_file} C C C NULL",
        ]
    else:
        styles = [reactive_style, reactive_style, "zbl 0.95 1.00", "kolmogorov/crespi/z 16.0"]
        if cross_lj:
            styles.append("table linear 20000")
        lines = ["pair_style hybrid/overlay " + " ".join(styles)]
        style_name = reactive_style.split()[0]
        lines += [
            f"pair_coeff * * {style_name} 1 {reactive_file} C C C NULL",
            f"pair_coeff * * {style_name} 2 {reactive_file} NULL NULL NULL C",
            "pair_coeff * * zbl 6.0 6.0",
            f"pair_coeff 1 2 kolmogorov/crespi/z {kc_file} C C C NULL",
            f"pair_coeff 2 3 kolmogorov/crespi/z {kc_file} C C C NULL",
            f"pair_coeff 1 3 kolmogorov/crespi/z {kc_file} C C C NULL",
        ]
    if cross_lj:
        lines += [
            f"pair_coeff 1 4 table {lj_table} CC_LJ_TAIL",
            f"pair_coeff 2 4 table {lj_table} CC_LJ_TAIL",
            f"pair_coeff 3 4 table {lj_table} CC_LJ_TAIL",
        ]
    lines += ["neighbor 2.0 bin", "neigh_modify every 1 delay 0 check yes one 10000"]
    return "\n".join(lines)


def geometry_block(spec: dict[str, Any], *, include_offset: bool = True) -> str:
    b = spec["baseline"]
    sheet, box, d = b["sheet_size_angstrom"], b["box_xy_angstrom"], 3.35
    radius, gap = b["projectile_radius_angstrom"], b["projectile_gap_angstrom"]
    ox, oy = spec["impact_offset_angstrom"] if include_offset else (0.0, 0.0)
    return f"""
variable a equal 2.46
variable b equal v_a*sqrt(3.0)
variable sheet equal {sheet:.12g}
variable boxxy equal {box:.12g}
variable zhalf equal {b['box_zhalf_angstrom']:.12g}
variable d equal {d:.12g}
variable sxlo equal (v_boxxy-v_sheet)/2.0
variable sxhi equal (v_boxxy+v_sheet)/2.0
region sim block 0 ${{boxxy}} 0 ${{boxxy}} -${{zhalf}} ${{zhalf}} units box
create_box 4 sim
lattice custom 1.0 a1 ${{a}} 0 0 a2 0 ${{b}} 0 a3 0 0 100.0 basis 0 0 0 basis 0.5 0.166667 0 basis 0.5 0.5 0 basis 0 0.666667 0
region layer block ${{sxlo}} ${{sxhi}} ${{sxlo}} ${{sxhi}} -0.5 0.5 units box
create_atoms 1 region layer
group l1 type 1
displace_atoms l1 move 0 0 -${{d}} units box
lattice custom 1.0 a1 ${{a}} 0 0 a2 0 ${{b}} 0 a3 0 0 100.0 basis 0 0 0 basis 0.5 0.166667 0 basis 0.5 0.5 0 basis 0 0.666667 0 origin 0.5 0.166667 0
create_atoms 2 region layer
lattice custom 1.0 a1 ${{a}} 0 0 a2 0 ${{b}} 0 a3 0 0 100.0 basis 0 0 0 basis 0.5 0.166667 0 basis 0.5 0.5 0 basis 0 0.666667 0
create_atoms 3 region layer
group l3 type 3
displace_atoms l3 move 0 0 ${{d}} units box
variable cx equal v_boxxy/2.0+{ox:.12g}
variable cy equal v_boxxy/2.0+{oy:.12g}
variable cz equal v_d+{gap + radius:.12g}
region projectile_sphere sphere ${{cx}} ${{cy}} ${{cz}} {radius:.12g} units box
lattice custom 1.0 a1 ${{a}} 0 0 a2 0 ${{b}} 0 a3 0 0 6.70 basis 0 0 0 basis 0.5 0.166667 0 basis 0.5 0.5 0 basis 0 0.666667 0 basis 0.5 0.166667 0.5 basis 0 0.333333 0.5 basis 0 0.666667 0.5 basis 0.5 0.833333 0.5
create_atoms 4 region projectile_sphere
lattice none 1.0
mass * 12.0107
group target type 1 2 3
group projectile type 4
"""


def base_input(spec: dict[str, Any], data_path: str) -> str:
    b = spec["baseline"]
    seed = int(spec["seed"])
    pseed = seed + 100003
    pair = pair_block(1, 0, spec["reactive_model"], spec["potentials"])
    edge = float(b["edge_width_angstrom"])
    return f"""# Generated base state; do not edit by hand. campaign={spec['campaign_id']}\nclear
units metal
dimension 3
boundary p p f
atom_style full
newton on
comm_modify cutoff 20.0
{geometry_block(spec, include_offset=False)}
{pair}
timestep {b['timestep_ps']:.12g}
velocity target create {b['temperature_K']:.12g} {seed} dist gaussian mom yes rot yes
velocity projectile create {b['temperature_K']:.12g} {pseed} dist gaussian mom yes rot yes
velocity projectile zero linear
fix eqt target nvt temp {b['temperature_K']:.12g} {b['temperature_K']:.12g} 0.02
fix eqp projectile nvt temp {b['temperature_K']:.12g} {b['temperature_K']:.12g} 0.02
thermo 100
thermo_style custom step time temp etotal
run {max(1, round(b['relax_ps']/b['timestep_ps']))}
unfix eqt
unfix eqp
variable ew equal {edge:.12g}
region ex1 block ${{sxlo}} $(v_sxlo+v_ew) INF INF INF INF units box
region ex2 block $(v_sxhi-v_ew) ${{sxhi}} INF INF INF INF units box
region ey1 block INF INF ${{sxlo}} $(v_sxlo+v_ew) INF INF units box
region ey2 block INF INF $(v_sxhi-v_ew) ${{sxhi}} INF INF units box
region edges union 4 ex1 ex2 ey1 ey2
group edge_raw region edges
group edge intersect edge_raw target
velocity edge set 0 0 0
write_data {data_path}
print "BASE_COMPLETE {spec['base_id']}" file generated/base/{spec['base_id']}.status
"""


def impact_input(spec: dict[str, Any]) -> str:
    b = spec["baseline"]
    nimpact = max(1, round(((b["projectile_gap_angstrom"] + 2*3.35 + 30.0)/(spec["velocity_mps"]/100.0))/b["timestep_ps"]))
    nrecover = max(1, round(b["recovery_ps"]/b["timestep_ps"]))
    pair = pair_block(spec["cross_rebo"], spec["cross_lj"], spec["reactive_model"], spec["potentials"])
    edge = float(b["edge_width_angstrom"])
    ox, oy = spec["impact_offset_angstrom"]
    trajectory = ""
    if b.get("trajectory_enabled", False):
        trajectory = f"""dump traj all custom {int(b['trajectory_every_steps'])} runs/{spec['run_id']}/trajectory.lammpstrj id type x y z vx vy vz c_pe_atom\ndump_modify traj sort id\nrun 0\n"""
    if spec["cross_rebo"]:
        rname = "rebo" if spec["reactive_model"] == "rebo" else "airebo"
        reactive_compute = f"compute e_reactive all pair {rname}\nvariable e_reactive equal c_e_reactive"
    else:
        rname = "rebo" if spec["reactive_model"] == "rebo" else "airebo"
        reactive_compute = f"compute e_reactive_1 all pair {rname} 1\ncompute e_reactive_2 all pair {rname} 2\nvariable e_reactive equal c_e_reactive_1+c_e_reactive_2"
    lj_compute = "compute e_lj all pair table\nvariable e_lj equal c_e_lj" if spec["cross_lj"] else "variable e_lj equal 0.0"
    return f"""# Generated P1-P4 impact input; run={spec['run_id']}\nclear
units metal
dimension 3
boundary p p f
atom_style full
newton on
comm_modify cutoff 20.0
read_data generated/base/{spec['base_id']}.data
group target type 1 2 3
group projectile type 4
group l1 type 1
group l2 type 2
group l3 type 3
displace_atoms projectile move {ox:.12g} {oy:.12g} 0 units box
variable sheet equal {b['sheet_size_angstrom']:.12g}
variable boxxy equal {b['box_xy_angstrom']:.12g}
variable sxlo equal (v_boxxy-v_sheet)/2.0
variable sxhi equal (v_boxxy+v_sheet)/2.0
variable ew equal {edge:.12g}
region ex1 block ${{sxlo}} $(v_sxlo+v_ew) INF INF INF INF units box
region ex2 block $(v_sxhi-v_ew) ${{sxhi}} INF INF INF INF units box
region ey1 block INF INF ${{sxlo}} $(v_sxlo+v_ew) INF INF units box
region ey2 block INF INF $(v_sxhi-v_ew) ${{sxhi}} INF INF units box
region edges union 4 ex1 ex2 ey1 ey2
group edge_raw region edges
group edge intersect edge_raw target
velocity edge set 0 0 0
fix freeze edge setforce 0 0 0
{pair}
timestep {b['timestep_ps']:.12g}
velocity projectile zero linear
velocity projectile set 0 0 -{spec['velocity_mps']/100.0:.12g} sum yes units box
compute ke_all all ke
compute pe_all all pe
compute ke_target target ke
compute ke_projectile projectile ke
compute pe_atom all pe/atom
compute pe_target target reduce sum c_pe_atom
compute pe_projectile projectile reduce sum c_pe_atom
compute zmin_target target reduce min z
compute zmax_target target reduce max z
compute zmin_projectile projectile reduce min z
compute zmax_projectile projectile reduce max z
{reactive_compute}
compute e_zbl all pair zbl
compute e_kc all pair kolmogorov/crespi/z
{lj_compute}
variable vcmx equal vcm(projectile,x)
variable vcmy equal vcm(projectile,y)
variable vcmz equal vcm(projectile,z)
variable projectile_zcm equal xcm(projectile,z)
variable nprojectile equal count(projectile)
variable ke_projectile_cm equal 0.5*v_nprojectile*12.0107*(v_vcmx*v_vcmx+v_vcmy*v_vcmy+v_vcmz*v_vcmz)*1.0364269e-4
variable ke_projectile_nontrans equal c_ke_projectile-v_ke_projectile_cm
variable etotal equal c_ke_all+c_pe_all
variable stage equal 0
variable time_out equal time
variable vz_mps equal v_vcmz*100.0
thermo 500
thermo_style custom step time c_ke_all c_pe_all v_etotal c_ke_target c_ke_projectile v_ke_projectile_cm v_ke_projectile_nontrans v_vcmz
fix nveall all nve
fix track all ave/time {int(b['track_every_steps'])} 1 {int(b['track_every_steps'])} v_time_out v_stage c_ke_all c_pe_all v_etotal c_ke_target c_ke_projectile v_ke_projectile_cm v_ke_projectile_nontrans c_pe_target c_pe_projectile v_e_reactive c_e_zbl c_e_kc v_e_lj c_zmin_target c_zmax_target c_zmin_projectile c_zmax_projectile v_projectile_zcm v_vcmz file runs/{spec['run_id']}/track.dat title1 "# P1-P4 track" title2 "# units: step ps stage eV eV eV eV eV eV eV eV eV eV eV eV eV A A A A A A/ps"
run 0
print "INITIAL {spec['run_id']} $(v_etotal) $(c_ke_projectile) $(c_ke_target) $(c_pe_all) $(v_ke_projectile_cm) $(v_ke_projectile_nontrans)" file runs/{spec['run_id']}/initial_state.txt
{trajectory}run {nimpact}
run {nrecover}
print "FINAL {spec['run_id']} $(v_etotal) $(v_vcmz) $(c_zmin_target) $(c_zmax_target) $(c_zmin_projectile) $(c_zmax_projectile)" file runs/{spec['run_id']}/final_state.txt
write_data runs/{spec['run_id']}/final.data
print "RUN_COMPLETE {spec['run_id']}" append runs/{spec['run_id']}/final_state.txt
"""


def base_id(spec: dict[str, Any]) -> str:
    b = spec["baseline"]
    raw = {k: b[k] for k in ["sheet_size_angstrom", "box_xy_angstrom", "box_zhalf_angstrom", "projectile_radius_angstrom", "projectile_gap_angstrom", "temperature_K", "timestep_ps", "relax_ps"]}
    raw.update({"seed": spec["seed"], "reactive_model": spec["reactive_model"]})
    digest = hashlib.sha256(json.dumps(raw, sort_keys=True).encode()).hexdigest()[:12]
    return f"base_{digest}"


def make_spec(campaign: dict[str, Any], stage: str, cross_rebo: int, cross_lj: int, velocity: float, seed: int, offset: list[float], overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    baseline = dict(campaign["baseline"])
    lj = dict(campaign["lj_tail"])
    reactive_model = "rebo"
    variation = {"name": "baseline", "value": None}
    if overrides:
        for key, value in overrides.items():
            if key in baseline:
                baseline[key] = value
            elif key in lj:
                lj[key] = value
            elif key == "reactive_model":
                reactive_model = value
            else:
                raise ValueError(f"Unknown override: {key}")
        variation = {"name": next(iter(overrides)), "value": next(iter(overrides.values()))}
    spec: dict[str, Any] = {
        "campaign_id": campaign["campaign_id"], "stage": stage, "cross_rebo": cross_rebo,
        "cross_lj": cross_lj, "velocity_mps": velocity, "seed": seed,
        "impact_offset_angstrom": offset, "baseline": baseline, "lj_tail": lj,
        "reactive_model": reactive_model, "potentials": dict(campaign["potentials"]), "variation": variation,
    }
    table_key = hashlib.sha256(json.dumps(lj, sort_keys=True).encode()).hexdigest()[:12]
    spec["potentials"]["lj_table"] = f"potentials/switched_cc_lj_{table_key}.table"
    spec["base_id"] = base_id(spec)
    label = f"{stage}_{'R'+str(cross_rebo)+'_L'+str(cross_lj)}_v{int(velocity)}_s{seed}_{spec['base_id'][-6:]}"
    if variation["name"] != "baseline":
        label += "_" + str(variation["name"]).replace("_angstrom", "").replace("_", "-") + "-" + str(variation["value"]).replace(".", "p")
    spec["run_id"] = label
    return spec


def stage_specs(campaign: dict[str, Any], stage: str, root: Path, smoke: bool) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    if stage == "screen":
        c = campaign["screen"]
        for rebo, lj, _ in FACTOR_LEVELS:
            for velocity in c["velocities_mps"]:
                for seed, offset in zip(c["seeds"], c["impact_offsets_angstrom"]):
                    specs.append(make_spec(campaign, stage, rebo, lj, velocity, seed, offset))
    elif stage == "threshold":
        plan_path = root / "results" / "threshold_plan.json"
        if not plan_path.exists():
            return []
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        c = campaign["threshold"]
        for cell in plan.get("cells", []):
            if cell.get("status") != "ready":
                continue
            for velocity in cell["velocities_mps"]:
                for i in range(c["replicates"]):
                    specs.append(make_spec(campaign, stage, cell["cross_rebo"], cell["cross_lj"], velocity, c["seed_start"]+i, c["impact_offsets_angstrom"][i]))
    elif stage == "mechanism":
        c = campaign["mechanism"]
        for rebo, lj, _ in FACTOR_LEVELS:
            for seed, offset in zip(c["seeds"], c["impact_offsets_angstrom"]):
                specs.append(make_spec(campaign, stage, rebo, lj, c["velocity_mps"], seed, offset))
    elif stage == "convergence":
        c = campaign["convergence"]
        dimensions = [
            ("sheet_size_angstrom", c["sheet_sizes_angstrom"]),
            ("edge_width_angstrom", c["edge_widths_angstrom"]),
            ("timestep_ps", c["timesteps_ps"]),
            ("recovery_ps", c["recovery_ps"]),
            ("epsilon_eV", c["lj_epsilons_eV"]),
            ("sigma_angstrom", c["lj_sigmas_angstrom"]),
            ("reactive_model", c["reactive_models"]),
        ]
        for dim, values in dimensions:
            for value in values:
                for seed, offset in zip(c["seeds"], c["impact_offsets_angstrom"]):
                    specs.append(make_spec(campaign, stage, 1, 1, c["velocity_mps"], seed, offset, {dim: value}))
    elif stage == "smoke":
        smoke_overrides = {"sheet_size_angstrom": 20.0, "box_xy_angstrom": 50.0, "box_zhalf_angstrom": 100.0, "projectile_radius_angstrom": 4.0, "projectile_gap_angstrom": 10.0, "timestep_ps": 0.0002, "relax_ps": 0.002, "recovery_ps": 0.003, "track_every_steps": 5, "trajectory_every_steps": 5, "trajectory_enabled": True}
        for rebo, lj, _ in FACTOR_LEVELS:
            specs.append(make_spec(campaign, stage, rebo, lj, 10000, 71001, [0.0, 0.0], smoke_overrides))
    return specs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", required=True, type=Path)
    parser.add_argument("--smoke", action="store_true", help="generate only the four short syntax tests")
    args = parser.parse_args()
    campaign_path = args.campaign.resolve()
    root = campaign_path.parent.parent
    campaign = json.loads(campaign_path.read_text(encoding="utf-8"))
    required = [root / value for value in campaign["potentials"].values()]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing versioned potential files: " + "; ".join(missing))
    generated, runs, potentials = root / "generated", root / "runs", root / "potentials"
    for folder in [generated / "base", generated / "base_inputs", generated / "inputs", generated / "metadata", runs, root / "results", potentials]:
        folder.mkdir(parents=True, exist_ok=True)
    stages = ["smoke"] if args.smoke else ["screen", "threshold", "mechanism", "convergence"]
    specs: list[dict[str, Any]] = []
    for stage in stages:
        specs.extend(stage_specs(campaign, stage, root, args.smoke))
    lj_audits: dict[str, dict[str, Any]] = {}
    for spec in specs:
        table_relative = spec["potentials"]["lj_table"]
        table_path = root / table_relative
        if table_relative not in lj_audits:
            lj_audits[table_relative] = generate_lj_table(spec["lj_tail"], table_path)
    write_json(root / "results" / "lj_tail_audit.json", lj_audits)
    bases: dict[str, dict[str, Any]] = {}
    for spec in specs:
        bases.setdefault(spec["base_id"], spec)
    for bid, spec in bases.items():
        base_path = generated / "base_inputs" / f"{bid}.in"
        base_path.write_text(base_input(spec, f"generated/base/{bid}.data"), encoding="utf-8")
    manifest_runs = []
    for spec in specs:
        input_path = generated / "inputs" / f"{spec['run_id']}.in"
        input_path.write_text(impact_input(spec), encoding="utf-8")
        run_dir = runs / spec["run_id"]
        run_dir.mkdir(parents=True, exist_ok=True)
        metadata = dict(spec)
        base_input_path = generated / "base_inputs" / (spec["base_id"] + ".in")
        metadata.update({"input_file": str(input_path.relative_to(root)), "input_sha256": sha256(input_path), "base_input_file": str(base_input_path.relative_to(root)), "track_columns": TRACK_COLUMNS})
        write_json(generated / "metadata" / f"{spec['run_id']}.json", metadata)
        manifest_runs.append(metadata)
    potential_hashes = {name: sha256(root / rel) for name, rel in campaign["potentials"].items()}
    potential_hashes.update({relative: audit["table_sha256"] for relative, audit in lj_audits.items()})
    manifest = {"campaign": campaign, "campaign_sha256": sha256(campaign_path), "potential_sha256": potential_hashes, "lj_tail_audit": lj_audits, "runs": manifest_runs, "base_inputs": sorted(bases)}
    write_json(generated / "run_manifest.json", manifest)
    print(f"Generated {len(bases)} base inputs and {len(specs)} impact inputs under {root}")
    if not specs and not args.smoke:
        print("No threshold inputs yet: run and analyse the screen stage first.")


if __name__ == "__main__":
    main()
