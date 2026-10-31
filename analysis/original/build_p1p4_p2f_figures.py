#!/usr/bin/env python3
"""Create provenance-linked, claim-safe figures from completed P1-P4 and P2-F data."""
from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
P1 = ROOT / "P1_P4_UNFINISHED_Workstation_Package_v2/q1_mechanism_pipeline/results/run_results_all.json"
P2F = ROOT / "P2_F_FRAGMENT_TRANSMISSION_WORKSTATION_PACKAGE_2026-07-29/q1_mechanism_pipeline/p2_fragment_transmission_pipeline/results/p2f_transmission_results.json"
OUT = ROOT / "q1_revision_assets/p1p4_p2f"

COLORS = {(0, 0): "#4C78A8", (0, 1): "#72B7B2", (1, 0): "#F58518", (1, 1): "#54A24B"}
LABELS = {(0, 0): "R0 L0", (0, 1): "R0 L1", (1, 0): "R1 L0", (1, 1): "R1 L1"}
ORDER = [(0, 0), (0, 1), (1, 0), (1, 1)]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tidy(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)


def mean_sd(values: list[float]) -> tuple[float, float]:
    mean = sum(values) / len(values)
    sd = math.sqrt(sum((x - mean) ** 2 for x in values) / (len(values) - 1)) if len(values) > 1 else 0.0
    return mean, sd


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.16, 1.08, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


def figure_factorial(p1: dict) -> Path:
    groups: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for row in p1["results"]:
        if row["stage"] == "mechanism":
            groups[(row["cross_rebo"], row["cross_lj"])].append(row)
    metrics = [
        ("Projectile COM kinetic-energy loss (eV)", lambda r: r["projectile_com_ke_loss_eV"]),
        ("Target kinetic energy at terminal time (eV)", lambda r: r["final"]["ke_target_eV"]),
        ("Net target undercoordination surrogate (atoms)", lambda r: r["structure"]["target_damage_edge_excluded_net"]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(7.25, 2.55), sharex=True)
    jitters = [-0.16, -0.08, 0.0, 0.08, 0.16]
    for ax, (ylabel, getter), letter in zip(axes, metrics, "ABC"):
        for x, key in enumerate(ORDER):
            values = [getter(row) for row in groups[key]]
            for j, value in zip(jitters, values):
                ax.scatter(x+j, value, s=22, marker="o", facecolor="white", edgecolor=COLORS[key], linewidth=1.1, zorder=3)
            mean, sd = mean_sd(values)
            ax.errorbar(x, mean, yerr=sd, fmt="D", color=COLORS[key], markeredgecolor="black", markersize=4.2, capsize=2.4, linewidth=1.1, zorder=4)
        ax.set_ylabel(ylabel)
        ax.set_xticks(range(4), ["R0\nL0", "R0\nL1", "R1\nL0", "R1\nL1"])
        ax.set_xlim(-0.55, 3.55)
        tidy(ax); panel_label(ax, letter)
    axes[0].set_title("Energy redirection")
    axes[1].set_title("Target kinetic response")
    axes[2].set_title("Structural-response surrogate")
    fig.text(0.5, -0.02, "Cross-interface reactive topology (R) and switched LJ-tail state (L)", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    png = OUT / "Fig1_factorial_energy_and_structure.png"
    pdf = OUT / "Fig1_factorial_energy_and_structure.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png


def wilson(successes: int, trials: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = successes / trials
    denominator = 1 + z*z/trials
    center = (p + z*z/(2*trials)) / denominator
    half = z * math.sqrt(p*(1-p)/trials + z*z/(4*trials*trials)) / denominator
    return max(0.0, center-half), min(1.0, center+half)


def figure_p2f(p2f: dict) -> Path:
    groups: dict[tuple[int, int, int], list[dict]] = defaultdict(list)
    for row in p2f["valid_results"]:
        groups[(row["cross_rebo"], row["cross_lj"], int(row["velocity_mps"]))].append(row)
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 2.8), gridspec_kw={"width_ratios": [1.18, 1.0]})
    ax = axes[0]
    jitters = [-120, -84, -48, -12, 12, 48, 84, 120, -30, 30]
    for key2 in ORDER:
        points = sorted((velocity, groups[(key2[0], key2[1], velocity)]) for r, l, velocity in groups if (r, l) == key2)
        xs, means = [], []
        for velocity, rows in points:
            values = [row["primary_fraction"] for row in rows]
            for jitter, value in zip(jitters, values):
                ax.scatter(velocity+jitter, value, s=16, facecolor="white", edgecolor=COLORS[key2], linewidth=0.85, zorder=3)
            mean, sd = mean_sd(values)
            ax.errorbar(velocity, mean, yerr=sd, fmt="o", color=COLORS[key2], markeredgecolor="black", markersize=4, capsize=2, linewidth=1.05, zorder=4)
            xs.append(velocity); means.append(mean)
        ax.plot(xs, means, color=COLORS[key2], linewidth=1.0, label=LABELS[key2])
    ax.axhline(0.50, color="#555555", linewidth=0.9, linestyle="--")
    ax.text(19000, 0.515, "pre-registered event boundary", fontsize=6.6, color="#444444")
    ax.set(xlim=(19500, 24500), ylim=(0, 0.75), xlabel="Impact velocity (km s$^{-1}$)", ylabel="Terminal projectile-origin downstream fraction")
    ax.set_xticks([20000, 22000, 24000], ["20", "22", "24"])
    ax.legend(frameon=False, fontsize=6.5, title="Factor cell", title_fontsize=6.5, ncol=2, loc="upper left")
    tidy(ax); panel_label(ax, "A")
    ax = axes[1]
    for key2 in ORDER:
        points = sorted((velocity, groups[(key2[0], key2[1], velocity)]) for r, l, velocity in groups if (r, l) == key2)
        xs, probs, lows, highs = [], [], [], []
        for velocity, rows in points:
            events = sum(bool(row["primary_event"]) for row in rows)
            lo, hi = wilson(events, len(rows))
            probability = events/len(rows)
            xs.append(velocity); probs.append(probability); lows.append(max(0.0, probability-lo)); highs.append(max(0.0, hi-probability))
        ax.errorbar(xs, probs, yerr=[lows, highs], marker="o", markersize=4, capsize=2.2, color=COLORS[key2], linewidth=1.05, label=LABELS[key2])
    ax.axhline(0.50, color="#555555", linewidth=0.9, linestyle="--")
    ax.set(xlim=(19500, 24500), ylim=(-0.04, 1.04), xlabel="Impact velocity (km s$^{-1}$)", ylabel="Observed event fraction (n = 10)")
    ax.set_xticks([20000, 22000, 24000], ["20", "22", "24"])
    tidy(ax); panel_label(ax, "B")
    fig.tight_layout()
    png = OUT / "Fig2_p2f_projectile_origin_mass_transmission.png"
    pdf = OUT / "Fig2_p2f_projectile_origin_mass_transmission.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png


def figure_sensitivity(p1: dict) -> Path:
    reference = []
    variations: dict[str, list[tuple[str, list[float]]]] = defaultdict(list)
    for row in p1["results"]:
        if row["stage"] != "convergence":
            continue
        name, value = row["variation"]["name"], str(row["variation"]["value"])
        variations[(name, value)].append(row["projectile_com_ke_loss_eV"])
        if name == "sheet_size_angstrom" and value == "60.0":
            reference.append(row["projectile_com_ke_loss_eV"])
    reference_mean, _ = mean_sd(reference)
    labels, values, errors = [], [], []
    for (name, value), rows in sorted(variations.items()):
        mean, sd = mean_sd(rows)
        labels.append(f"{name.replace('_angstrom', '').replace('_eV', '').replace('_ps', '')}: {value}")
        values.append(100*(mean-reference_mean)/reference_mean)
        errors.append(100*sd/reference_mean)
    fig, ax = plt.subplots(figsize=(6.75, 4.45))
    ypos = list(range(len(labels)))
    ax.errorbar(values, ypos, xerr=errors, fmt="o", color="#4C78A8", markeredgecolor="black", markersize=4, capsize=2.3, linewidth=1.0)
    ax.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    ax.set_yticks(ypos, labels, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("Change in mean projectile COM kinetic-energy loss relative to 60 Å reference (%)")
    ax.set_title("P4 numerical and model-parameter sensitivity (n = 3 per setting)")
    tidy(ax)
    fig.tight_layout()
    png = OUT / "Fig3_p4_sensitivity.png"
    pdf = OUT / "Fig3_p4_sensitivity.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png


def graphical_abstract() -> Path:
    fig, ax = plt.subplots(figsize=(7.25, 2.2))
    ax.set_axis_off(); ax.set_xlim(0, 10); ax.set_ylim(0, 3)
    for y in (0.7, 0.98, 1.26):
        ax.plot([0.6, 3.55], [y, y], color="#3B7EA1", linewidth=3.4, solid_capstyle="round")
    circle = plt.Circle((2.1, 2.2), 0.30, facecolor="#5F6B7A", edgecolor="black", linewidth=0.8)
    ax.add_patch(circle); ax.annotate("", xy=(2.1, 1.48), xytext=(2.1, 1.93), arrowprops=dict(arrowstyle="-|>", lw=1.4, color="#333333"))
    ax.text(0.55, 2.72, "Multilayer graphene + carbon projectile", fontsize=9, weight="bold")
    ax.text(0.55, 0.22, "2 × 2 cross-interface topology", fontsize=8)
    ax.annotate("", xy=(4.65, 1.5), xytext=(3.75, 1.5), arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#333333"))
    ax.add_patch(plt.Rectangle((4.75, 0.55), 2.1, 1.9, fill=False, edgecolor="#555555", linewidth=1.1))
    ax.text(5.8, 2.15, "Matched MD", ha="center", fontsize=9, weight="bold")
    ax.text(5.8, 1.68, "R: off/on", ha="center", fontsize=7.8)
    ax.text(5.8, 1.38, "L: off/on", ha="center", fontsize=7.8)
    ax.text(5.8, 0.88, "5 matched replicates", ha="center", fontsize=7.3)
    ax.annotate("", xy=(7.65, 1.5), xytext=(6.95, 1.5), arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#333333"))
    ax.text(8.7, 2.28, "Energy redirection", ha="center", fontsize=8.5, weight="bold")
    ax.text(8.7, 1.78, "structural-response\nsurrogates", ha="center", fontsize=7.3)
    ax.text(8.7, 0.96, "P2-F: projectile-origin\nmass transmission", ha="center", fontsize=7.3)
    ax.text(8.7, 0.22, "No chemical-bond or V50 claim", ha="center", fontsize=6.7, color="#7A3E00")
    fig.tight_layout(pad=0.2)
    png = OUT / "Graphical_Abstract_P1P4_P2F.png"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return png


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 8.3, "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "pdf.fonttype": 42, "ps.fonttype": 42})
    p1, p2f = json.loads(P1.read_text(encoding="utf-8")), json.loads(P2F.read_text(encoding="utf-8"))
    figures = [graphical_abstract(), figure_factorial(p1), figure_p2f(p2f), figure_sensitivity(p1)]
    manifest = {"p1_p4_source": str(P1.relative_to(ROOT)), "p1_p4_sha256": sha256(P1), "p2f_source": str(P2F.relative_to(ROOT)), "p2f_sha256": sha256(P2F), "figures": [str(path.relative_to(ROOT)) for path in figures], "claims": ["P1 is a 2x2 factorial comparison with matched seed/offset base states.", "P2-F reports terminal projectile-origin mass fractions and observed event proportions, not a conventional V50.", "Distance-based structural surrogates are not bond-order evidence."]}
    (OUT / "figure_data_script_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"created": [path.name for path in figures], "directory": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
