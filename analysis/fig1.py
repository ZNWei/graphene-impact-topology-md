"""Recreate manuscript figure using the original plotting functions.
Only portable paths, imports and the single-figure entry point were adapted.
Requires matplotlib; numeric re-analysis uses the Python standard library.
"""
from __future__ import annotations
import json,math
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/recreated_figures"
COLORS = {(0, 0): "#4C78A8", (0, 1): "#72B7B2", (1, 0): "#F58518", (1, 1): "#54A24B"}
LABELS = {(0, 0): "R0 L0", (0, 1): "R0 L1", (1, 0): "R1 L0", (1, 1): "R1 L1"}
ORDER = [(0, 0), (0, 1), (1, 0), (1, 1)]

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

if __name__=="__main__":
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"Arial","font.size":7.5,"axes.labelsize":7.5,"axes.titlesize":8.3,"xtick.labelsize":6.8,"ytick.labelsize":6.8,"pdf.fonttype":42,"ps.fonttype":42})
    data=json.loads((ROOT/"results/source/p1_p4_results_all.json").read_text(encoding="utf-8"))
    print(figure_factorial(data))
