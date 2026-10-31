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

if __name__=="__main__":
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"Arial","font.size":7.5,"axes.labelsize":7.5,"axes.titlesize":8.3,"xtick.labelsize":6.8,"ytick.labelsize":6.8,"pdf.fonttype":42,"ps.fonttype":42})
    data=json.loads((ROOT/"results/source/p1_p4_results_all.json").read_text(encoding="utf-8"))
    print(figure_sensitivity(data))
