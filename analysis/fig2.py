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
ORDER = ["R0_L0", "R0_L1", "R1_L0", "R1_L1"]
LABELS = {key: key.replace("_", " ") for key in ORDER}
COLORS = {
    "R0_L0": "#0072B2",
    "R0_L1": "#56B4E9",
    "R1_L0": "#E69F00",
    "R1_L1": "#D55E00",
}
MARKERS = {"R0_L0": "o", "R0_L1": "s", "R1_L0": "^", "R1_L1": "D"}

def tidy(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)

def logistic(x_kms: float, midpoint_kms: float, slope_per_1000_mps: float) -> float:
    eta = slope_per_1000_mps * (x_kms - midpoint_kms)
    if eta >= 0:
        return 1.0 / (1.0 + math.exp(-eta))
    exp_eta = math.exp(eta)
    return exp_eta / (1.0 + exp_eta)

def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.14, 1.06, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")

def threshold_figure(results: dict) -> tuple[Path, Path]:
    estimates = results["threshold_estimates"]
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.05), gridspec_kw={"width_ratios": [1.55, 1.0]})

    ax = axes[0]
    for key in ORDER:
        estimate = estimates[key]
        groups = estimate["event_proportions"]
        xs = [g["velocity_mps"] / 1000.0 for g in groups]
        ys = [g["event_probability"] for g in groups]
        lows = [max(0.0, y - g["wilson_95_ci"][0]) for y, g in zip(ys, groups)]
        highs = [max(0.0, g["wilson_95_ci"][1] - y) for y, g in zip(ys, groups)]
        ax.errorbar(
            xs,
            ys,
            yerr=[lows, highs],
            fmt=MARKERS[key],
            markersize=4.8,
            markerfacecolor="white",
            markeredgewidth=1.1,
            capsize=2.4,
            color=COLORS[key],
            linewidth=1.0,
            label=LABELS[key],
            zorder=3,
        )
        grid = [min(xs) + (max(xs) - min(xs)) * i / 160.0 for i in range(161)]
        midpoint = estimate["v_t_50_mps"] / 1000.0
        slope = estimate["slope_per_1000_mps"]
        ax.plot(grid, [logistic(x, midpoint, slope) for x in grid], color=COLORS[key], linewidth=1.35)

    ax.axhline(0.5, color="#555555", linestyle="--", linewidth=0.9)
    ax.text(23.93, 0.525, "event probability = 0.5", ha="right", va="bottom", fontsize=6.5, color="#444444")
    ax.set_xlim(19.75, 24.25)
    ax.set_ylim(-0.04, 1.04)
    ax.set_xticks([20, 20.5, 21, 21.5, 22, 24])
    ax.set_xticklabels(["20", "20.5", "21", "21.5", "22", "24"])
    ax.set_xlabel("Impact velocity (km s$^{-1}$)")
    ax.set_ylabel("Projectile-origin mass-transmission event fraction")
    ax.legend(frameon=False, ncol=2, fontsize=6.7, loc="lower right", columnspacing=1.0, handletextpad=0.4)
    ax.set_title("Observed proportions and pooled grouped-logistic fits")
    tidy(ax)
    panel_label(ax, "A")

    ax = axes[1]
    y_positions = list(reversed(range(len(ORDER))))
    for y, key in zip(y_positions, ORDER):
        estimate = estimates[key]
        value = estimate["v_t_50_mps"] / 1000.0
        low, high = [v / 1000.0 for v in estimate["wald_95_ci_mps"]]
        ax.errorbar(
            value,
            y,
            xerr=[[value - low], [high - value]],
            fmt=MARKERS[key],
            markersize=5.3,
            markerfacecolor=COLORS[key],
            markeredgecolor="black",
            markeredgewidth=0.6,
            color=COLORS[key],
            capsize=3.0,
            linewidth=1.4,
            zorder=3,
        )
        ax.text(high + 0.035, y, f"{value:.3f}", va="center", fontsize=6.8, color="#222222")
    ax.set_yticks(y_positions, [LABELS[key] for key in ORDER])
    ax.set_xlim(20.0, 22.45)
    ax.set_ylim(-0.65, 3.65)
    ax.set_xlabel("$V_{T,50}$ (km s$^{-1}$)")
    ax.set_title("Point estimates and Wald 95% CIs")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.7)
    tidy(ax)
    panel_label(ax, "B")

    fig.tight_layout(w_pad=2.2)
    png = OUT / "Fig2_p2f_vt50_threshold_complete.png"
    pdf = OUT / "Fig2_p2f_vt50_threshold_complete.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf

if __name__=="__main__":
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"Arial","font.size":7.5,"axes.labelsize":7.5,"axes.titlesize":8.3,"xtick.labelsize":6.8,"ytick.labelsize":6.8,"pdf.fonttype":42,"ps.fonttype":42})
    data=json.loads((ROOT/"results/source/p2f_round2_results.json").read_text(encoding="utf-8"))
    print(threshold_figure(data))
