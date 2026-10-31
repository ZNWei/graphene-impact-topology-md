#!/usr/bin/env python3
"""Build final, provenance-linked P2-F threshold figures without rerunning MD."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
RESULTS = (
    ROOT
    / "P2F_VT50_ADAPTIVE_RESOLUTION_WORKSTATION_PACKAGE_2026-07-30"
    / "q1_mechanism_pipeline"
    / "P2F_VT50_ANALYSIS_CORRECTION_2026-07-31"
    / "q1_mechanism_pipeline"
    / "p2f_vt50_round2_pipeline"
    / "results"
)
ROUND2 = RESULTS / "p2f_vt50_round2_results.json"
GATE = RESULTS / "p2f_vt50_round2_gate.json"
SUMMARY = RESULTS / "p2f_vt50_round2_summary.csv"
OUT = ROOT / "q1_revision_assets" / "p1p4_p2f"

ORDER = ["R0_L0", "R0_L1", "R1_L0", "R1_L1"]
LABELS = {key: key.replace("_", " ") for key in ORDER}
COLORS = {
    "R0_L0": "#0072B2",
    "R0_L1": "#56B4E9",
    "R1_L0": "#E69F00",
    "R1_L1": "#D55E00",
}
MARKERS = {"R0_L0": "o", "R0_L1": "s", "R1_L0": "^", "R1_L1": "D"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def graphical_abstract(results: dict) -> Path:
    estimates = results["threshold_estimates"]
    fig, ax = plt.subplots(figsize=(7.25, 2.35))
    ax.set_axis_off()
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 3.4)

    box_kw = dict(facecolor="#FAFAFA", edgecolor="#6B6B6B", linewidth=1.0)
    ax.add_patch(plt.Rectangle((0.25, 0.62), 2.35, 2.15, **box_kw))
    ax.add_patch(plt.Rectangle((3.15, 0.62), 2.35, 2.15, **box_kw))
    ax.add_patch(plt.Rectangle((6.05, 0.62), 2.55, 2.15, **box_kw))
    ax.add_patch(plt.Rectangle((9.15, 0.62), 2.6, 2.15, **box_kw))
    for x0, x1 in [(2.6, 3.15), (5.5, 6.05), (8.6, 9.15)]:
        ax.annotate("", xy=(x1 - 0.05, 1.70), xytext=(x0 + 0.05, 1.70), arrowprops=dict(arrowstyle="-|>", lw=1.25, color="#333333"))

    ax.text(1.425, 2.48, "Hypervelocity MD", ha="center", fontsize=9.2, fontweight="bold")
    for y in (1.25, 1.47, 1.69):
        ax.plot([0.62, 2.22], [y, y], color="#3B7EA1", linewidth=2.7, solid_capstyle="round")
    ax.add_patch(plt.Circle((1.42, 2.05), 0.18, facecolor="#5F6B7A", edgecolor="black", linewidth=0.6))
    ax.annotate("", xy=(1.42, 1.80), xytext=(1.42, 1.94), arrowprops=dict(arrowstyle="-|>", lw=1.0, color="#333333"))
    ax.text(1.425, 0.88, "3-layer graphene +\ncarbon projectile", ha="center", fontsize=6.8)

    ax.text(4.325, 2.48, "Matched 2 × 2 topology", ha="center", fontsize=9.0, fontweight="bold")
    ax.text(4.325, 2.04, "R: cross-interface reactive\nstate off/on", ha="center", fontsize=6.8)
    ax.text(4.325, 1.38, "L: switched LJ-tail\nstate off/on", ha="center", fontsize=6.8)
    ax.text(4.325, 0.84, "identical seed/offset pairs", ha="center", fontsize=6.5, color="#555555")

    ax.text(7.325, 2.48, "Energy redirection", ha="center", fontsize=9.0, fontweight="bold")
    ax.text(7.325, 2.02, "L effect at R0: −407.2 eV", ha="center", fontsize=7.0, color=COLORS["R0_L0"])
    ax.text(7.325, 1.60, "L effect at R1: −0.6 eV", ha="center", fontsize=7.0, color=COLORS["R1_L1"])
    ax.text(7.325, 1.12, "non-additive response", ha="center", fontsize=7.1, fontweight="bold")
    ax.text(7.325, 0.82, "interaction contrast: +406.6 eV", ha="center", fontsize=6.4, color="#555555")

    ax.text(10.45, 2.48, "Origin-resolved $V_{T,50}$", ha="center", fontsize=9.0, fontweight="bold")
    y0 = 2.04
    for i, key in enumerate(ORDER):
        value = estimates[key]["v_t_50_mps"] / 1000.0
        ax.text(9.42, y0 - 0.34 * i, LABELS[key], ha="left", va="center", fontsize=6.8, color=COLORS[key], fontweight="bold")
        ax.text(11.48, y0 - 0.34 * i, f"{value:.3f} km s$^{{-1}}$", ha="right", va="center", fontsize=6.8)

    ax.text(6.0, 0.18, "Fragmentation-mediated projectile-origin mass transmission; not chemical-bond evidence and not conventional ballistic V50", ha="center", fontsize=6.8, color="#7A3E00", fontweight="bold")
    fig.tight_layout(pad=0.15)
    png = OUT / "Graphical_Abstract_P1P4_P2F_threshold_complete.png"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return png


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.4,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.2,
            "xtick.labelsize": 6.7,
            "ytick.labelsize": 6.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    results = json.loads(ROUND2.read_text(encoding="utf-8"))
    gate = json.loads(GATE.read_text(encoding="utf-8"))
    if gate.get("status") != "pass":
        raise RuntimeError("P2-F pooled V_T,50 gate is not pass")
    if results.get("pooled_run_count") != 150:
        raise RuntimeError("Expected 150 pooled P2-F trajectories")
    if len(results.get("new_valid_results", [])) != 30 or results.get("new_invalid") or results.get("new_incomplete"):
        raise RuntimeError("Round-2 completion checks failed")
    if any(results["threshold_estimates"][key].get("status") != "estimated" for key in ORDER):
        raise RuntimeError("At least one factor cell lacks an estimated V_T,50")

    fig2_png, fig2_pdf = threshold_figure(results)
    graphical_png = graphical_abstract(results)
    manifest = {
        "round2_results_source": str(ROUND2.relative_to(ROOT)),
        "round2_results_sha256": sha256(ROUND2),
        "round2_gate_source": str(GATE.relative_to(ROOT)),
        "round2_gate_sha256": sha256(GATE),
        "round2_summary_source": str(SUMMARY.relative_to(ROOT)),
        "round2_summary_sha256": sha256(SUMMARY),
        "script": Path(__file__).name,
        "script_sha256": sha256(Path(__file__)),
        "outputs": [str(path.relative_to(ROOT)) for path in (fig2_png, fig2_pdf, graphical_png)],
        "evidence_boundary": [
            "The passed gate applies only to pooled P2-F V_T,50 estimability.",
            "V_T,50 is a fragmentation-mediated projectile-origin downstream mass-transmission threshold.",
            "It is not conventional ballistic V50, intact-projectile perforation, or chemical-bond evidence.",
        ],
    }
    manifest_path = OUT / "figure_data_script_manifest_v1_3.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"created": manifest["outputs"], "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
