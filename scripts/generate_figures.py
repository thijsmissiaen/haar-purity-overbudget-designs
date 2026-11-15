"""Generate publication-quality figures strictly from scripts/verified_results.json."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA_PATH = Path(__file__).resolve().parent / "verified_results.json"
FIG_DIR = ROOT / "paper" / "figures"


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
        "mathtext.fontset": "cm",
        "font.size": 9.5,
        "axes.titlesize": 10.0,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 300,
        "savefig.dpi": 300,
    })


def make_figure1(data: dict) -> None:
    sep12 = data["lp_and_witt_certificates"]["separation_certificates"]["12"]
    opt12 = data["lp_and_witt_certificates"]["optimal_n12_families"]

    ms_sep = [row["m"] for row in sep12]
    r2 = [row["ratio_d2"] for row in sep12]
    r3 = [row["ratio_d3"] for row in sep12]
    r10 = [row["ratio_d10"] for row in sep12]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)

    # Panel (a): Separation ratio S_{12,d}(F_Q) / S_{12,d}(F_H)
    x = np.arange(len(ms_sep))
    w = 0.25
    ax1.bar(x - w, r2, width=w, color="#2b5c8f", edgecolor="#1a3857", linewidth=0.5, label=r"$d=2$")
    ax1.bar(x, r3, width=w, color="#4f8fc0", edgecolor="#2a5577", linewidth=0.5, label=r"$d=3$")
    ax1.bar(x + w, r10, width=w, color="#c85a32", edgecolor="#7d3319", linewidth=0.5, label=r"$d=10$")
    ax1.axhline(1.0, color="#333333", linestyle=":", linewidth=0.9)
    # Mark boundary between m=60 (index 3) and m=140 (index 4) -> m=66 Witt threshold
    ax1.axvline(3.5, color="#555555", linestyle="--", linewidth=1.0)
    ax1.text(3.58, 9.5, r"$m=66$" "\n" r"barrier", ha="left", va="top", fontsize=8.0, color="#444444")

    ax1.set_yscale("log")
    ax1.set_ylim(0.92, 13.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels([str(m) for m in ms_sep])
    ax1.set_xlabel(r"Cut budget $m$ ($n=12$)")
    ax1.set_ylabel(r"Ratio $S_{12,d}(F_Q)\,/\,S_{12,d}(F_H)$ (log scale)")
    ax1.set_title(r"(a) $E(s^2)$ vs. Haar separation across local dimensions $d$")
    ax1.legend(frameon=False, loc="upper left", ncol=3)
    ax1.grid(axis="y", linestyle=":", alpha=0.35)

    # Annotate peak at m=60, d=10
    ax1.annotate(
        f"{r10[3]:.2f}$\\times$",
        xy=(3 + w, r10[3]),
        xytext=(3 + w - 0.05, r10[3] * 1.12),
        ha="center",
        va="bottom",
        fontsize=7.8,
        color="#7d3319",
    )

    # Panel (b): Adjacent pair fraction c_8 / [m(m-1)/2] vs Delsarte LP floor
    m_grid = np.arange(11, 463)
    lp_frac = np.where(m_grid <= 66, 0.0, (m_grid * (m_grid - 66)) / (22.0 * m_grid * (m_grid - 1) / 2.0))
    ax2.plot(m_grid, lp_frac, color="#1f4e79", linewidth=1.8, label=r"Delsarte LP floor $\max(0,\,(m-66)/(11(m-1)))$")

    # W_1 union S linear wall on (66, 132)
    m_wall = np.arange(66, 133)
    wall_frac = 6.0 * (m_wall - 66) / (m_wall * (m_wall - 1) / 2.0)
    ax2.plot(m_wall, wall_frac, color="#c85a32", linestyle="--", linewidth=1.4, label=r"$W_1 \cup S$ extension floor $12(m-66)/(m(m-1))$")

    # Exact optimal families
    opt_ms = sorted(int(k) for k in opt12.keys())
    opt_fracs = [
        opt12[str(m)]["counts"].get("8", 0) / (m * (m - 1) / 2.0)
        for m in opt_ms
    ]
    ax2.scatter(opt_ms, opt_fracs, color="#1a7f37", edgecolor="#0d441d", zorder=5, s=46, marker="o", label="Certified universal optima")
    for m_pt, f_pt in zip(opt_ms, opt_fracs):
        if m_pt in (132, 330, 396, 462):
            ax2.text(m_pt, f_pt + 0.0045, str(m_pt), ha="center", va="bottom", fontsize=7.5, color="#0d441d")
        elif m_pt == 66:
            ax2.text(m_pt + 10, f_pt + 0.0035, "66", ha="left", va="bottom", fontsize=7.5, color="#0d441d")

    ax2.axvline(66, color="#555555", linestyle=":", linewidth=0.9)
    ax2.set_xlim(0, 480)
    ax2.set_ylim(-0.006, 0.095)
    ax2.set_xlabel(r"Cut budget $m$ ($n=12$, $M_{12}=462$)")
    ax2.set_ylabel(r"Adjacent-pair fraction $2c_8(F)\,/\,(m(m-1))$")
    ax2.set_title(r"(b) $S(5,6,12)/\mathbb{Z}_2$ phase transition and certified optima")
    ax2.legend(frameon=False, loc="upper left", fontsize=7.8)
    ax2.grid(axis="both", linestyle=":", alpha=0.35)

    fig.savefig(FIG_DIR / "fig1_es2_vs_haar_n12.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig1_es2_vs_haar_n12.png", bbox_inches="tight")
    plt.close(fig)


def make_figure2(data: dict) -> None:
    spec = data["folded_johnson_spectral"]["spectral_summary"]
    asymp = data["overbudget_asymptotics"]["asymptotic_prefactor_convergence"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)

    # Panel (a): Normalized eigenvalues Lambda_i(n,2) / Lambda_0(n,2) for n=12, 14, 16, 18, 20
    styles = [
        ("12", "#1f4e79", "o", "-"),
        ("16", "#2c8259", "s", "--"),
        ("20", "#c85a32", "^", "-."),
    ]
    for n_str, col, mk, ls in styles:
        evens = sorted(int(k) for k in spec[n_str]["lambda_d2"].keys())
        lam0 = spec[n_str]["lambda_d2"]["0"]
        vals = [spec[n_str]["lambda_d2"][str(i)] / lam0 for i in evens]
        ax1.plot(evens, vals, marker=mk, linestyle=ls, color=col, linewidth=1.6, markersize=5.0, label=f"$n={n_str}$")

    ax1.set_yscale("log")
    ax1.set_xticks([0, 2, 4, 6, 8, 10])
    ax1.set_xlabel(r"Folded Johnson eigenspace index $i$ (even)")
    ax1.set_ylabel(r"Normalized eigenvalue $\Lambda_i(n,2)\,/\,\Lambda_0(n,2)$")
    ax1.set_title(r"(a) Spectral decay on $\overline{J}(2k,k)$ ($d=2$)")
    ax1.legend(frameon=False, loc="upper right")
    ax1.grid(axis="both", linestyle=":", alpha=0.35)

    # Panel (b): Parseval asymptotic prefactor convergence to (d+1)/sqrt(d)
    d_styles = [
        ("2", "#1f4e79", "o", "-"),
        ("3", "#2c8259", "s", "--"),
        ("4", "#c85a32", "^", "-."),
    ]
    for d_str, col, mk, ls in d_styles:
        rows = asymp[d_str]["rows"]
        ns = [r["n"] for r in rows]
        errs = [r["relative_error_to_limit"] for r in rows]
        ax2.plot(ns, errs, marker=mk, linestyle=ls, color=col, linewidth=1.6, markersize=5.0, label=f"$d={d_str}$")

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xticks([20, 40, 80, 160, 320])
    ax2.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax2.minorticks_off()
    ax2.set_xlabel(r"System size $n$ (log scale)")
    ax2.set_ylabel(r"Relative excess $(R_{n,d}/M_n)/L_{n,d} - 1$")
    ax2.set_title(r"(b) $O(1/n)$ convergence of $R_{n,d}/M_n$ to $\frac{d+1}{\sqrt{d}}(\frac{d+1}{2})^n$")
    ax2.legend(frameon=False, loc="upper right")
    ax2.grid(axis="both", linestyle=":", alpha=0.35)

    fig.savefig(FIG_DIR / "fig2_asymptotic_collapse.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig2_asymptotic_collapse.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    configure_style()
    data = json.loads(DATA_PATH.read_text(encoding="utf-8"))
    make_figure1(data)
    make_figure2(data)
    print("Generated fig1_es2_vs_haar_n12 and fig2_asymptotic_collapse cleanly.")


if __name__ == "__main__":
    main()
