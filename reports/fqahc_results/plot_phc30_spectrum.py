#!/usr/bin/env python3
"""Plot the saved 30-site PHC low-energy spectrum used by the optical note."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPORT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA = REPORT_DIR / "data" / "optical_response" / "phc30_spectrum_V1_100.dat"
DEFAULT_OUTPUT = (
    REPORT_DIR / "figures" / "optical_response" / "phc30_manybody_spectrum.png"
)


def load_spectrum(path: Path) -> np.ndarray:
    """Load and validate the two-column ``sector, E-E0`` spectrum file."""
    rows = np.loadtxt(path, comments="#", dtype=float)
    if rows.ndim != 2 or rows.shape[1] != 2:
        raise ValueError(f"expected two spectrum columns, got shape {rows.shape}")

    sectors = rows[:, 0]
    if not np.allclose(sectors, np.rint(sectors)):
        raise ValueError("momentum sectors must be integers")
    sector_ids = sectors.astype(int)
    if set(sector_ids) != set(range(15)):
        raise ValueError("expected momentum sectors 0 through 14")
    counts = np.bincount(sector_ids, minlength=15)
    if not np.all(counts == 8):
        raise ValueError(f"expected eight eigenvalues per sector, got {counts.tolist()}")
    return rows


def summarize_low_energy(
    rows: np.ndarray, manifold_size: int = 15
) -> dict[str, float | int]:
    """Return the width and isolation of the lowest ``manifold_size`` levels."""
    energies = np.sort(rows[:, 1])
    if len(energies) <= manifold_size:
        raise ValueError("spectrum does not contain a level above the manifold")
    ground = energies[0]
    manifold_top = energies[manifold_size - 1]
    next_level = energies[manifold_size]
    return {
        "manifold_size": manifold_size,
        "manifold_width": float(manifold_top - ground),
        "sixteenth_energy": float(next_level - ground),
        "separation": float(next_level - manifold_top),
    }


def format_summary_text(summary: dict[str, float | int]) -> str:
    """Format the low-energy annotation with real line breaks."""
    return (
        "15 states, one in every $k$ sector\n"
        f"$E_{{15}}-E_1={summary['manifold_width']:.6f}$\n"
        f"$E_{{16}}-E_{{15}}={summary['separation']:.6f}$"
    )


def make_figure(rows: np.ndarray, output: Path) -> None:
    """Create the full-spectrum and low-energy-zoom panels."""
    summary = summarize_low_energy(rows)
    order = np.argsort(rows[:, 1], kind="stable")
    manifold_mask = np.zeros(len(rows), dtype=bool)
    manifold_mask[order[: int(summary["manifold_size"])]] = True

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "legend.fontsize": 9,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.0), constrained_layout=True)

    for ax in axes:
        ax.scatter(
            rows[~manifold_mask, 0],
            rows[~manifold_mask, 1],
            s=24,
            color="#607d9b",
            alpha=0.78,
            linewidths=0,
            label="higher saved levels",
            zorder=2,
        )
        ax.scatter(
            rows[manifold_mask, 0],
            rows[manifold_mask, 1],
            s=48,
            color="#d62728",
            edgecolors="white",
            linewidths=0.5,
            label="lowest 15 states",
            zorder=3,
        )
        ax.set_xlim(-0.45, 14.45)
        ax.set_xticks(range(15))
        ax.set_xlabel("momentum sector $k$")
        ax.grid(True, color="0.85", linewidth=0.7)

    axes[0].set_title("all eight saved levels per sector")
    axes[0].set_ylabel(r"$E-E_0$")
    axes[0].set_ylim(-0.04, 2.70)
    axes[0].legend(loc="upper left", frameon=True)

    axes[1].set_title("low-energy zoom")
    axes[1].set_ylim(-0.008, 0.30)
    axes[1].axhline(
        summary["manifold_width"],
        color="#d62728",
        linestyle="--",
        linewidth=1.1,
        label=r"$E_{15}-E_1$",
        zorder=1,
    )
    axes[1].axhline(
        summary["sixteenth_energy"],
        color="#1f77b4",
        linestyle=":",
        linewidth=1.4,
        label=r"$E_{16}-E_1$",
        zorder=1,
    )
    axes[1].text(
        0.03,
        0.95,
        format_summary_text(summary),
        transform=axes[1].transAxes,
        va="top",
        ha="left",
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.92},
    )
    handles, labels = axes[1].get_legend_handles_labels()
    axes[1].legend(handles[-2:], labels[-2:], loc="center right", frameon=True)

    fig.suptitle(
        r"Tilted 30-site ED: $N_s=30$, $N_{\rm uc}=15$, $N_p=12$; "
        r"$t_1=1$, $t_3=0.2$, $V_1=100$, $V_2=V_3=0$",
        fontsize=13,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    rows = load_spectrum(args.data)
    summary = summarize_low_energy(rows)
    make_figure(rows, args.output)
    print(f"saved {args.output}")
    print(f"points={len(rows)} sectors={len(set(rows[:, 0].astype(int)))}")
    print(
        "manifold_size={manifold_size} manifold_width={manifold_width:.10f} "
        "sixteenth_energy={sixteenth_energy:.10f} separation={separation:.10f}".format(
            **summary
        )
    )


if __name__ == "__main__":
    main()
