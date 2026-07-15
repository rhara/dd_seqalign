"""Sequence-coverage track visualization: one row per structure, one
column per canonical UniProt position, colored by whether that position is
resolved and matches canonical in that structure. Built from each
structure's `coverage_string` (see `pipeline._coverage_string`) rather than
re-deriving it from the full per-residue alignment, since that compact
string is already exactly the character-per-position representation this
needs.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

STATUS_MISSING = 0
STATUS_MATCH = 1
STATUS_MISMATCH = 2

_CMAP_COLORS = ["#e8e8e8", "#4c78a8", "#d62728"]  # missing / match / mismatch


def coverage_matrix(structures: Sequence[dict]) -> np.ndarray:
    """`(n_structures, canonical_length)` int array of STATUS_* codes,
    read off each structure dict's `coverage_string`."""
    n = len(structures[0]["coverage_string"]) if structures else 0
    mat = np.zeros((len(structures), n), dtype=int)
    for i, s in enumerate(structures):
        for j, ch in enumerate(s["coverage_string"]):
            mat[i, j] = STATUS_MISSING if ch == "-" else (STATUS_MATCH if ch == "." else STATUS_MISMATCH)
    return mat


def plot_coverage(
    structures: Sequence[dict], *, site_canonical_positions: Optional[Sequence[int]] = None, fig_width: float = 14.0,
):
    """A matplotlib figure: one horizontal row per structure (label on the
    y-axis), colored by residue status at each canonical position;
    `site_canonical_positions` (1-indexed) are marked with thin vertical
    lines across all rows so the active site's location relative to each
    structure's coverage gaps is visible at a glance."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    mat = coverage_matrix(structures)
    labels = [s["label"] for s in structures]

    fig_height = max(2.0, 0.35 * len(structures) + 1.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.imshow(mat, aspect="auto", cmap=ListedColormap(_CMAP_COLORS), vmin=0, vmax=2, interpolation="nearest")

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Canonical UniProt residue position")

    if site_canonical_positions:
        for pos in site_canonical_positions:
            ax.axvline(pos - 1, color="black", linewidth=0.6, alpha=0.5)

    from matplotlib.patches import Patch

    legend_items = [
        Patch(facecolor=_CMAP_COLORS[STATUS_MATCH], label="match"),
        Patch(facecolor=_CMAP_COLORS[STATUS_MISMATCH], label="mismatch"),
        Patch(facecolor=_CMAP_COLORS[STATUS_MISSING], label="not resolved"),
    ]
    if site_canonical_positions:
        legend_items.append(Patch(facecolor="none", edgecolor="black", label="active site position"))
    ax.legend(handles=legend_items, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=4, frameon=False)

    fig.tight_layout()
    return fig
