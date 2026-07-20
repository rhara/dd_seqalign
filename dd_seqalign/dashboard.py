"""Summary table (one row per structure) from a dd_seqalign report.json --
the tabular counterpart to `scene.py`'s 3D overlay and `seqplot.py`'s
coverage track.
"""
from __future__ import annotations

from typing import Dict

import pandas as pd


def summary_dataframe(report: Dict) -> pd.DataFrame:
    rows = []
    for s in report["structures"]:
        rows.append(
            {
                "label": s["label"],
                "chain": s["chain"],
                "method": s["method"],
                "resolution": s["resolution"],
                "coverage": round(s["coverage"], 3),
                "mismatches": s["n_mismatch"],
                "rmsd": round(s["rmsd"], 3) if s["rmsd"] is not None else None,
                "site_atoms": s["n_site_atoms"],
                "title": s["title"],
            }
        )
    return pd.DataFrame(rows).sort_values("label").reset_index(drop=True)
