"""End-to-end orchestration in two phases, mirroring dd_prep/dd_af's own
fetch-then-process split:

- `fetch_all`: discover and download every known structure of a protein
  (every cross-referenced PDB entry, plus the AlphaFold DB model) and the
  UniProt canonical sequence, writing `manifest.json`.
- `analyze`: for every downloaded structure, align its sequence to the
  canonical sequence (`sequence.py`), define an active site once on one
  "site source" structure and map it into every other structure's own
  numbering (`activesite.py`), then superpose everything onto one
  reference structure (`structalign.py`), writing `report.json` plus the
  superposed coordinate files under `aligned/`.
"""
from __future__ import annotations

import json
import urllib.error
from pathlib import Path
from typing import Dict, List, Optional, Union

import dd_prep.fetch as dd_prep_fetch

from . import activesite, structalign
from .fetch import fetch_entry_metadata, fetch_uniprot_fasta, list_pdb_ids_for_uniprot
from .sequence import align_to_canonical, extract_chain_sequences, pick_target_chain

AFDB_LABEL = "AFDB"


def fetch_all(uniprot_id: str, out_dir: Union[str, Path]) -> dict:
    """Download every RCSB entry cross-referenced to `uniprot_id`, its
    AlphaFold DB model, and the UniProt canonical sequence. Cached like
    `dd_prep.fetch`: re-running with the same `out_dir` skips files
    already on disk. `out_dir` is resolved to an absolute path before any
    path gets written into manifest.json/report.json, so those stay valid
    regardless of which directory a later process (e.g. the Streamlit app)
    happens to be run from."""
    out_dir = Path(out_dir).resolve()
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    canonical = fetch_uniprot_fasta(uniprot_id)
    (out_dir / "canonical.fasta").write_text(f">{uniprot_id}\n{canonical}\n")

    entries: List[dict] = []
    skipped: List[dict] = []
    for pdb_id in list_pdb_ids_for_uniprot(uniprot_id):
        meta = fetch_entry_metadata(pdb_id)
        dest = raw_dir / f"{pdb_id}_raw.pdb"
        try:
            dd_prep_fetch.download_pdb(pdb_id, dest)
        except urllib.error.HTTPError as e:
            # A handful of very recently released entries have no legacy
            # .pdb file generated yet (mmCIF-only) -- skip rather than
            # abort the whole batch fetch over one entry.
            skipped.append({"label": pdb_id, "reason": f"HTTP {e.code} fetching legacy .pdb format"})
            continue
        entries.append(
            {
                "label": pdb_id, "kind": "pdb", "path": str(dest),
                "method": meta.method, "resolution": meta.resolution, "title": meta.title,
            }
        )

    afdb_dest = raw_dir / f"{AFDB_LABEL}_raw.pdb"
    dd_prep_fetch.download_afdb(uniprot_id, afdb_dest)
    entries.append(
        {
            "label": AFDB_LABEL, "kind": "afdb", "path": str(afdb_dest),
            "method": "AlphaFold", "resolution": None, "title": "AlphaFold DB predicted model",
        }
    )

    manifest = {"uniprot_id": uniprot_id.upper(), "canonical_length": len(canonical), "entries": entries, "skipped": skipped}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def _read_canonical(out_dir: Path) -> str:
    lines = (out_dir / "canonical.fasta").read_text().splitlines()
    return "".join(lines[1:])


def _coverage_string(aln) -> str:
    """One character per canonical position: '.' = matches canonical,
    '-' = not resolved in this structure, else the observed (mismatching)
    one-letter code -- compact enough to embed directly in report.json and
    drive a coverage-track visualization straight off the string."""
    return "".join("-" if r.status == "missing" else ("." if r.status == "match" else r.structure_code) for r in aln.residues)


def analyze(
    out_dir: Union[str, Path], *,
    site_mode: str = "ligand",
    reference: Optional[str] = None,
    site_source: Optional[str] = None,
    ligand_cutoff: float = 5.0,
    pocket_rank: int = 1,
) -> dict:
    """Run sequence alignment + active-site structural alignment across
    every structure `fetch_all` downloaded into `out_dir`.

    `reference` defaults to the AlphaFold model (`AFDB_LABEL`): it's the
    one structure guaranteed to be full-length and ligand-free, so it's a
    safe universal superposition target regardless of `site_mode` --the
    active site itself, when needed, still comes from `site_source` (a
    real co-crystal for `site_mode="ligand"`), not from the reference.
    `site_source` defaults to `reference` for `"pocket"` (fpocket needs no
    ligand) or to the first ligand-bearing entry found for `"ligand"`.
    """
    out_dir = Path(out_dir).resolve()
    manifest = json.loads((out_dir / "manifest.json").read_text())
    canonical = _read_canonical(out_dir)
    reference = reference or AFDB_LABEL

    per_structure: Dict[str, dict] = {}
    for entry in manifest["entries"]:
        chains = extract_chain_sequences(entry["path"])
        alignments = {cid: align_to_canonical(cs, canonical) for cid, cs in chains.items()}
        target_chain = pick_target_chain(alignments)
        aln = alignments[target_chain]
        per_structure[entry["label"]] = {
            "path": entry["path"], "chain": target_chain, "alignment": aln,
            "coverage": aln.coverage, "n_mismatch": aln.n_mismatch, "coverage_string": _coverage_string(aln),
            "method": entry["method"], "resolution": entry["resolution"], "title": entry["title"],
        }

    if reference not in per_structure:
        raise ValueError(f"reference {reference!r} not found among fetched structures: {sorted(per_structure)}")

    site_canonical: Optional[List[int]] = None
    resolved_site_source = None
    if site_mode == "pocket":
        resolved_site_source = site_source or reference
        src = per_structure[resolved_site_source]
        raw_site = activesite.site_from_pocket(src["path"], chain_id=src["chain"], pocket_rank=pocket_rank)
        site_canonical = sorted({src["alignment"].canonical_for_resseq(resseq) for _c, resseq in raw_site} - {None})
    elif site_mode == "ligand":
        candidates = [site_source] if site_source else list(per_structure)
        for label in candidates:
            src = per_structure[label]
            raw_site = activesite.site_from_ligand(src["path"], chain_id=src["chain"], cutoff=ligand_cutoff)
            if raw_site:
                resolved_site_source = label
                site_canonical = sorted({src["alignment"].canonical_for_resseq(resseq) for _c, resseq in raw_site} - {None})
                break
        if site_canonical is None:
            raise ValueError(
                "site_mode='ligand': no structure with a bound ligand was found "
                "(pass --site-source explicitly, or use --site-mode pocket instead)"
            )
    elif site_mode != "none":
        raise ValueError(f"site_mode must be one of {structalign.SITE_MODES}, got {site_mode!r}")

    structures = []
    site_resseqs_by_label: Dict[str, Optional[List[int]]] = {}
    for label, info in per_structure.items():
        if site_canonical is not None:
            resseqs = sorted({info["alignment"].resseq_for_canonical(pos) for pos in site_canonical} - {None})
        else:
            resseqs = None
        site_resseqs_by_label[label] = resseqs
        structures.append(
            structalign.StructureInput(label=label, pdb_path=info["path"], chain_id=info["chain"], site_resseqs=resseqs)
        )

    align_results = {r.label: r for r in structalign.align_structures(
        structures, reference_label=reference, out_dir=out_dir / "aligned", site_mode=site_mode,
    )}

    report = {
        "uniprot_id": manifest["uniprot_id"],
        "canonical_length": manifest["canonical_length"],
        "reference": reference,
        "site_mode": site_mode,
        "site_source": resolved_site_source,
        "site_canonical_positions": site_canonical,
        "structures": [
            {
                "label": label,
                "chain": info["chain"],
                "method": info["method"],
                "resolution": info["resolution"],
                "title": info["title"],
                "coverage": info["coverage"],
                "n_mismatch": info["n_mismatch"],
                "coverage_string": info["coverage_string"],
                "rmsd": align_results[label].rmsd,
                "n_site_atoms": align_results[label].n_atoms,
                "aligned_pdb": align_results[label].aligned_pdb,
                "align_error": align_results[label].error,
                "site_resseqs": site_resseqs_by_label[label],
            }
            for label, info in per_structure.items()
        ],
    }
    (out_dir / "report.json").write_text(json.dumps(report, indent=2))
    return report
