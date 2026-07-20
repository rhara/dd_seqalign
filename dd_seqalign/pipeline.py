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


def fetch_all(uniprot_id: str, out_dir: Union[str, Path], *, show_progress: bool = True) -> dict:
    """Download every RCSB entry cross-referenced to `uniprot_id`, its
    AlphaFold DB model, and the UniProt canonical sequence. Cached like
    `dd_prep.fetch`: re-running with the same `out_dir` skips files
    already on disk. `out_dir` is resolved to an absolute path before any
    path gets written into manifest.json/report.json, so those stay valid
    regardless of which directory a later process (e.g. the Streamlit app)
    happens to be run from.

    `show_progress` prints one line per completed item as it happens
    (`print(..., flush=True)`, same convention as dd_prep/dd_af) -- there
    is no `--n-jobs` here (fetches are sequential, one entry at a time),
    so unlike those sibling projects' `parallel_map`-driven CLI loops, the
    progress printing lives directly in this function rather than in
    `cli.py`.
    """
    out_dir = Path(out_dir).resolve()
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    canonical_dest = out_dir / "canonical.fasta"
    if canonical_dest.exists():
        canonical = _read_canonical(out_dir)
        if show_progress:
            print(f"[fetch] {uniprot_id.upper()}: canonical sequence already downloaded, skipping -> {canonical_dest.name}", flush=True)
    else:
        canonical = fetch_uniprot_fasta(uniprot_id)
        canonical_dest.write_text(f">{uniprot_id}\n{canonical}\n")
        if show_progress:
            print(f"[fetch] {uniprot_id.upper()}: canonical sequence ({len(canonical)} aa) -> {canonical_dest.name}", flush=True)

    pdb_ids = list_pdb_ids_for_uniprot(uniprot_id)
    if show_progress:
        print(f"[fetch] {uniprot_id.upper()}: {len(pdb_ids)} PDB entr{'y' if len(pdb_ids) == 1 else 'ies'} found on RCSB", flush=True)

    entries: List[dict] = []
    skipped: List[dict] = []
    for i, pdb_id in enumerate(pdb_ids, start=1):
        meta = fetch_entry_metadata(pdb_id)
        dest = raw_dir / f"{pdb_id}_raw.pdb"
        already_had_it = dest.exists()
        if not already_had_it:
            try:
                dd_prep_fetch.download_pdb(pdb_id, dest)
            except urllib.error.HTTPError as e:
                # A handful of very recently released entries have no legacy
                # .pdb file generated yet (mmCIF-only) -- skip rather than
                # abort the whole batch fetch over one entry.
                reason = f"HTTP {e.code} fetching legacy .pdb format"
                skipped.append({"label": pdb_id, "reason": reason})
                if show_progress:
                    print(f"[fetch] ({i}/{len(pdb_ids)}) {pdb_id}: skipped ({reason})", flush=True)
                continue
        entries.append(
            {
                "label": pdb_id, "kind": "pdb", "path": str(dest),
                "method": meta.method, "resolution": meta.resolution, "title": meta.title,
            }
        )
        if show_progress:
            if already_had_it:
                print(f"[fetch] ({i}/{len(pdb_ids)}) {pdb_id}: already downloaded, skipping -> {dest.name}", flush=True)
            else:
                res = f", {meta.resolution:.2f}A" if meta.resolution else ""
                print(f"[fetch] ({i}/{len(pdb_ids)}) {pdb_id}: {meta.method}{res} -> {dest.name}", flush=True)

    afdb_dest = raw_dir / f"{AFDB_LABEL}_raw.pdb"
    afdb_already_had_it = afdb_dest.exists()
    dd_prep_fetch.download_afdb(uniprot_id, afdb_dest)
    entries.append(
        {
            "label": AFDB_LABEL, "kind": "afdb", "path": str(afdb_dest),
            "method": "AlphaFold", "resolution": None, "title": "AlphaFold DB predicted model",
        }
    )
    if show_progress:
        if afdb_already_had_it:
            print(f"[fetch] {AFDB_LABEL}: already downloaded, skipping -> {afdb_dest.name}", flush=True)
        else:
            print(f"[fetch] {AFDB_LABEL}: AlphaFold DB model -> {afdb_dest.name}", flush=True)

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
    show_progress: bool = True,
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

    `show_progress` prints one line per completed item as it happens (see
    `fetch_all`'s docstring for why this lives here rather than in
    `cli.py`).
    """
    out_dir = Path(out_dir).resolve()
    manifest = json.loads((out_dir / "manifest.json").read_text())
    canonical = _read_canonical(out_dir)
    reference = reference or AFDB_LABEL

    if show_progress:
        print(f"[align] sequence-aligning {len(manifest['entries'])} structure(s) against the canonical sequence...", flush=True)

    per_structure: Dict[str, dict] = {}
    for i, entry in enumerate(manifest["entries"], start=1):
        chains = extract_chain_sequences(entry["path"])
        alignments = {cid: align_to_canonical(cs, canonical) for cid, cs in chains.items()}
        target_chain = pick_target_chain(alignments)
        aln = alignments[target_chain]
        per_structure[entry["label"]] = {
            "path": entry["path"], "chain": target_chain, "alignment": aln,
            "coverage": aln.coverage, "n_mismatch": aln.n_mismatch, "coverage_string": _coverage_string(aln),
            "method": entry["method"], "resolution": entry["resolution"], "title": entry["title"],
        }
        if show_progress:
            print(
                f"[align] ({i}/{len(manifest['entries'])}) {entry['label']}: chain {target_chain}, "
                f"coverage={aln.coverage:.2f}, mismatches={aln.n_mismatch}", flush=True,
            )

    if reference not in per_structure:
        raise ValueError(f"reference {reference!r} not found among fetched structures: {sorted(per_structure)}")

    site_canonical: Optional[List[int]] = None
    resolved_site_source = None
    if site_mode == "pocket":
        resolved_site_source = site_source or reference
        if show_progress:
            print(f"[site] site_mode=pocket: running fpocket on {resolved_site_source}...", flush=True)
        src = per_structure[resolved_site_source]
        raw_site = activesite.site_from_pocket(src["path"], chain_id=src["chain"], pocket_rank=pocket_rank)
        site_canonical = sorted({src["alignment"].canonical_for_resseq(resseq) for _c, resseq in raw_site} - {None})
        if show_progress:
            print(f"[site] pocket rank {pocket_rank} on {resolved_site_source}: {len(site_canonical)} residue(s)", flush=True)
    elif site_mode == "ligand":
        candidates = [site_source] if site_source else list(per_structure)
        for label in candidates:
            src = per_structure[label]
            if show_progress:
                print(f"[site] site_mode=ligand: checking {label} for a bound ligand...", flush=True)
            raw_site = activesite.site_from_ligand(src["path"], chain_id=src["chain"], cutoff=ligand_cutoff)
            if raw_site:
                resolved_site_source = label
                site_canonical = sorted({src["alignment"].canonical_for_resseq(resseq) for _c, resseq in raw_site} - {None})
                if show_progress:
                    print(f"[site] ligand-proximal site on {label}: {len(site_canonical)} residue(s)", flush=True)
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
        show_progress=show_progress,
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
