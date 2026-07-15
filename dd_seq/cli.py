"""Command-line entry points:
  dd_seq-fetch  UNIPROT -o out_dir
  dd_seq-align  out_dir --site-mode {pocket,ligand,none}
  dd_seq-run    UNIPROT -o out_dir --site-mode {pocket,ligand,none}
"""
from __future__ import annotations

import argparse

from . import pipeline
from .structalign import SITE_MODES


def build_fetch_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_seq-fetch",
        description="Download every known structure (all cross-referenced PDB entries + the AlphaFold DB model) of a protein, given its UniProt accession.",
    )
    parser.add_argument("uniprot", help="UniProt accession, e.g. P06493")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument("--no-progress", action="store_true", help="Suppress the one-line-per-item progress output")
    return parser


def _add_align_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--site-mode", choices=SITE_MODES, default="ligand",
        help="'ligand' (default): fit on residues near a bound ligand; 'pocket': fit on an fpocket-auto-detected "
             "druggable pocket; 'none': no active-site restriction, whole-chain CE structural alignment",
    )
    parser.add_argument("--reference", default=None, help="Label to superpose everything onto (default: AFDB, the AlphaFold model)")
    parser.add_argument("--site-source", default=None, help="Label of the structure the active site is defined on (default: auto -- see README)")
    parser.add_argument("--ligand-cutoff", type=float, default=5.0, help="Angstrom cutoff for --site-mode ligand (default: 5.0)")
    parser.add_argument("--pocket-rank", type=int, default=1, help="Druggability-ranked pocket to use for --site-mode pocket (default: 1, top-ranked)")
    parser.add_argument("--no-progress", action="store_true", help="Suppress the one-line-per-item progress output")


def build_align_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_seq-align",
        description="Align every fetched structure's sequence to the UniProt canonical sequence, and superpose every structure onto one reference by active site (or whole-chain).",
    )
    parser.add_argument("out_dir", help="Directory previously populated by dd_seq-fetch")
    _add_align_args(parser)
    return parser


def build_run_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_seq-run",
        description="dd_seq-fetch followed by dd_seq-align in one step.",
    )
    parser.add_argument("uniprot", help="UniProt accession, e.g. P06493")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    _add_align_args(parser)
    return parser


def main_fetch(argv=None) -> None:
    args = build_fetch_parser().parse_args(argv)
    manifest = pipeline.fetch_all(args.uniprot, args.out_dir, show_progress=not args.no_progress)
    n_pdb = sum(1 for e in manifest["entries"] if e["kind"] == "pdb")
    print(f"\n[done] {n_pdb} PDB entr{'y' if n_pdb == 1 else 'ies'} + 1 AlphaFold model -> {args.out_dir}")
    for s in manifest.get("skipped", []):
        print(f"  [skipped] {s['label']}: {s['reason']}")


def main_align(argv=None) -> None:
    args = build_align_parser().parse_args(argv)
    report = pipeline.analyze(
        args.out_dir, site_mode=args.site_mode, reference=args.reference, site_source=args.site_source,
        ligand_cutoff=args.ligand_cutoff, pocket_rank=args.pocket_rank, show_progress=not args.no_progress,
    )
    _print_report(report, args.out_dir)


def main_run(argv=None) -> None:
    args = build_run_parser().parse_args(argv)
    pipeline.fetch_all(args.uniprot, args.out_dir, show_progress=not args.no_progress)
    report = pipeline.analyze(
        args.out_dir, site_mode=args.site_mode, reference=args.reference, site_source=args.site_source,
        ligand_cutoff=args.ligand_cutoff, pocket_rank=args.pocket_rank, show_progress=not args.no_progress,
    )
    _print_report(report, args.out_dir)


def _print_report(report: dict, out_dir: str) -> None:
    print(f"\n[{report['uniprot_id']}] site_mode={report['site_mode']} site_source={report['site_source']} reference={report['reference']}")
    for s in report["structures"]:
        res = f"{s['resolution']:.2f}A" if s["resolution"] else "-"
        rmsd = f"{s['rmsd']:.3f} ({s['n_site_atoms']} atoms)" if s["rmsd"] is not None else f"SKIPPED: {s['align_error']}"
        print(f"  {s['label']:<8} {s['method']:<10} {res:<7} coverage={s['coverage']:.2f} mismatch={s['n_mismatch']:<3} rmsd={rmsd}")
    print(f"\n[done] report -> {out_dir}/report.json")
