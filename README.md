# dd_seq

Compares every known structure of a protein -- every PDB entry cross-
referenced to a UniProt accession (X-ray, EM, any oligomeric state or
fragment), plus the AlphaFold DB predicted model -- on two axes: sequence
coverage against the canonical UniProt sequence, and active-site-based
structural (RMSD) alignment. Designed as a reusable package, not tied to
any specific target (same philosophy as `dd_prep`/`dd_af`/`dd_viewer`/etc.
-- every example below uses human CDK1, UniProt `P06493`, but any
accession works). Reuses `dd_prep` (structure download, HETATM
classification) and `dd_af` (fpocket-based pocket detection) directly
rather than reimplementing either.

- **Fetch (`dd_seq-fetch`)**: `list_pdb_ids_for_uniprot` (RCSB Search API)
  finds every PDB entry cross-referenced to the accession; each is
  downloaded via `dd_prep.fetch.download_pdb`, plus the AlphaFold DB model
  via `dd_prep.fetch.download_afdb` and the canonical sequence via the
  UniProt REST API. Re-running against the same `-o` directory skips
  anything already on disk (canonical.fasta, each PDB entry, the AlphaFold
  model) rather than re-downloading it, printing `already downloaded,
  skipping` for each -- `list_pdb_ids_for_uniprot` itself is still queried
  fresh every run, so a newly-released entry gets picked up on a re-run
  without re-fetching everything else. A handful of very recently released
  entries have no legacy `.pdb` file yet (mmCIF-only) -- these are skipped,
  not fatal, and recorded in `manifest.json`'s `"skipped"` list.
- **Align (`dd_seq-align`)**: for every fetched structure, extracts each
  chain's sequence (`sequence.py`, Biopython) and glocal-aligns it against
  the canonical sequence (free end gaps -- every input is a fragment/
  isoform of the *same* protein, not a set of divergent homologs, so a
  full MSA tool is unnecessary; a single reference is enough). The chain
  actually corresponding to the protein of interest is picked by identity
  (`pick_target_chain`, ranked by matching-residue count, not raw
  coverage -- necessary to avoid picking a homologous partner chain, e.g.
  CDK7 in a CAK-CDK1-cyclinB1 assembly, which can have *higher* coverage
  than the true target chain but is mostly mismatches).

  An active site is then defined once on one "site source" structure
  (`activesite.py`, two modes -- `--site-mode ligand`: residues near the
  auto-picked bound ligand; `--site-mode pocket`: fpocket's top-ranked
  druggable pocket via `dd_af.pocket`) and translated into every other
  structure's own residue numbering by round-tripping through canonical
  UniProt positions (`map_site_to_structure`) -- this is what makes the
  site comparable across structures with completely different numbering/
  chain layouts. Every structure is then superposed onto one reference
  (default: the AlphaFold model, since it's always full-length and
  ligand-free) via PyMOL (`structalign.py`): `cmd.pair_fit` on the known
  site-residue correspondence for `ligand`/`pocket` mode, or `cmd.cealign`
  (topology-independent CE structural alignment, no residue
  correspondence needed) for `--site-mode none`. A structure that doesn't
  resolve the site at all (e.g. a co-complex crystallized around an
  unrelated fragment of the protein, not its folded domain) is skipped
  with a recorded reason rather than aborting the whole batch.
- **Run (`dd_seq-run`)**: fetch + align in one step.
- **App (`streamlit run app.py -- --report-dir DIR`)**: three tabs --
  Overview (per-structure method/resolution/coverage/RMSD table),
  Sequence coverage (a match/mismatch/not-resolved track per structure
  across canonical positions), Structure overlay (py3Dmol, every
  structure's target chain superposed and colored distinctly, active site
  highlighted, ligands optional).

## Installation

Requires Biopython, pandas, numpy, PyMOL (`pymol2`, importable as a
library -- not the GUI), the `fpocket` CLI, and the `dd_prep`/`dd_af`
packages. The `mpro` conda env already has everything:

```bash
cd dd_prep && pip install -e . && cd ..   # if not already installed
cd dd_af && pip install -e . && cd ..     # if not already installed
cd dd_seq && pip install -e ".[app]"      # [app] adds streamlit/py3Dmol/matplotlib
```

This installs three console commands: `dd_seq-fetch`, `dd_seq-align`,
`dd_seq-run`.

## Usage

```bash
dd_seq-run P06493 -o data --site-mode ligand
streamlit run app.py -- --report-dir data
```

`--site-mode` (default `ligand`): `ligand` (fit on residues near a bound
ligand), `pocket` (fit on an fpocket-auto-detected druggable pocket, works
on apo/AlphaFold structures too), or `none` (no active-site restriction,
whole-chain CE alignment). `--reference`/`--site-source` override the
defaults described above; `--ligand-cutoff`/`--pocket-rank` tune site
detection.

All three commands print one line per completed item as it happens
(fetch/skip per structure, sequence-alignment result per structure,
structural-fit result or skip reason per structure) -- pass
`--no-progress` to suppress this and only print the final summary table.

## Design notes

- **Why not a real MSA tool**: mafft/clustalo aren't in the `mpro` env,
  and aren't the right tool anyway -- every structure here is the same
  protein, so a reference-based pairwise glocal alignment (Biopython
  `PairwiseAligner`, BLOSUM62, free end gaps) against the UniProt
  canonical sequence gives a more directly useful result (a per-canonical-
  position coverage/mismatch table across every structure at once) than a
  generic multiple sequence alignment would.
- **Canonical UniProt position as the common coordinate system**: every
  structure has its own author residue numbering (offset, insertion
  codes, gaps from missing density); rather than trying to reconcile
  those numbering schemes pairwise, everything (active-site residues,
  coverage tracks) is expressed in canonical UniProt positions and
  translated into a given structure's own numbering only at the point of
  use (`ChainAlignment.resseq_for_canonical`/`canonical_for_resseq`).
- **`pair_fit` over `align`/`cealign` for site-mode fitting**: the site
  residue correspondence is already known exactly (both sides are the
  same canonical positions), so `cmd.pair_fit` (direct Kabsch
  superposition on given atom pairs) is used instead of `cealign`/`align`
  (which do their own internal structural/sequence re-matching) -- no
  risk of PyMOL silently pairing the wrong residues.

## Known limitations

- A co-crystal where the protein of interest contributes only a small
  unrelated peptide fragment (not its folded domain bound to a partner
  protein's own site) will have no chain that meaningfully covers the
  active site region -- `dd_seq-align` correctly skips these (see
  `report.json`'s `"align_error"` per structure) rather than fitting them
  incorrectly.
- `site_from_pocket`/fpocket needs a single, isolated chain to detect a
  sensible pocket; it is run on the target chain stripped of every other
  chain, so an inter-chain-only pocket (e.g. a groove that only exists at
  a protein-protein interface) will not be found this way.
