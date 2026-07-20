"""Discover and download every known structure of a protein: all PDB
entries cross-referenced to a UniProt accession (RCSB Search API) plus its
AlphaFold DB predicted model, and the UniProt canonical sequence itself
(the reference every structure's own sequence gets compared against in
`sequence.py`). Actual structure file download is not reimplemented here --
`dd_prep.fetch.download_pdb`/`download_afdb` already do that (with local-file
caching), so `pipeline.py` calls those directly; this module only adds the
"which PDB entries exist for this UniProt accession" lookup and PDB entries'
experimental metadata, neither of which dd_prep needs for its own job of
prepping a single already-chosen structure.
"""
from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass
from typing import List, Optional

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"


def fetch_uniprot_fasta(accession: str) -> str:
    """The canonical (isoform 1) amino-acid sequence for a UniProt
    accession, as one contiguous string (FASTA header line stripped)."""
    with urllib.request.urlopen(UNIPROT_FASTA.format(accession=accession.upper())) as fh:
        text = fh.read().decode()
    lines = text.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"UniProt has no entry for {accession!r}")
    return "".join(lines[1:])


def list_pdb_ids_for_uniprot(accession: str) -> List[str]:
    """Every RCSB PDB entry ID cross-referenced (via SIFTS) to this UniProt
    accession, covering X-ray, EM, and any other experimental method --
    RCSB's search index doesn't distinguish structure determination method
    here, so `fetch_entry_metadata` is what tells them apart afterwards."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession.upper(),
            },
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    req = urllib.request.Request(
        RCSB_SEARCH, data=json.dumps(query).encode(), headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req) as fh:
        result = json.load(fh)
    return [hit["identifier"] for hit in result.get("result_set", [])]


@dataclass
class EntryMetadata:
    pdb_id: str
    method: str
    resolution: Optional[float]
    title: str


def fetch_entry_metadata(pdb_id: str) -> EntryMetadata:
    """Experimental method (X-ray/EM/NMR/...), resolution (None for
    methods that don't report one, e.g. NMR), and title, straight from
    RCSB's entry-level summary. Chain/ligand composition is *not* fetched
    here -- once the structure file itself is downloaded, `sequence.py`/
    `activesite.py` read that directly instead of trusting a second,
    possibly-inconsistent API description of the same file."""
    with urllib.request.urlopen(RCSB_ENTRY.format(pdb_id=pdb_id.upper())) as fh:
        entry = json.load(fh)
    info = entry.get("rcsb_entry_info", {})
    resolution_list = info.get("resolution_combined") or []
    return EntryMetadata(
        pdb_id=pdb_id.upper(),
        method=info.get("experimental_method", "unknown"),
        resolution=resolution_list[0] if resolution_list else None,
        title=entry.get("struct", {}).get("title", ""),
    )
