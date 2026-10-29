"""Provenance: app version, data-source releases, methods and citations for every export.

So that student work is reproducible and citable, exports and saved items record *when*, *with which data
release* and *by which method* a number was produced. Check each source's own "how to cite" page for its latest
reference before publishing.
"""
from __future__ import annotations

import datetime as dt

from . import __version__

SOURCES = {
    "uniprot": {"name": "UniProtKB", "url": "https://www.uniprot.org", "licence": "CC BY 4.0",
                "cite": "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res. 2025;53:D609–D617.",
                "how_to_cite": "https://www.uniprot.org/help/publications"},
    "alphafold": {"name": "AlphaFold Protein Structure Database", "url": "https://alphafold.ebi.ac.uk", "licence": "CC BY 4.0",
                  "cite": "Jumper J, et al. Nature. 2021;596:583–589; Varadi M, et al. Nucleic Acids Res. 2024;52:D368–D375.",
                  "how_to_cite": "https://alphafold.ebi.ac.uk/faq"},
    "alphamissense": {"name": "AlphaMissense", "url": "https://github.com/google-deepmind/alphamissense", "licence": "CC BY 4.0",
                      "cite": "Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 2023;381:eadg7492.",
                      "how_to_cite": "https://github.com/google-deepmind/alphamissense"},
    "string": {"name": "STRING", "url": "https://string-db.org", "licence": "CC BY 4.0",
               "cite": "Szklarczyk D, et al. The STRING database in 2023. Nucleic Acids Res. 2023;51:D638–D646.",
               "how_to_cite": "https://string-db.org/cgi/about?footer_active_subpage=references"},
    "reactome": {"name": "Reactome", "url": "https://reactome.org", "licence": "CC BY 4.0",
                 "cite": "Milacic M, et al. The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res. 2024;52:D672–D678.",
                 "how_to_cite": "https://reactome.org/cite"},
    "interpro": {"name": "InterPro / Pfam", "url": "https://www.ebi.ac.uk/interpro", "licence": "CC0",
                 "cite": "Paysan-Lafosse T, et al. InterPro in 2022. Nucleic Acids Res. 2023;51:D418–D427.",
                 "how_to_cite": "https://www.ebi.ac.uk/interpro/about/interpro/"},
    "ebi_blast": {"name": "EMBL-EBI Job Dispatcher (NCBI BLAST+)", "url": "https://www.ebi.ac.uk/jdispatcher/sss/ncbiblast",
                  "licence": "EMBL-EBI terms of use",
                  "cite": "Madeira F, et al. The EMBL-EBI Job Dispatcher sequence analysis tools framework in 2024. Nucleic Acids Res. 2024;52:W521–W525; Camacho C, et al. BMC Bioinformatics. 2009;10:421.",
                  "how_to_cite": "https://www.ebi.ac.uk/jdispatcher/sss/ncbiblast"},
}

METHODS = {
    "molecular_weight": "Average isotopic residue masses + water (same values as ExPASy ProtParam; Gasteiger et al. 2005).",
    "pI": "Isoelectric point by bisection with Bjellqvist pK values, matching ExPASy Compute pI/Mw (Bjellqvist et al. 1993, Electrophoresis 14:1023).",
    "charge": "Net charge by Henderson–Hasselbalch with the same Bjellqvist pK values.",
    "gravy": "GRAVY: mean Kyte–Doolittle hydropathy (Kyte & Doolittle 1982, J Mol Biol 157:105).",
    "tm": "Transmembrane segments: Kyte–Doolittle, 19-residue window, threshold 1.6. A simple heuristic: use DeepTMHMM/Phobius for research.",
    "aromaticity": "Aromaticity: fraction of F+W+Y (Lobry & Gautier 1994, Nucleic Acids Res 22:3174).",
    "pairwise": "Pairwise alignment: dynamic programming with affine gaps (Gotoh 1982), BLOSUM62 (Henikoff & Henikoff 1992), gap open 11 / extend 1. Identity = identical columns / alignment length.",
    "msa": "Multiple alignment: center-star from global pairwise alignments. Tree: neighbour-joining (Saitou & Nei 1987) on Kimura-corrected protein distances, midpoint-rooted, no bootstrap support.",
    "secondary_structure": "Secondary structure: DSSP hydrogen-bond rules (Kabsch & Sander 1983) via a port of PyDSSP (Minami, MIT), 3 states, applied to AlphaFold models. Disorder: pLDDT < 50.",
    "superposition": "Superposition: residues paired by sequence alignment, Kabsch least-squares fit (Kabsch 1976) refined on pairs < 4 Å; TM-score formula of Zhang & Skolnick (2004) on those pairs (approximate; TM-align optimises the pairing).",
    "mutation": "Substitutions scored with BLOSUM62; pathogenicity from AlphaMissense (Cheng et al. 2023); known variants from UniProtKB.",
    "classifier": "Function classifier trained on reviewed UniProtKB entries (time split at 2018-12-31, taxonomically stratified). ESM-2 embeddings (Lin et al. 2023, Science 379:1123) + logistic regression, or 434 hand-crafted features + random forest. Predictions are not ground truth.",
    "explain": "Explanations: occlusion (region masking) and, for the random forest, exact tree-path feature contributions.",
    "blast": "BLAST: NCBI BLAST+ blastp run by EMBL-EBI (BLOSUM62, E-value ≤ 1e-3).",
    "digest": "Protease digest: simplified cleavage rules (e.g. trypsin after K/R not before P); monoisotopic residue masses; optional carbamidomethyl-Cys (+57.021464).",
    "a280": "Extinction coefficient at 280 nm: Pace et al. 1995 (Protein Sci 4:2411), as in ExPASy ProtParam.",
    "translation": "Translation: NCBI genetic code tables 1, 2 and 11.",
}

KIND_METHODS = {
    "analysis": ["molecular_weight", "pI", "charge", "gravy", "aromaticity", "tm"],
    "protein": ["molecular_weight", "pI", "charge", "gravy", "aromaticity", "tm", "secondary_structure"],
    "comparison": ["pairwise"], "alignment": ["msa", "pairwise"], "mutation": ["mutation", "pI", "gravy"],
    "classification": ["classifier", "explain"], "blast": ["blast"], "superposition": ["superposition", "secondary_structure"],
    "digest": ["digest"], "a280": ["a280"], "orf": ["translation"], "sequence": [], "note": [],
}
KIND_SOURCES = {
    "analysis": [], "protein": ["uniprot", "alphafold", "interpro", "string", "reactome"], "comparison": ["uniprot"],
    "alignment": ["uniprot"], "mutation": ["uniprot", "alphamissense"], "classification": ["uniprot"],
    "blast": ["ebi_blast", "uniprot"], "superposition": ["alphafold", "uniprot"], "digest": [], "a280": [], "orf": [],
    "sequence": [], "note": [],
}


def stamp(kind: str | None = None, uniprot_release: dict | None = None) -> dict:
    """Small provenance record attached to exports and saved items."""
    rel = uniprot_release or {}
    return {
        "app": f"Protein Analyzer {__version__}",
        "generated": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "uniprot_release": rel.get("release"), "uniprot_release_date": rel.get("date"),
        "methods": [METHODS[m] for m in KIND_METHODS.get(kind, [])] if kind else [],
        "sources": [f"{SOURCES[s]['name']} ({SOURCES[s]['licence']})" for s in KIND_SOURCES.get(kind, [])] if kind else [],
    }


def full() -> dict:
    return {"app": f"Protein Analyzer {__version__}", "version": __version__, "methods": METHODS, "sources": SOURCES,
            "kind_methods": KIND_METHODS, "kind_sources": KIND_SOURCES}
