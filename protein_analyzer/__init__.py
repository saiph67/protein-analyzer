"""Protein Analyzer: sequence analysis, UniProt access, comparison and function prediction."""
__version__ = "1.5.0"

from .analyzer import analyze, read_fasta, SequenceError
from .compare import align, kmer_similarity, nearest_neighbor_similarity
from .uniprot import UniProtClient, parse_entry, save_fasta

__all__ = ["analyze", "read_fasta", "SequenceError", "align", "kmer_similarity",
           "nearest_neighbor_similarity", "UniProtClient", "parse_entry", "save_fasta"]
