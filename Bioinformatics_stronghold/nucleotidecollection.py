"""
Base class for biological sequence manipulation (DNA, RNA).
"""

from typing import Dict

from nucleotide import Nucleotide
from biosequencecollection import BioSequenceCollection

class NucleotideCollection(BioSequenceCollection):
    """Class for manipulating biological sequence collections (DNA, RNA)"""

    def __init__(self, seq_dict: Dict[str, Nucleotide] = None, seq_type: str = 'DNA'):

        if seq_dict is None:
            seq_dict = {'seq1': 'ATGC'}

        self.seq_type: str = seq_type.upper()

        self.seq_dict: Dict[str, Nucleotide] = {
            label: Nucleotide(seq, seq_type) for label, seq in seq_dict.items()
        }

    def gc_content_dict(self) -> Dict[str, float]:
        """Returns a dictionary of the GC content % of a dictionary of sequences"""
        return {key: nuc_seq.gc_content() for key, nuc_seq in self.seq_dict.items()}
