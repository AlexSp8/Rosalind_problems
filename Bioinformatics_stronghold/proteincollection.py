"""
Base class for biological sequence manipulation (DNA, RNA).
"""

from typing import List, Dict

from biosequencecollection import BioSequenceCollection
from protein import Protein

class ProteinCollection(BioSequenceCollection):
    """Class for manipulating biological sequence collections (DNA, RNA)"""

    def __init__(self, seq_dict: Dict[str, Protein] = None):
        if seq_dict is None:
            seq_dict = {'seq1': 'M'}

        self.seq_dict: Dict[str, Protein] = {
            label: Protein(seq) for label, seq in seq_dict.items()
        }

    def motif_locations_dict(self, motif_rules: List[tuple] = None) -> Dict[str, List[int]]:
        """Returns a dictionary of locations in a dictionary of sequence
        where a motif's rules are satisfied
        """
        locations_dict = {}
        for key, prot_seq in self.seq_dict.items():
            motif_locations_list = prot_seq.motif_locations(motif_rules)
            if motif_locations_list:
                locations_dict[key] = motif_locations_list

        return locations_dict

    @staticmethod
    def all_proteins_from_all_reading_frames(
        reading_frames: List[str]) -> List[Protein]:
        """Returns all protein sequences from a list of reading frames"""

        proteins_dict = {}
        for rf in reading_frames:
            proteins = ProteinCollection.all_proteins_from_reading_frame(rf)
            for p in proteins:
                proteins_dict[p.seq] = p
        return list(proteins_dict.values())

    @staticmethod
    def all_proteins_from_reading_frame(rf: str) -> List[Protein]:
        """Returns all protein sequences from a single reading frame"""

        proteins = []
        for i in range(len(rf)):
            protein = Protein.first_protein_from_reading_frame(rf[i:])
            if protein:
                proteins.append(protein)
        return proteins
