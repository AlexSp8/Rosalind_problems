"""
Base class for biological sequence manipulation (DNA, RNA).
"""

from typing import List, Dict

from biosequence import BioSequence

class BioSequenceCollection:
    """Class for manipulating biological sequence collections (DNA, RNA)"""

    def __init__(self, seq_dict: Dict[str, str] = None):
        if seq_dict is None:
            seq_dict = {'seq1': 'ATGC'}

        self.seq_dict: Dict[str, BioSequence] = {
            label: BioSequence(seq) for label, seq in seq_dict.items()
        }

    def monomer_content_dict(self, m: str) -> Dict[str, float]:
        """Returns a dictionary of the monomer (m) content % of a dictionary of sequences"""
        return {key: bio_seq.monomer_content(m) for key, bio_seq in self.seq_dict.items()}

    def profile_matrix(self, monomers: List[str] = None) -> List[List[int]]:
        """Returns the profile matrix of a dictionary of sequences
        given the monomers alphabet of the sequences"""

        bio_seq_list = list(self.seq_dict.values())

        monomer_index = {n: i for i, n in enumerate(monomers)}

        rows = len(monomers)
        cols = max(len(bio_seq.seq) for bio_seq in bio_seq_list)

        profile_matrix = [[0]*cols for _ in range(rows)]
        for bio_seq in bio_seq_list:
            for j, m in enumerate(bio_seq.seq):
                profile_matrix[monomer_index[m]][j] += 1
        return profile_matrix

    @staticmethod
    def consensus_string(profile_matrix: List[List[int]], monomers: List[str]) -> str:
        """Returns the consensus string of a profile matrix
        given the monomers alphabet of the sequences"""

        rows = len(monomers)
        cols = len(profile_matrix[0])

        consensus = []
        for j in range(cols):
            row = max(range(rows), key=lambda i: profile_matrix[i][j])
            consensus.append(monomers[row])
        return "".join(consensus)

    def adjacency_list(self, k: int) -> List[str]:
        """Returns the list of k-th order adjacent edges of a dictionary of sequences"""

        suffix = {key: bio_seq.seq[-k:] for key, bio_seq in self.seq_dict.items()}
        prefix = {key: bio_seq.seq[:k] for key, bio_seq in self.seq_dict.items()}

        edges = []
        for s in self.seq_dict.keys():
            for t in self.seq_dict.keys():
                if s != t and suffix[s] == prefix[t]:
                    edges.append(f"{s[1:]} {t[1:]}")
        return edges

    def longest_shared_motif(self) -> str:
        """Returns the longest common substring of a dictionary of sequences"""

        seq_list = [bio_seq.seq for bio_seq in self.seq_dict.values()]

        shortest_seq = min(seq_list, key=len)
        max_l = len(shortest_seq)
        for l in range(max_l, 0, -1):
            for i in range(max_l - l + 1):
                motif = shortest_seq[i:i + l]
                if all(motif in seq for seq in seq_list):
                    return motif
