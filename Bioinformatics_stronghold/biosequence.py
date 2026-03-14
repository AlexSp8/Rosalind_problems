"""
Base class for biological sequence manipulation (DNA, RNA).
"""

from typing import List, Dict

from collections import Counter
from itertools import product


class BioSequence:
    """Class for manipulating biological sequences"""

    def __init__(self, seq: str = 'ATGC'):
        self.seq: str = seq.upper()

    def monomers_count(self) -> Dict[str, int]:
        """Returns the count of each monomer in a sequence"""
        return Counter(self.seq)

    def reverse(self) -> BioSequence:
        """Returns the reverse of a sequence"""
        return BioSequence(self.seq[::-1], self.seq_type)

    def monomer_content(self, m: str) -> float:
        """Returns the monomer (m) content % of a sequence"""
        if not self.seq:
            return 0.0
        content = self.seq.count(m)
        return content*100/len(self.seq)

    def point_mutations(self, mutation_seq: str) -> int:
        """Returns the number of point mutations between two sequences"""
        if len(self.seq) != len(mutation_seq):
            raise ValueError("Sequences must be the same length!")
        return sum(c1 != c2 for c1, c2 in zip(self.seq, mutation_seq))

    @staticmethod
    def motif_rules(motif: str) -> List[tuple]:
        """Returns a list of rules for a motif"""

        rules = []
        i = 0
        while i < len(motif):
            if motif[i] == '{':
                jend = motif.index('}', i)
                rules.append(("not", motif[i+1:jend]))
                i = jend + 1
            elif motif[i] == '[':
                jend = motif.index(']', i)
                rules.append(("one of", motif[i+1:jend]))
                i = jend + 1
            else:
                rules.append(("exact", motif[i]))
                i += 1
        return rules

    @staticmethod
    def motif_window_match(seq_window: str, rules: List[tuple] = None) -> bool:
        """Checks if a window of a sequence matches the given motif rules"""

        for i in range(len(rules)):
            rule_type, letters = rules[i]
            if rule_type == "not":
                if seq_window[i] in letters:
                    return False
            elif rule_type == "one of":
                if seq_window[i] not in letters:
                    return False
            elif rule_type == "exact":
                if seq_window[i] != letters:
                    return False
        return True

    def motif_locations(self, motif_rules: List[tuple] = None) -> List[int]:
        """
        Returns a list of locations in a sequence where a motif's rules are satisfied
        """

        if motif_rules is None:
            motif_rules = []

        len_m = len(motif_rules)
        motif_locations = []
        for i in range(len(self.seq) - len_m + 1):
            window = self.seq[i:i + len_m]
            if BioSequence.motif_window_match(window, motif_rules):
                motif_locations.append(i + 1)
        return motif_locations

    def splice_introns(self, introns: List[str] = None) -> BioSequence:
        """Splice out introns from a sequence and return the exons"""

        if introns is None:
            introns = []

        exons_seq = self.seq
        for i in introns:
            exons_seq = exons_seq.replace(i, '')
        return BioSequence(exons_seq)

    def k_mers_from_sequence(self, n_length: int) -> List[str]:
        """Returns all k-mers of length n in a sequence"""

        if n_length is None:
            n_length = 1

        k_mers = []
        for i in range(len(self.seq) - n_length + 1):
            k_mer = self.seq[i:i+n_length]
            k_mers.append(k_mer)
        return k_mers

    @staticmethod
    def k_mers_from_collection(collection: str, n_length: int) -> List[str]:
        """Returns all possible k-mers of length n from a collection of symbols"""

        if n_length is None:
            n_length = 1

        collection_sorted = ''.join(sorted(collection))
        k_mers = [''.join(p) for p in product(collection_sorted, repeat=n_length)]
        return k_mers
