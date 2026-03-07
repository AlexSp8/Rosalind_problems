"""
Base class for biological sequence manipulation (DNA, RNA).
"""

from typing import Optional, List, Dict
import sys
from pathlib import Path

# Add parent directory to Python path for imports
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from config import BASE_NUCLEOTIDES, RNA_CODON_TO_AMINO, DNA_CODON_TO_AMINO, RNA_AMINO_TO_CODON


class Biosequence:
    """Class for manipulating biological sequences (DNA, RNA)"""

    def __init__(self, seq: str = 'ATGC', seq_type: str = 'DNA', seq_list: Optional[List[str]] = None, seq_dict: Optional[Dict[str, str]] = None):
        if seq_list is None:
            seq_list = ['ATGC']
        if seq_dict is None:
            seq_dict = {'seq1': 'ATGC'}

        self.seq: str = seq.upper()
        self.seq_type: str = seq_type
        self.seq_list: List[str] = seq_list
        self.seq_dict: Dict[str, str] = seq_dict
        self.is_valid: bool = self.__validate_sequence()
        assert self.is_valid, f"{self.seq} data is not a correct {self.seq_type} sequence!"

    def __validate_sequence(self) -> bool:
        """Check nucleotides and uppercase all characters"""
        return set(BASE_NUCLEOTIDES[self.seq_type]).issuperset(self.seq)

    def return_target_seq(self, seq: Optional[str] = None) -> str:
        """Get target sequence or default to instance sequence."""
        return seq if seq is not None else self.seq

    def return_target_seq_type(self, seq_type: Optional[str] = None) -> str:
        """Get target sequence type or default to instance type."""
        return seq_type if seq_type is not None else self.seq_type

    def get_nucleotides_count(self, seq: Optional[str] = None) -> Dict[str, int]:
        """Returns the number of nucleotides in sequence"""
        target_seq = self.return_target_seq(seq)
        from collections import Counter
        return dict(Counter(target_seq))

    def transcribe_sequence(self, seq: Optional[str] = None, seq_type: Optional[str] = None) -> str:
        """Returns transcribed RNA from DNA sequence or starting DNA from RNA sequence"""
        target_seq = self.return_target_seq(seq)
        target_seq_type = self.return_target_seq_type(seq_type)
        if target_seq_type == "DNA":
            return target_seq.replace("T", "U")
        elif target_seq_type == 'RNA':
            return target_seq.replace("U", "T")
        else:
            raise ValueError("Invalid sequence type! Transcription not performed!")

    def complement_sequence(self, seq: Optional[str] = None, seq_type: Optional[str] = None) -> str:
        """Returns complement of sequence"""
        target_seq = self.return_target_seq(seq)
        target_seq_type = self.return_target_seq_type(seq_type)
        if target_seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        elif target_seq_type == "RNA":
            mapping = str.maketrans('AUCG', 'UAGC')
        else:
            mapping = ''
        return target_seq.translate(mapping)

    def reverse_complement_sequence(self, seq: Optional[str] = None, seq_type: Optional[str] = None) -> str:
        """Returns reverse complement of sequence"""
        return ''.join(self.complement_sequence(seq, seq_type))[::-1]

    def get_gc_content(self, seq: Optional[str] = None) -> float:
        """Returns GC content % of sequence"""
        target_seq = self.return_target_seq(seq)
        gc_content = target_seq.count('C') + target_seq.count('G')
        return gc_content * 100 / len(target_seq)

    def get_point_mutations(self, mutation_seq: str, seq: Optional[str] = None) -> int:
        """Count point mutations between sequences"""
        target_seq = self.return_target_seq(seq)
        if len(target_seq) != len(mutation_seq):
            raise ValueError("Sequences must be the same length!")
        return sum(c1 != c2 for c1, c2 in zip(target_seq, mutation_seq))

    def translate_sequence(self, istart: int = 0, seq: Optional[str] = None, seq_type: Optional[str] = None) -> str:
        """Returns amino-acid sequence from translation of sequence"""
        target_seq = self.return_target_seq(seq)
        target_seq_type = self.return_target_seq_type(seq_type)
        codon_len = 3
        protein = []
        for i in range(istart, len(target_seq) - codon_len + 1, codon_len):
            codon = target_seq[i:i + codon_len]
            if target_seq_type == 'DNA':
                amino = DNA_CODON_TO_AMINO[codon]
            elif target_seq_type == 'RNA':
                amino = RNA_CODON_TO_AMINO[codon]
            else:
                raise ValueError("Invalid sequence type! Translation not performed!")
            if amino == '_':
                break
            protein.append(amino)
        return ''.join(protein)

    def motif_locations_in_sequence(self, motif: str, seq: Optional[str] = None) -> str:
        """Returns all motif locations in a sequence (1-indexed)"""
        target_seq = self.return_target_seq(seq)
        m = len(motif)
        ilocs = [i + 1 for i in range(len(target_seq) - m + 1) if target_seq[i:i + m] == motif]
        return " ".join(map(str, ilocs))

    def profile_matrix(self, seq_list: Optional[List[str]] = None) -> List[List[int]]:
        """Returns the profile matrix of a list of sequences"""
        target_seq_list = seq_list if seq_list is not None else self.seq_list
        nucleotides = ["A", "C", "G", "T"]
        nuc_index = {n: i for i, n in enumerate(nucleotides)}
        rows = len(nucleotides)
        cols = max(len(seq) for seq in target_seq_list)
        profile_matrix = [[0] * cols for _ in range(rows)]
        for seq in target_seq_list:
            for j, nuc in enumerate(seq):
                profile_matrix[nuc_index[nuc]][j] += 1
        return profile_matrix

    @staticmethod
    def consensus_string(profile_matrix: List[List[int]]) -> str:
        """Returns the consensus string of a profile matrix"""
        nucleotides = ["A", "C", "G", "T"]
        rows = len(nucleotides)
        cols = len(profile_matrix[0])
        consensus = []
        for j in range(cols):
            row = max(range(rows), key=lambda i: profile_matrix[i][j])
            consensus.append(nucleotides[row])
        return "".join(consensus)

    def adjacency_list(self, seq_dict: Optional[Dict[str, str]] = None, k: int = 3) -> List[str]:
        """Returns a list of adjacent edges of a dictionary of sequences"""
        target_seq_dict = seq_dict if seq_dict is not None else self.seq_dict
        suffix = {key: seq[-k:] for key, seq in target_seq_dict.items()}
        prefix = {key: seq[:k] for key, seq in target_seq_dict.items()}
        edges = []
        for s in target_seq_dict:
            for t in target_seq_dict:
                if s != t and suffix[s] == prefix[t]:
                    edges.append(f"{s[1:]} {t[1:]}")
        return edges

    def longest_shared_motif(self, seq_list: Optional[List[str]] = None) -> Optional[str]:
        """Returns the longest common substring of a list of sequences"""
        target_seq_list = seq_list if seq_list is not None else self.seq_list
        shortest_seq = min(target_seq_list, key=len)
        max_l = len(shortest_seq)
        for l in range(max_l, 0, -1):
            for i in range(max_l - l + 1):
                motif = shortest_seq[i:i + l]
                if all(motif in seq for seq in target_seq_list):
                    return motif
        return None

    @staticmethod
    def get_protein_seq_from_uniprot(uid: str) -> str:
        """Returns the protein sequence for a given UniProt ID"""
        import urllib.request
        accession_id = uid.split('_')[0]
        url = f'http://www.uniprot.org/uniprotkb/{accession_id}.fasta'

        try:
            with urllib.request.urlopen(url) as response:
                fasta_data = response.read().decode('utf-8')
                lines_list = fasta_data.splitlines()
                return ''.join(lines_list[1:])
        except urllib.error.URLError as e:
            return f"Failed to fetch sequence: {e}"

    @staticmethod
    def protein_motif_rules(protein_motif: str) -> List[tuple]:
        """Returns a list of rules for a protein motif"""
        rules = []
        i = 0
        while i < len(protein_motif):
            if protein_motif[i] == '{':
                jend = protein_motif.index('}', i)
                rules.append(("not", protein_motif[i + 1:jend]))
                i = jend + 1
            elif protein_motif[i] == '[':
                jend = protein_motif.index(']', i)
                rules.append(("one of", protein_motif[i + 1:jend]))
                i = jend + 1
            else:
                rules.append(("exact", protein_motif[i]))
                i += 1
        return rules

    @staticmethod
    def protein_motif_match_window(protein_subseq: str, rules: List[tuple]) -> bool:
        """Checks if a protein subsequence matches the given protein motif rules"""
        for i in range(len(rules)):
            rule_type, letters = rules[i]
            if rule_type == "not":
                if protein_subseq[i] in letters:
                    return False
            elif rule_type == "one of":
                if protein_subseq[i] not in letters:
                    return False
            elif rule_type == "exact":
                if protein_subseq[i] != letters:
                    return False
        return True

    def protein_matching_motif_locations(self, protein_seq: Optional[str] = None, rules_list: Optional[List[tuple]] = None) -> List[int]:
        """Returns a list of locations in a protein sequence where the motif is found"""
        if rules_list is None:
            rules_list = []
        target_seq = self.return_target_seq(protein_seq)
        len_motif = len(rules_list)
        locations_list = []
        for i in range(len(target_seq) - len_motif + 1):
            protein_subseq = target_seq[i:i + len_motif]
            if Biosequence.protein_motif_match_window(protein_subseq, rules_list):
                locations_list.append(i + 1)
        return locations_list

    def get_mRNA_seqs_from_protein(self, protein_seq: Optional[str] = None) -> int:
        """Returns the total number of mRNA sequences that can possibly
        translate a protein sequence (modulo 10**6)"""
        full_protein_seq = self.return_target_seq(protein_seq) + '_'
        mRNA_seq_count = 1
        for aa in full_protein_seq:
            mRNA_seq_count = (mRNA_seq_count * len(RNA_AMINO_TO_CODON[aa])) % (10 ** 6)
        return mRNA_seq_count

    @staticmethod
    def get_reading_frame(seq: str, seq_type: str, istart: int = 0) -> str:
        """Returns a reading frame from a sequence"""
        codon_len = 3
        rframe = []
        for i in range(istart, len(seq) - codon_len + 1, codon_len):
            codon = seq[i:i + codon_len]
            if seq_type == 'DNA':
                amino = DNA_CODON_TO_AMINO[codon]
            elif seq_type == 'RNA':
                amino = RNA_CODON_TO_AMINO[codon]
            else:
                break
            rframe.append(amino)
        return ''.join(rframe)

    def get_sequence_reading_frames(self, seq: Optional[str] = None, seq_type: Optional[str] = None) -> List[str]:
        """Get all 6 reading frames (3 forward + 3 reverse complement)"""
        target_seq = self.return_target_seq(seq)
        target_seq_type = self.return_target_seq_type(seq_type)
        reverse_complement = self.reverse_complement_sequence(target_seq, target_seq_type)
        rframes = []
        for i in range(3):
            rframes.append(Biosequence.get_reading_frame(target_seq, target_seq_type, i))
        for i in range(3):
            rframes.append(Biosequence.get_reading_frame(reverse_complement, target_seq_type, i))
        return rframes

    @staticmethod
    def get_proteins_from_reading_frame(rframe: str) -> List[str]:
        """Extract proteins from a reading frame"""
        proteins = []
        for i in range(len(rframe)):
            if rframe[i] == 'M':
                protein = ''
                for amino in rframe[i:]:
                    if amino == '_':
                        proteins.append(protein)
                        break
                    protein += amino
        return proteins
