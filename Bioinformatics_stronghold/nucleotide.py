"""
Base class for biological sequence manipulation (DNA, RNA).
"""

from typing import List
import sys
from pathlib import Path

# Add parent directory to Python path for imports
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from config import BASE_NUCLEOTIDES, DNA_CODON_TO_AMINO, RNA_CODON_TO_AMINO
from biosequence import BioSequence
from protein import Protein


class Nucleotide(BioSequence):
    """Class for manipulating nucleotide sequences (DNA, RNA)"""

    def __init__(self, seq: str = 'ATGC', seq_type: str = 'DNA'):
        self.seq: str = seq.upper()
        self.seq_type: str = seq_type.upper()
        self._validate()

    def _validate(self) -> None:
        """Ensures the sequence matches the allowed nucleotides for its type."""
        allowed_nucleotides = BASE_NUCLEOTIDES.get(self.seq_type)
        if not allowed_nucleotides:
            raise ValueError(f"Unsupported sequence type: {self.seq_type}")

        seq_set = set(self.seq)
        if not set(seq_set).issubset(allowed_nucleotides):
            invalid_chars = set(seq_set) - set(allowed_nucleotides)
            raise ValueError(
                f"Invalid characters {invalid_chars} for {self.seq_type} sequence.")

    def transcribe(self) -> Nucleotide:
        """Returns the transcribed RNA sequence object from DNA sequence
        or the starting DNA object from RNA sequence"""
        if self.seq_type == 'DNA':
            new_seq = self.seq.replace("T", "U")
            new_seq_type = 'RNA'
        elif self.seq_type == 'RNA':
            new_seq = self.seq.replace("U", "T")
            new_seq_type = 'DNA'
        else:
            raise ValueError("Invalid sequence type! Transcription not performed!")
        return Nucleotide(new_seq, new_seq_type)

    def complement(self) -> Nucleotide:
        """Returns the complement of a sequence"""
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        elif self.seq_type == "RNA":
            mapping = str.maketrans('AUCG', 'UAGC')
        else:
            raise ValueError("Invalid sequence type! Complement not performed!")
        new_seq = self.seq.translate(mapping)
        return Nucleotide(new_seq, self.seq_type)

    def reverse_complement(self) -> Nucleotide:
        """Returns the reverse complement of a sequence"""
        complement_nuc_seq = self.complement()
        reverse_seq = complement_nuc_seq.seq[::-1]
        return Nucleotide(reverse_seq, self.seq_type)

    def gc_content(self) -> float:
        """Returns the GC content % of a sequence"""
        return self.monomer_content('C') + self.monomer_content('G')

    def translate(self, i_start: int = 0) -> Protein:
        """Returns amino-acid sequence from translation of sequence"""
        if self.seq_type == 'DNA':
            codon_dict = DNA_CODON_TO_AMINO
        elif self.seq_type == 'RNA':
            codon_dict = RNA_CODON_TO_AMINO
        else:
            raise ValueError("Invalid sequence type! Translation not performed!")
        codon_len = 3
        protein = []
        for i in range(i_start, len(self.seq) - codon_len + 1, codon_len):
            codon = self.seq[i:i + codon_len]
            amino = codon_dict[codon]
            if amino == '_':
                break
            protein.append(amino)
        protein_seq = ''.join(protein)
        return Protein(protein_seq)

    def reading_frame(self, i_start: int = 0) -> str:
        """Returns a reading frame from a sequence"""

        if self.seq_type == 'DNA':
            codon_dict = DNA_CODON_TO_AMINO
        elif self.seq_type == 'RNA':
            codon_dict = RNA_CODON_TO_AMINO
        else:
            return ''

        codon_len = 3

        rf = []
        for i in range(i_start, len(self.seq) - codon_len + 1, codon_len):
            codon = self.seq[i:i + codon_len]
            amino = codon_dict[codon]
            rf.append(amino)
        return ''.join(rf)

    def all_reading_frames(self) -> List[str]:
        """Returns all 6 reading frames (3 forward + 3 reverse complement)"""

        reverse_complement = self.reverse_complement()

        reading_frames = []
        for i in range(3):
            reading_frames.append(self.reading_frame(i))
        for i in range(3):
            reading_frames.append(reverse_complement.reading_frame(i))
        return reading_frames

    def reverse_palindromes(self, l_min = 4, l_max = 12) -> List[tuple]:
        """Find locations of reverse palindromes in a DNA sequence"""

        palindromes = []
        for i in range(len(self.seq)):
            for l in range(l_min, l_max+1):
                if i+l > len(self.seq):
                    break
                sub_seq = Nucleotide(self.seq[i:i+l], self.seq_type)
                reverse_complement = sub_seq.reverse_complement()
                if sub_seq.seq == reverse_complement.seq:
                    palindromes.append((i+1, l))
        return palindromes
