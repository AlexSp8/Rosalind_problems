"""
Base class for biological sequence manipulation (DNA, RNA).
"""

import sys
from pathlib import Path

# Add parent directory to Python path for imports
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from config import RNA_AMINO_TO_CODON, MONOISOTROPIC_AMINO_MASS_TABLE
from biosequence import BioSequence


class Protein(BioSequence):
    """Class for manipulating biological sequences"""

    def __init__(self, seq: str = 'M'):
        self.seq: str = seq.upper()
        self._validate()

    def _validate(self) -> None:
        """Ensures the sequence matches the allowed amino acids."""

        allowed_aminos = RNA_AMINO_TO_CODON.keys()

        seq_set = set(self.seq)
        if not set(seq_set).issubset(allowed_aminos):
            invalid_chars = set(seq_set) - set(allowed_aminos)
            raise ValueError(f"Invalid characters {invalid_chars} for protein sequence.")

    @staticmethod
    def get_from_uniprot(uid: str) -> Protein:
        """Returns the protein sequence for a given UniProt ID"""

        import urllib.request

        accession_id = uid.split('_')[0]
        url = f'http://www.uniprot.org/uniprotkb/{accession_id}.fasta'

        try:
            with urllib.request.urlopen(url) as response:
                fasta_data = response.read().decode('utf-8')
                lines_list = fasta_data.splitlines()
                seq = ''.join(lines_list[1:])
                return Protein(seq)
        except urllib.error.URLError as e:
            return f"Failed to fetch sequence: {e}"

    def mRNA_sequences_count(self) -> int:
        """Returns the total number of mRNA sequences that can possibly
        translate a protein sequence (modulo 10**6)"""

        full_protein_seq = self.seq + '_'

        mRNA_seq_count = 1
        for aa in full_protein_seq:
            mRNA_seq_count = (mRNA_seq_count * len(RNA_AMINO_TO_CODON[aa])) % (10 ** 6)
        return mRNA_seq_count, int(mRNA_seq_count/(10**6))

    @staticmethod
    def first_protein_from_reading_frame(rf: str) -> Protein | None:
        """Returns the first possible protein from a reading frame"""

        if rf[0] != 'M':
            return None

        end_pos = rf.find('_')
        if end_pos != -1:
            protein = rf[:end_pos]
            return Protein(protein)

    def get_protein_mass(self) -> float:
        """Calculate total mass of a protein sequence"""

        total_mass = sum(MONOISOTROPIC_AMINO_MASS_TABLE[aa] for aa in self.seq)
        return total_mass
