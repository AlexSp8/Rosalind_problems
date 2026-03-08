"""
BioStronghold solver - Solutions for Rosalind Bioinformatics Stronghold problems.
"""

from typing import List
import sys
from pathlib import Path

# Add parent directory to Python path for imports
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from core import fasta_path_to_dict
from biosequence import Biosequence


class BioStronghold(Biosequence):
    """Class that contains functions to solve all
      Bioinformatics BioStronghold problems in Rosalind"""

    def __init__(self, seq: str = "ATGC", seq_type: str = "DNA"):
        """Initialize BioStronghold solver with a sequence"""
        super().__init__(seq, seq_type)

    def counting_DNA_nucleotides(self) -> str:
        """Count DNA nucleotides and return formatted result"""
        count_dict = super().get_nucleotides_count()
        rosalind_output = [f"{key}: {value}"
                            for key, value in count_dict.items()]
        return ' '.join(rosalind_output)

    def transcribing_DNA_into_RNA(self) -> str:
        """Transcribe DNA to RNA"""
        return super().transcribe_sequence()

    def complementing_a_strand_of_DNA(self) -> str:
        """Get reverse complement of DNA sequence"""
        return super().reverse_complement_sequence()

    @staticmethod
    def rabbits_and_recurrence_relations(n: int, k: int) -> int:
        """Calculate Fibonacci with rabbits given reproduction factor k"""
        f2, f1 = 1, 1
        for _ in range(3, n + 1):
            f1, f2 = f1 + k * f2, f1
        return f1

    def computing_GC_content(self, file_path: str) -> str:
        """Find sequence with highest GC content"""
        fasta_dict = fasta_path_to_dict(file_path)
        gc_content_dict = {key: super().get_gc_content(seq)
                            for key, seq in fasta_dict.items()}
        max_key = max(gc_content_dict, key=gc_content_dict.get)
        return f"{max_key}\n{gc_content_dict[max_key]}"

    def counting_point_mutations(self, mutation_seq: str) -> int:
        """Count point mutations between sequences"""
        return super().get_point_mutations(mutation_seq)

    @staticmethod
    def mendel_first_law(k: int, m: int, n: int) -> float:
        """Calculate probability of dominant phenotype offspring"""
        t = k + m + n
        pr_yy = ( (n*(n-1)) + (n*m) + (m*(m-1)/4) )/( t*(t-1) )
        return 1 - pr_yy

    def translating_RNA_into_protein(self) -> str:
        """Translate RNA to protein"""
        return super().translate_sequence(istart=0)

    def finding_a_motif_in_DNA(self, motif: str) -> str:
        """Find all locations of motif in DNA sequence"""
        return super().motif_locations_in_sequence(motif)

    def consensus_and_profile(self, file_path: str) -> str:
        """Generate consensus string and nucleotide profile"""
        fasta_dict = fasta_path_to_dict(file_path)
        seq_list = list(fasta_dict.values())
        profile_matrix = super().profile_matrix(seq_list)
        consensus = super().consensus_string(profile_matrix)

        nucleotides = ["A", "C", "G", "T"]
        rosalind_output = [consensus]
        for nuc, counts in zip(nucleotides, profile_matrix):
            rosalind_output.append(
                f"{nuc}: " + " ".join(str(x) for x in counts))
        return '\n'.join(rosalind_output)

    @staticmethod
    def mortal_fibonacci_rabbits(n_months: int, m_ages: int) -> int:
        """Calculate mortal rabbit population growth"""
        ages = [0]*m_ages
        ages[0] = 1
        for _ in range(1, n_months):
            newborns = sum(ages[1:])
            ages[1:] = ages[:-1]
            ages[0] = newborns
        return sum(ages)

    def overlap_graphs(self, file_path: str, k: int = 3) -> str:
        """Find overlap graph edges between sequences"""
        seq_dict = fasta_path_to_dict(file_path)
        edges = super().adjacency_list(seq_dict, k=k)
        return "\n".join(edges)

    @staticmethod
    def calculating_expected_offsprings(couples: List[int], n_offs: int = 2) -> float:
        """Calculate expected dominant offspring from couples"""
        # order: AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa
        dom_prop = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
        e_dom_offsprings = 0
        for i in range(len(couples)):
            e_dom_offsprings += couples[i] * dom_prop[i]
        return n_offs * e_dom_offsprings

    def finding_a_shared_motif(self, file_path: str) -> str:
        """Find longest shared motif in multiple sequences"""
        seq_list = list(fasta_path_to_dict(file_path).values())
        return super().longest_shared_motif(seq_list)

    @staticmethod
    def independent_alleles(gen: int, n_min: int, prop: float) -> float:
        """Calculate probability of at least n_min dominant alleles in generation"""
        from math import comb
        n = 2 ** gen
        total_p = 0
        for s in range(n_min, n + 1):
            total_p += comb(n, s) * (prop ** s) * ((1 - prop) ** (n - s))
        return total_p

    def finding_a_protein_motif(self, uniprot_ids: List[str], protein_motif: str) -> str:
        """Find protein motif locations in UniProt sequences"""
        motif_rules_list = super().protein_motif_rules(protein_motif)
        locations_dict = super().get_protein_seq_from_uniprot(uid)
        locations_dict = {}
        for uid in uniprot_ids:
            protein_str = super().get_protein_seq_from_uniprot(uid)
            locations_list = super().protein_matching_motif_locations(protein_str, motif_rules_list)
            if locations_list:
                locations_dict[uid] = locations_list

        rosalind_output = []
        for key, value in locations_dict.items():
            rosalind_output.append(key)
            rosalind_output.append(' '.join(str(x) for x in value))
        return '\n'.join(rosalind_output)

    def inferring_mRNA_from_protein(self, protein_seq: str) -> int:
        """Count possible mRNA sequences for a protein"""
        return super().get_mRNA_seqs_from_protein(protein_seq)

    def open_reading_frames(self, fasta_path: str) -> str:
        """Find all open reading frames in DNA sequence"""
        seq = list(fasta_path_to_dict(fasta_path).values())[0]
        rframes = super().get_sequence_reading_frames(seq, seq_type='DNA')
        total_proteins = super().get_proteins_from_all_reading_frames(rframes)
        return '\n'.join(total_proteins)

    @staticmethod
    def enumerating_gene_orders(n: int) -> List[tuple]:
        """Calculate number of gene orders for n genes"""
        import itertools
        numbers = [x for x in range(1, n + 1)]
        permutations = itertools.permutations(numbers)
        return list(permutations)

    def calculating_protein_mass(self, protein_seq: str) -> float:
        """Calculate the total mass of a protein sequence"""
        return super().get_protein_mass(protein_seq)
