
from utilities import fasta_path_to_dict
import biostructures
from biosequence import Biosequence

class BioStronghold(Biosequence):

    """Class that contains functions to solve all Bioinformatics BioStronghold problems in Rosalind"""

    def __init__(self, seq="ATGC", seq_type="DNA"):
        super().__init__(seq, seq_type)

    def counting_DNA_nucleotides(self):
        count_dict = super().get_nucleotides_count()
        rosalind_output = [f"{key}: {value}" for key, value in count_dict.items()]
        return ' '.join(rosalind_output)

    def transcribing_DNA_into_RNA(self):
        return super().transcribe_sequence()

    def complementing_a_strand_of_DNA(self):
        return super().complement_sequence()[::-1]

    @staticmethod
    def rabbits_and_recurrence_relations(n, k):
        f2, f1 = 1, 1
        for _ in range(3, n+1):
            # fn = f1 + k*f2
            # f2 = f1
            # f1 = fn
            f1, f2 = f1+k*f2, f1
        return f1

    def computing_GC_content(self, file_path):
        fasta_dict = fasta_path_to_dict(file_path)
        gc_content_dict = {key: super().get_gc_content(seq) for key, seq in fasta_dict.items()}
        max_key = max(gc_content_dict, key=gc_content_dict.get)
        return f"{max_key}\n{gc_content_dict[max_key]}"

    def counting_point_mutations(self, mutation_seq):
        return super().get_point_mutations(mutation_seq)

    @staticmethod
    def mendel_first_law(k, m, n):
        t = k+m+n
        pr_yy = ( (n*(n-1))+(n*m)+(m*(m-1)/4) )/(t*(t-1))
        return 1-pr_yy

    def translating_RNA_into_protein(self):
        return super().translate_sequence(istart=0)

    def finding_a_motif_in_DNA(self, motif):
        return super().motif_locations_in_sequence(motif)

    def consensus_and_profile(self, file_path):
        fasta_dict = fasta_path_to_dict(file_path)
        sequences = list(fasta_dict.values())
        profile_matrix = super().profile_matrix(sequences)
        consensus = super().consensus_string(profile_matrix)

        nucleotides = ["A", "C", "G", "T"]
        rosalind_output = [consensus]
        for nuc, counts in zip(nucleotides, profile_matrix):
            rosalind_output.append(f"{nuc}: "+" ".join(str(x) for x in counts))
        return '\n'.join(rosalind_output)

    @staticmethod
    def mortal_fibonacci_rabbits(n_months, m_ages):
        ages = [0]*m_ages
        ages[0] = 1
        for _ in range(1,n_months):
            newborns = sum(ages[1:])
            ages[1:] = ages[:-1]
            ages[0] = newborns
        return sum(ages)

    def overlap_graphs(self, file_path, k=3):
        seqs_dict = fasta_path_to_dict(file_path)
        edges = super().adjacency_list(seqs_dict, k=k)
        return "\n".join(edges)

    @staticmethod
    def calculating_expected_offsprings(couples, n_offs=2):
        #order: AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa
        dom_prop = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
        e_dom_offsprings = 0
        for i in range(len(couples)):
            e_dom_offsprings += couples[i]*dom_prop[i]
        return n_offs*e_dom_offsprings

    def finding_a_shared_motif(self, file_path):
        sequences = list(fasta_path_to_dict(file_path).values())
        return super().longest_shared_motif(sequences)

    @staticmethod
    def independent_alleles(gen,n_min,prop):
        from math import comb
        n = 2**gen
        total_p = 0
        for s in range(n_min,n+1):
            total_p += comb(n, s)*(prop**s)*((1-prop)**(n-s))
        return total_p

    def finding_a_protein_motif(self, uniprot_ids,protein_motif):
        motif_rules_list = super().protein_motif_rules(protein_motif)
        locations_dict = {}
        for uid in uniprot_ids:
            protein_str = super().get_protein_seq_from_uniprot(uid)
            locations_list = super().protein_matching_motif_locations(protein_str,motif_rules_list)
            if locations_list:
                locations_dict[uid] = locations_list

        rosalind_output = []
        for key, value in locations_dict.items():
            rosalind_output.append(key)
            rosalind_output.append( ' '.join(str(x) for x in value) )
        return '\n'.join(rosalind_output)

    def inferring_mRNA_from_protein(self, protein_seq):
        return super().get_mRNA_seqs_from_protein(protein_seq)

    @staticmethod
    def proteins_from_reading_frame(reading_frame):
        proteins = []
        start_indices = [i for i, aa in enumerate(reading_frame) if aa == 'M']
        for istart in start_indices:
            remaining_rf = reading_frame[istart:]
            if '_' in remaining_rf:
                istop = remaining_rf.find('_')
                protein = remaining_rf[:istop]
                proteins.append(protein)
        return proteins

    def open_reading_frames(self, fasta_path):
        seq = list(fasta_path_to_dict(fasta_path).values())[0]
        rframes = super().get_sequence_reading_frames(seq, seq_type='DNA')
        # print(rframes)
        total_proteins = []
        for rf in rframes:
            rf_proteins = super().get_proteins_from_reading_frame(rf)
            if rf_proteins:
                for p in rf_proteins:
                    if p not in total_proteins:
                        total_proteins.append(p)
        return '\n'.join(total_proteins)
