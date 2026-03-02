
from utilities import fasta_to_dict
from structures import Base_Nucleotides, RNA_Aminos

class BioinformaticsStronghold():

    """Class that contains functions to solve all Bioinformatics Stronghold problems in Rosalind"""

    def __init__(self, seq="ATGC", seq_type="DNA"):
        self.seq = seq.upper()
        self.seq_type = seq_type
        self.is_valid = self.__validateSequence()
        assert self.is_valid, f"{self.seq} data is not a correct {self.seq_type} sequence!"

    #private method
    def __validateSequence(self):
        """Check and uppercase all characters"""
        return set(Base_Nucleotides[self.seq_type]).issuperset(self.seq)

    def counting_DNA_nucleotides(self):
        from collections import Counter
        return dict(Counter(self.seq))
        # seq_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
        # for nuc in self.seq:
        #     seq_dict[nuc] += 1
        # return seq_dict

    def transcribing_DNA_into_RNA(self):
        return self.seq.replace("T", "U")

    def complementing_a_strand_of_DNA(self):
        mapping = str.maketrans('ATCG', 'TAGC')
        return self.seq.translate(mapping)[::-1]

    @staticmethod
    def rabbits_and_recurrence_relations(n, k):
        f2, f1 = 1, 1
        for _ in range(3, n+1):
            # fn = f1 + k*f2
            # f2 = f1
            # f1 = fn
            f1, f2 = f1+k*f2, f1
        return f1

    @staticmethod
    def seq_GC_content(seq):
        gc_content = seq.count('C') + seq.count('G')
        gc_content = 100*gc_content/(len(seq))
        return gc_content

    @staticmethod
    def gc_content_to_dict(fasta_dict):
        gc_content_dict = {}
        for key, seq in fasta_dict.items():
            gc_content_dict[key] = BioinformaticsStronghold.seq_GC_content(seq)
        return gc_content_dict

    @staticmethod
    def max_gc_content(gc_content_dict):
        max_gc_key = max(gc_content_dict, key=gc_content_dict.get)
        max_gc_value = gc_content_dict[max_gc_key]
        return max_gc_key, max_gc_value

    @staticmethod
    def computing_GC_content(fasta_path):
        fasta_dict = fasta_to_dict(fasta_path)
        # print(fasta_dict)
        gc_content_dict = BioinformaticsStronghold.gc_content_to_dict(fasta_dict)
        # print(gc_content_dict)
        return BioinformaticsStronghold.max_gc_content(gc_content_dict)

    def counting_point_mutations(self, mutation_seq):
        if (len(self.seq) != len(mutation_seq)):
            raise ValueError("Sequences must be the same length!")
        point_mutations = 0
        for i in range(len(self.seq)):
            if (self.seq[i] != mutation_seq[i]):
                point_mutations += 1
        return point_mutations
        # return sum(c1 != c2 for c1, c2 in zip(self.seq,mutation_seq))

    @staticmethod
    def mendel_first_law(k, m, n):
        t = k+m+n
        pr_yy = ( (n*(n-1))+(n*m)+(m*(m-1)/4) )/(t*(t-1))
        return 1-pr_yy

    def translating_RNA_into_protein(self):
        protein = []
        for i in range(0, len(self.seq)-2, 3):
            codon = self.seq[i:i+3]
            amino = RNA_Aminos[codon]
            if amino == '_':
                break
            protein.append(amino)
        return ''.join(protein)

    def finding_a_motif_in_DNA(self, motif):
        m = len(motif)
        ilocs = []
        for i in range(0, len(self.seq)-m+1):
            subseq = self.seq[i:i+m]
            if subseq == motif:
                ilocs.append(i+1)
        return " ".join(map(str, ilocs))
        # return " ".join(str(i) for i in ilocs)

    @staticmethod
    def profile_matrix(sequences):
        nucleotides = ["A", "C", "G", "T"]
        nuc_index = {n: i for i, n in enumerate(nucleotides)}

        rows = len(nucleotides)
        cols = max(len(seq) for seq in sequences)

        profile_matrix = [ [0]*cols for _ in range(rows) ]
        for seq in sequences:
            for j, nuc in enumerate(seq):
                profile_matrix[nuc_index[nuc]][j] += 1
        return profile_matrix

    @staticmethod
    def consensus_string(profile_matrix):
        nucleotides = ["A", "C", "G", "T"]
        rows = len(nucleotides)
        cols = len(profile_matrix[0])
        consensus = []
        for j in range(cols):
            row = max(range(rows), key=lambda i: profile_matrix[i][j])
            consensus.append(nucleotides[row])
        return "".join(consensus)

    @staticmethod
    def consensus_and_profile(fasta_path):
        fasta_dict = fasta_to_dict(fasta_path)

        sequences = list(fasta_dict.values())

        profile_matrix = BioinformaticsStronghold.profile_matrix(sequences)

        consensus = BioinformaticsStronghold.consensus_string(profile_matrix)

        nucleotides = ["A", "C", "G", "T"]
        rosalind_output = [consensus]
        for nuc, counts in zip(nucleotides, profile_matrix):
            rosalind_output.append(f"{nuc}: "+" ".join(str(x) for x in counts))

        return '\n'.join(rosalind_output)

    @staticmethod
    def mortal_fibonacci_rabbits(n, m):
        ages = [0]*m
        ages[0] = 1
        for _ in range(1,n):
            newborns = sum(ages[1:])
            ages[1:] = ages[:-1]
            ages[0] = newborns
        return sum(ages)

    @staticmethod
    def overlap_graphs(fasta_path, k):
        fasta_dict = fasta_to_dict(fasta_path)

        suffix = {key: seq[-k:] for key, seq in fasta_dict.items()}
        prefix = {key: seq[:k] for key, seq in fasta_dict.items()}

        edges = []
        for s in fasta_dict:
            for t in fasta_dict:
                if s != t and suffix[s] == prefix[t]:
                    edges.append(f"{s[1:]} {t[1:]}")

        return "\n".join(edges)

    @staticmethod
    def calculating_expected_offsprings(couples):
        #order: AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa
        dominant_prop = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
        e_offsprings = 0
        for i in range(len(couples)):
            e_offsprings += couples[i]*dominant_prop[i]
        return 2*e_offsprings

    @staticmethod
    def finding_a_shared_motif(fasta_path):
        sequences = list(fasta_to_dict(fasta_path).values())
        shortest_seq = min(sequences, key=len)
        max_l = len(shortest_seq)
        for l in range(max_l,0,-1):
            for i in range(max_l-l+1):
                motif = shortest_seq[i:i+l]
                if all(motif in seq for seq in sequences):
                    return motif
