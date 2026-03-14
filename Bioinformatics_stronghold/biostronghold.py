"""
BioStronghold solver - Solutions for Rosalind BioInformatics Stronghold problems.
"""

from typing import List
import sys
from pathlib import Path

# Add parent directory to Python path for imports
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from core import file_to_str_list, fasta_path_to_dict, file_to_int_list
from biosequence import BioSequence
from nucleotide import Nucleotide
from protein import Protein
from biosequencecollection import BioSequenceCollection
from nucleotidecollection import NucleotideCollection
from proteincollection import ProteinCollection


def counting_DNA_nucleotides(file_path: str) -> str:
    """Count nucleotides in a DNA sequence"""

    dna_seq = file_to_str_list(file_path)[0]
    bio_seq = BioSequence(dna_seq)

    count_dict = bio_seq.monomers_count()

    target_order = ['A', 'C', 'G', 'T']
    count_dict = {k: count_dict[k] for k in target_order}
    output = [f"{key}: {value}" for key, value in count_dict.items()]
    # output = [f" {value}" for value in count_dict.values()]
    return ' '.join(output)

def transcribing_DNA_into_RNA(file_path: str) -> str:
    """Transcribe a DNA sequence to RNA"""

    dna_seq = file_to_str_list(file_path)[0]
    nuc_seq = Nucleotide(dna_seq, 'DNA')

    return nuc_seq.transcribe().seq

def complementing_a_strand_of_DNA(file_path: str) -> str:
    """Returns the reverse complement of DNA sequence"""

    dna_seq = file_to_str_list(file_path)[0]
    nuc_seq = Nucleotide(dna_seq, 'DNA')

    return nuc_seq.reverse_complement().seq

def rabbits_and_recurrence_relations(file_path: str) -> int:
    """Returns the number of Fibonacci rabbits
    at month n given reproduction factor k"""

    file_list = file_to_int_list(file_path)
    n, k = file_list[0][0], file_list[0][1]

    f2, f1 = 1, 1
    for _ in range(3, n + 1):
        f1, f2 = f1 + k * f2, f1
    return f1

def computing_GC_content(fasta_path: str) -> str:
    """Returns the DNA sequence with the highest GC content"""

    dna_seq_dict = fasta_path_to_dict(fasta_path)
    nuc_seq_col = NucleotideCollection(dna_seq_dict, 'DNA')

    gc_content_dict = nuc_seq_col.gc_content_dict()

    max_key = max(gc_content_dict, key=gc_content_dict.get)
    return f"{max_key}\n{gc_content_dict[max_key]}"

def counting_point_mutations(file_path: str) -> int:
    """Returns the number of point mutations between two sequences"""

    dna_seq = file_to_str_list(file_path)[0]
    mutation_seq = file_to_str_list(file_path)[1]
    bio_seq = BioSequence(dna_seq)

    return bio_seq.point_mutations(mutation_seq)

def mendel_first_law(file_path: str) -> float:
    """Calculate the probability of dominant phenotype offspring"""

    file_list = file_to_int_list(file_path)
    k, m, n = file_list[0][0], file_list[0][1], file_list[0][2]

    t = k + m + n
    pr_yy = ( (n*(n-1)) + (n*m) + (m*(m-1)/4) )/( t*(t-1) )
    return 1 - pr_yy

def translating_RNA_into_protein(file_path: str) -> str:
    """Translate an RNA sequence into protein"""

    rna_seq = file_to_str_list(file_path)[0]
    nuc_seq = Nucleotide(rna_seq, 'RNA')

    protein_seq = nuc_seq.translate(i_start=0)
    return protein_seq.seq

def finding_a_motif_in_DNA(file_path: str) -> str:
    """Find all locations of a motif in a DNA sequence"""

    dna_seq = file_to_str_list(file_path)[0]
    bio_seq = BioSequence(dna_seq)

    dna_motif = file_to_str_list(file_path)[1]
    motif_rules = BioSequence.motif_rules(dna_motif)

    motif_locations = bio_seq.motif_locations(motif_rules)
    return " ".join(str(x) for x in motif_locations)

def consensus_and_profile(fasta_path: str) -> str:
    """Generate the consensus string and profile matrix of a set of DNA sequences"""

    dna_seq_dict = fasta_path_to_dict(fasta_path)
    bio_seq_col = BioSequenceCollection(dna_seq_dict)

    monomers = ["A", "C", "G", "T"]
    profile_matrix = bio_seq_col.profile_matrix(monomers)
    consensus = bio_seq_col.consensus_string(profile_matrix, monomers)

    output = [consensus]
    for nuc, counts in zip(monomers, profile_matrix):
        output.append(f"{nuc}: " + " ".join(str(x) for x in counts))
    return '\n'.join(output)

def mortal_fibonacci_rabbits(file_path: str) -> int:
    """Calculate mortal rabbit population at the nth month
    considering that rabbits live for m months"""

    file_list = file_to_int_list(file_path)
    n_months, m_ages = file_list[0][0], file_list[0][1]

    ages = [0]*m_ages
    ages[0] = 1
    for _ in range(1, n_months):
        newborns = sum(ages[1:])
        ages[1:] = ages[:-1]
        ages[0] = newborns
    return sum(ages)

def overlap_graphs(fasta_path: str) -> str:
    """Find overlap graph edges of order k=3 of a dictionary sequences"""

    dna_seq_dict = fasta_path_to_dict(fasta_path)
    bio_seq_col = BioSequenceCollection(dna_seq_dict)

    edges = bio_seq_col.adjacency_list(k=3)
    return "\n".join(edges)

def calculating_expected_offsprings(file_path: str) -> float:
    """Calculate the expected number of dominant offspring (AA or Aa)
    from pairs of organisms which reproduce 2 offsprings each"""

    # order: AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa
    dom_prop = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
    n_offs = 2

    file_list = file_to_int_list(file_path)
    couples = file_list[0]

    e_dom_offs = 0
    for i in range(len(couples)):
        e_dom_offs += couples[i] * dom_prop[i]
    return n_offs * e_dom_offs

def finding_a_shared_motif(fasta_path: str) -> str:
    """Find the longest shared motif in a set of sequences"""

    dna_seq_dict = fasta_path_to_dict(fasta_path)
    bio_seq_col = BioSequenceCollection(dna_seq_dict)

    return bio_seq_col.longest_shared_motif()

def independent_alleles(file_path: str) -> float:
    """Calculate the probability of at least n_min dominant alleles
    in the kth generation considering 2 offsprings per generation"""

    from math import comb

    file_list = file_to_int_list(file_path)
    gen, n_min, n_offs = file_list[0][0], file_list[0][1], 2
    prop = 0.25
    n_total = n_offs ** gen

    total_p = 0
    for s in range(n_min, n_total + 1):
        total_p += comb(n_total, s) * (prop ** s) * ((1 - prop) ** (n_total - s))
    return total_p

def finding_a_protein_motif(file_path: str) -> str:
    """Find protein motif locations in a sequence of proteins
    obtained from UniProt"""

    protein_motif = 'N{P}[ST]{P}'
    motif_rules = BioSequence.motif_rules(protein_motif)

    uniprot_ids = file_to_str_list(file_path)
    seq_dict = {}
    for uid in uniprot_ids:
        seq_dict[uid] = Protein.get_from_uniprot(uid).seq
    protein_seq_col = ProteinCollection(seq_dict)

    motif_locations_dict = protein_seq_col.motif_locations_dict(motif_rules)

    rosalind_output = []
    for key, value in motif_locations_dict.items():
        rosalind_output.append(key)
        rosalind_output.append(' '.join(str(x) for x in value))
    return '\n'.join(rosalind_output)

def inferring_mRNA_from_protein(file_path: str) -> int:
    """Count possible mRNA sequences for a protein"""

    seq = file_to_str_list(file_path)[0]
    protein_seq = Protein(seq)

    return protein_seq.mRNA_sequences_count()

def open_reading_frames(fasta_path: str) -> str:
    """Find all open reading frames in DNA sequence"""

    dna_seq = list(fasta_path_to_dict(fasta_path).values())[0]
    nuc_seq = Nucleotide(dna_seq, 'DNA')

    reading_frames = nuc_seq.all_reading_frames()

    proteins = ProteinCollection.all_proteins_from_all_reading_frames(reading_frames)
    return '\n'.join( [p.seq for p in proteins] )

def enumerating_gene_orders(file_path: str) -> List[tuple]:
    """Calculate number of gene orders for n genes"""
    import itertools

    n = file_to_int_list(file_path)[0][0]

    numbers = [x for x in range(1, n + 1)]

    permutations = itertools.permutations(numbers)
    permutations_list = list(permutations)

    print(len(permutations_list))
    # permutations_str = ''
    # for p in permutations_list:
    #     for x in p:
    #         permutations_str += str(x)
    #         permutations_str += ' '
    #     permutations_str += '\n'
    permutations_str = '\n'.join(' '.join(map(str, p)) for p in permutations_list)
    return permutations_str

def calculating_protein_mass(file_path: str) -> float:
    """Calculate the total mass of a protein sequence"""

    protein = file_to_str_list(file_path)[0]
    protein_seq = Protein(protein)

    return protein_seq.get_protein_mass()

def locating_restriction_sites(fasta_path: str) -> str:
    """Find locations of restriction sites in DNA sequence"""

    l_min, l_max = 4, 12

    seq_dict = fasta_path_to_dict(fasta_path)
    dna_seq = list(seq_dict.values())[0]
    nuc_seq = Nucleotide(dna_seq, 'DNA')

    palindromes = nuc_seq.reverse_palindromes(l_min, l_max)

    rosalind_output = [f"{p[0]} {p[1]}" for p in palindromes]
    return '\n'.join(rosalind_output)

def rna_splicing(fasta_path: str) -> str:
    """Splice introns from DNA sequence and translate to protein"""

    fasta_dict = fasta_path_to_dict(fasta_path)
    seqs = list(fasta_dict.values())

    dna_seq = seqs[0]
    bio_seq = BioSequence(dna_seq)

    introns = seqs[1:]

    bio_seq_spliced = bio_seq.splice_introns(introns)

    dna_seq_spliced = Nucleotide(bio_seq_spliced.seq, 'DNA')

    return dna_seq_spliced.translate(i_start=0).seq

def enumerating_k_mers_lexicographically(file_path: str) -> str:
    """Enumerate all k-mers of given length
    from a collection of symbols in lexicographic order."""

    collection_list = file_to_str_list(file_path)
    symbols = collection_list[0].replace(' ','')
    n_length = int(collection_list[1])

    k_mers = BioSequence.k_mers_from_collection(symbols, n_length)
    return '\n'.join(k_mers)

def longest_increasing_subsequence(file_path: str):

    file_list = file_to_int_list(file_path)
    perm = file_list[1]

    n = len(perm)
    # inc_len[i] stores the length of the longest increasing subseq ending at i
    inc_len = [1] * n
    # inc_prev[i] stores the index of the previous element so we can reconstruct it
    inc_prev = [-1] * n

    dec_len = [1] * n
    dec_prev = [-1] * n

    for i in range(n):
        for j in range(i):
            if perm[j] < perm[i] and inc_len[j] + 1 > inc_len[i]:
                inc_len[i] = inc_len[j] + 1
                inc_prev[i] = j
            if perm[j] > perm[i] and dec_len[j] + 1 > dec_len[i]:
                dec_len[i] = dec_len[j] + 1
                dec_prev[i] = j

    def reconstruct(lengths, predecessors):
        max_idx = lengths.index(max(lengths))
        path = []
        while max_idx != -1:
            path.append(str(perm[max_idx]))
            max_idx = predecessors[max_idx]
        return " ".join(path[::-1])

    return reconstruct(inc_len, inc_prev) + "\n" + reconstruct(dec_len, dec_prev)
