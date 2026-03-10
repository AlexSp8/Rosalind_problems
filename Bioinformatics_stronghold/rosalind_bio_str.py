
"""Rosalind Bioinformatics Stronghold - Problem solver runner"""

import sys
from pathlib import Path

# Add parent directory to Python path for imports
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from core import file_to_list, file_to_int_list
from biostronghold import BioStronghold

# Get the path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / "Datasets"

# seq = "GAGCCTACTAACGGGAT"
seq = "GAGCGT"
bio_strong = BioStronghold(seq, 'DNA')

# print("Counting DNA Nucleotides")
# bio_strong.seq = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
# bio_strong.seq_type = 'DNA'
# print(bio_strong.counting_DNA_nucleotides())

# print("\nTranscribing DNA into RNA")
# bio_strong.seq = 'GATGGAACTTGACTACGTAAATT'
# bio_strong.seq_type = 'DNA'
# print(bio_strong.transcribing_DNA_into_RNA())

# print("\nComplementing a strand of DNA")
# bio_strong.seq = 'AAAACCCGGT'
# bio_strong.seq_type = 'DNA'
# print(bio_strong.complementing_a_strand_of_DNA())

# print("\nRabbits and Recurrence Relations")
# print(bio_strong.rabbits_and_recurrence_relations(n=5, k=3))

# print("\nComputing GC Content")
# file_path = DATASETS_DIR / "rosalind_gc.txt"
# print(bio_strong.computing_GC_content(file_path))

# print("\nCounting Point Mutations")
# bio_strong.seq = 'GAGCCTACTAACGGGAT'
# bio_strong.seq_type = 'DNA'
# mutation_seq = "CATCGTAATGACGGCCT"
# print(bio_strong.counting_point_mutations(mutation_seq))

# print("\nMendel's First Law")
# k_hom_dom, m_het, n_hom_rec = 2, 2, 2
# print(bio_strong.mendel_first_law(k_hom_dom, m_het, n_hom_rec))

# print("\nTranslating RNA into Protein")
# bio_strong.seq = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
# bio_strong.seq_type = 'RNA'
# print(bio_strong.translating_RNA_into_protein())

# print("\nFinding a motif in DNA")
# bio_strong.seq = 'GATATATGCATATACTT'
# bio_strong.seq_type = 'DNA'
# motif = "ATAT"
# print(bio_strong.finding_a_motif_in_DNA(motif))

# print("\nConsensus and profile")
# file_path = DATASETS_DIR / "rosalind_cons.txt"
# print(bio_strong.consensus_and_profile(file_path))

# print("\nMortal Fibonacci Rabbits")
# print(bio_strong.mortal_fibonacci_rabbits(n_months=6, m_ages=3))

# print("\nOverlap Graph")
# file_path = DATASETS_DIR / "rosalind_grph.txt"
# print(bio_strong.overlap_graphs(file_path, k=3))

# print("\nCalculating Expected Offsprings")
# couples = [1, 0, 0, 1, 0, 1]
# print(bio_strong.calculating_expected_offsprings(couples, n_offs=2))

# print("\nFinding a Shared Motif")
# file_path = DATASETS_DIR / "rosalind_lcsm.txt"
# print(bio_strong.finding_a_shared_motif(file_path))

# print("\nIndependent Alleles")
# print(bio_strong.independent_alleles(gen=2, n_min=1, prop=0.25))

# print("\nFinding a Protein Motif")
# uniprot_ids = file_to_list(DATASETS_DIR/'rosalind_mprt.txt')
# protein_motif = 'N{P}[ST]{P}'
# print(bio_strong.finding_a_protein_motif(uniprot_ids,protein_motif))

# print("\nInferring mRNA from Protein")
# print(bio_strong.inferring_mRNA_from_protein('MA'))

# print("\nOpen Reading Frames")
# fasta_path = DATASETS_DIR / "rosalind_orf.txt"
# print(bio_strong.open_reading_frames(fasta_path))

# print("\nEnumerating Gene Orders")
# permutations = bio_strong.enumerating_gene_orders(n=5)
# print(len(permutations))
# for p in permutations:
#     print(" ".join(str(x) for x in p))

# print("\nCalculating Protein Mass")
# protein_seq = "SKADYEK"
# print(bio_strong.calculating_protein_mass(protein_seq))

# print("\nLocating Restriction Sites")
# fasta_path = DATASETS_DIR / 'rosalind_revp.txt'
# print(bio_strong.locating_restriction_sites(fasta_path, seq_type='DNA', l_min=4, l_max=12))

# print("\nRNA Splicing")
# fasta_path = DATASETS_DIR / 'rosalind_splc.txt'
# print(bio_strong.rna_splicing(fasta_path, seq_type='DNA'))

# print("\nEnumerating k-mers Lexicographically")
# file_path = DATASETS_DIR / 'rosalind_lexf.txt'
# collection_list = file_to_list(file_path)
# symbols = collection_list[0].replace(' ','')
# n_length = int(collection_list[1])
# print(bio_strong.enumerating_k_mers_lexicographically(symbols, n_length))

print("\nLongest Increasing Subsequence")
file_path = DATASETS_DIR / 'rosalind_lgis.txt'
file_list = file_to_int_list(file_path)
permutation = file_list[1]
print(bio_strong.longest_increasing_subsequence(permutation))
