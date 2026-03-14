
"""Rosalind Bioinformatics Stronghold - Problem solver runner"""

from pathlib import Path

import biostronghold

# Path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / 'Datasets'

print('Counting DNA Nucleotides')
file_path = DATASETS_DIR / 'rosalind_dna.txt'
print(biostronghold.counting_DNA_nucleotides(file_path))

print('\nTranscribing DNA into RNA')
file_path = DATASETS_DIR / 'rosalind_rna.txt'
print(biostronghold.transcribing_DNA_into_RNA(file_path))

print("\nComplementing a strand of DNA")
file_path = DATASETS_DIR / 'rosalind_revc.txt'
print(biostronghold.complementing_a_strand_of_DNA(file_path))

print("\nRabbits and Recurrence Relations")
file_path = DATASETS_DIR / 'rosalind_fib.txt'
print(biostronghold.rabbits_and_recurrence_relations(file_path))

print("\nComputing GC Content")
fasta_path = DATASETS_DIR / "rosalind_gc.txt"
print(biostronghold.computing_GC_content(fasta_path))

print("\nCounting Point Mutations")
file_path = DATASETS_DIR / "rosalind_hamm.txt"
print(biostronghold.counting_point_mutations(file_path))

print("\nMendel's First Law")
file_path = DATASETS_DIR / 'rosalind_iprb.txt'
print(biostronghold.mendel_first_law(file_path))

print("\nTranslating RNA into Protein")
file_path = DATASETS_DIR / "rosalind_prot.txt"
print(biostronghold.translating_RNA_into_protein(file_path))

print("\nFinding a motif in DNA")
file_path = DATASETS_DIR / "rosalind_subs.txt"
print(biostronghold.finding_a_motif_in_DNA(file_path))

print("\nConsensus and profile")
fasta_path = DATASETS_DIR / "rosalind_cons.txt"
print(biostronghold.consensus_and_profile(fasta_path))

print("\nMortal Fibonacci Rabbits")
file_path = DATASETS_DIR / "rosalind_fibd.txt"
print(biostronghold.mortal_fibonacci_rabbits(file_path))

print("\nOverlap Graph")
fasta_path = DATASETS_DIR / "rosalind_grph.txt"
print(biostronghold.overlap_graphs(fasta_path))

print("\nCalculating Expected Offsprings")
file_path = DATASETS_DIR / "rosalind_iev.txt"
print(biostronghold.calculating_expected_offsprings(file_path))

print("\nFinding a Shared Motif")
filfastath = DATASETS_DIR / "rosalind_lcsm.txt"
print(biostronghold.finding_a_shared_motif(fasta_path))

print("\nIndependent Alleles")
file_path = DATASETS_DIR / "rosalind_lia.txt"
print(biostronghold.independent_alleles(file_path))

# print("\nFinding a Protein Motif")
# file_path = DATASETS_DIR/'rosalind_mprt.txt'
# print(biostronghold.finding_a_protein_motif(file_path))

print("\nInferring mRNA from Protein")
file_path = DATASETS_DIR/'rosalind_mrna.txt'
print(biostronghold.inferring_mRNA_from_protein(file_path))

print("\nOpen Reading Frames")
fasta_path = DATASETS_DIR / "rosalind_orf.txt"
print(biostronghold.open_reading_frames(fasta_path))

print("\nEnumerating Gene Orders")
file_path = DATASETS_DIR / "rosalind_perm.txt"
print(biostronghold.enumerating_gene_orders(file_path))

print("\nCalculating Protein Mass")
file_path = DATASETS_DIR / "rosalind_prtm.txt"
print(biostronghold.calculating_protein_mass(file_path))

print("\nLocating Restriction Sites")
fasta_path = DATASETS_DIR / 'rosalind_revp.txt'
print(biostronghold.locating_restriction_sites(fasta_path))

print("\nRNA Splicing")
fasta_path = DATASETS_DIR / 'rosalind_splc.txt'
print(biostronghold.rna_splicing(fasta_path))

print("\nEnumerating k-mers Lexicographically")
file_path = DATASETS_DIR / 'rosalind_lexf.txt'
print(biostronghold.enumerating_k_mers_lexicographically(file_path))

print("\nLongest Increasing Subsequence")
file_path = DATASETS_DIR / 'rosalind_lgis.txt'
print(biostronghold.longest_increasing_subsequence(file_path))
