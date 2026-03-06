
from biostronghold import BioStronghold
from utilities import file_to_list

seq = "GAGCCTACTAACGGGAT"
bio_strong = BioStronghold(seq, 'DNA')
file_path0 = "Bioinformatics_stronghold/Datasets/"

print("Counting DNA Nucleotides")
bio_strong.seq = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
bio_strong.seq_type = 'DNA'
print(bio_strong.counting_DNA_nucleotides())

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
# file_path = file_path0 + "rosalind_gc.txt"
# print(bio_strong.computing_GC_content(file_path))

# print("\nCounting Point Mutations")
# bio_strong.seq = 'GAGCCTACTAACGGGAT'
# bio_strong.seq_type = 'DNA'
# mutation_seq = "CATCGTAATGACGGCCT"
# print(bio_strong.counting_point_mutations(mutation_seq))

# print("\nMendel's First Law")
# k_hom_dom, m_het, n_hom_rec = 2, 2, 2
# print(bio_strong.mendel_first_law(k_hom_dom, m_het, n_hom_rec))

print("\nTranslating RNA into Protein")
bio_strong.seq = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
bio_strong.seq_type = 'RNA'
print(bio_strong.translating_RNA_into_protein())

# print("\nFinding a motif in DNA")
# bio_strong.seq = 'GATATATGCATATACTT'
# bio_strong.seq_type = 'DNA'
# motif = "ATAT"
# print(bio_strong.finding_a_motif_in_DNA(motif))

# print("\nConsensus and profile")
# file_path = file_path0 + "rosalind_cons.txt"
# print(bio_strong.consensus_and_profile(file_path))

# print("\nMortal Fibonacci Rabbits")
# print(bio_strong.mortal_fibonacci_rabbits(n_months=6, m_ages=3))

# print("\nOverlap Graph")
# file_path = file_path0+'rosalind_grph.txt'
# print(bio_strong.overlap_graphs(file_path, k=3))

# print("\nCalculating Expected Offsprings")
# couples = [1, 0, 0, 1, 0, 1]
# print(bio_strong.calculating_expected_offsprings(couples, n_offs=2))

# print("\nFinding a Shared Motif")
# file_path = file_path0+'rosalind_lcsm.txt'
# print(bio_strong.finding_a_shared_motif(file_path))

# print("\nIndependent Alleles")
# print(bio_strong.independent_alleles(gen=2, n_min=1, prop=0.25))

# print("\nFinding a Protein Motif")
# uniprot_ids = file_to_list(file_path0+'rosalind_mprt.txt')
# protein_motif = 'N{P}[ST]{P}'
# print(bio_strong.finding_a_protein_motif(uniprot_ids,protein_motif))

# print("\nInferring mRNA from Protein")
# print(bio_strong.inferring_mRNA_from_protein('MA'))

print("\nOpen Reading Frames")
fasta_path = file_path0+'rosalind_orf.txt'
print(bio_strong.open_reading_frames(fasta_path))
