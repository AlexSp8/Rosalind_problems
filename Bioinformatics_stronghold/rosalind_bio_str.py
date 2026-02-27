
from bioinformaticsstronghold import BioinformaticsStronghold

seq = "AAGAAGATTAGTGATTAGTTGGATTAGTGATTAGTTGGATTAGTGTACGATTAGTATGAATGCCCCAGATTAGTAGATTAGTGATTAGTTGATTAGTTGGATTAGTGTGAGCTAGATTAGTCGATTAGTCCTGATTAGTTTTGATTGATTAGTTGAGTGGCGATTAGTAGATTAGTCGATTAGTGATTAGTTCGATTAGTGATTAGTTCGGATTAGTTCCGTCTGATTAGTAGACCCCACCCCGATTAGTATTGATTAGTGGGATTAGTGGGATTAGTACTGACTCACGATTAGTAGCTCGATTAGTGCCCTTGGATTAGTCCTTGATTAGTCAAGCAGAAGCGTTGATTAGTTAAGATTAGTCGGATTAGTTGCCGGATTAGTTAATTATCGATTAGTTCGCGATTAGTTCACCATGATTAGTGATTAGTGATTAGTCGATTAGTATGATTAGTAATAGAGATTAGTTAGATTAGTTAAGTGATTAGTCTGCGATTAGTCTGATTAGTTTGGATTAGTACGAGGATTAGTGGGATTAGTGATTAGTATCGATTAGTTGATTAGTATGCGGAGATTAGTATGATTAGTCATTACCTATTGATTAGTCGAATGTGGTGATTAGTGATTAGTGATTAGTCGATTAGTGATTAGTGATTAGTATTGATTAGTCCCTGATTAGTAGATTAGTGATTAGTGGATTAGTGACAGAGGATTAGTGGTGATTAGTATAGGATTAGTCTGATTAGTGCGGATTAGTTCGATTAGTGATTAGTCGGTGATTAGTCGGATTAGTGGATTAGTGATTAGTAAGATTAGTGCTACCAGATTAGTTGGGATTAGTGTTGGATTAGTTGGAGATTAGTCGCCAGCATCTAGATTAGTCGCCTGGGGTGGATTAGTCTTCGATTAGTCGATTAGTTACAGGATTAGTGATTAGTTGATTAGTACTTGATTAGTTGGATTAGT"
bio_str = BioinformaticsStronghold(seq, "DNA")

# print("Counting DNA Nucleotides")
# count_dict = bio_str.counting_DNA_nucleotides()
# for key, value in count_dict.items():
#     print(f"{key}: {value}", end=' ')

# print("\n\nTranscribing DNA into RNA")
# print(bio_str.transcribing_DNA_into_RNA())

# print("\nComplementing a strand of DNA")
# print(bio_str.complementing_a_strand_of_DNA())
# print(bio_str.seq)

# print("\nRabbits and Recurrence Relations")
# print(bio_str.rabbits_and_recurrence_relations(n=36, k=5))

# print("\nComputing GC Content")
# FASTA_path = "Bioinformatics_stronghold/Datasets/rosalind_gc.txt"
# print(bio_str.computing_GC_content(FASTA_path))

# print("\nCounting Point Mutations")
# mutation_seq = "CATCGTAATGACGGCCT"
# print(bio_str.counting_point_mutations(mutation_seq))

# print("\nMendel's First Law")
# k_hom_dom, m_het, n_hom_rec = 27, 28, 16
# print(bio_str.mendel_first_law(k_hom_dom, m_het, n_hom_rec))

# print("\nTranslating RNA into Protein")
# print(bio_str.translating_RNA_into_protein())

print("\nFinding a motif in DNA")
motif = "GATTAGTGA"
print(bio_str.finding_a_motif_in_DNA(motif))
