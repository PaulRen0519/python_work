def most_frequent_nucleotide(gene_sequence):
    stop_codons = {'UAA', 'UAG', 'UGA'}
    for i in range(0, len(gene_sequence), 3):
        codon = gene_sequence[i:i+3]
        if codon in stop_codons:
            actual_sequence = gene_sequence[:i]    
            break
    codon_count = {}
    for i in range(0, len(actual_sequence), 3):
        codon = actual_sequence[i:i+3]
        if len(codon) == 3:
            codon_count[codon] = codon_count.get(codon, 0) + 1
    max_count = 0
    most_frequent_codon = []
    for codon, count in codon_count.items():
        if count > max_count:
            max_count = count
    for codon, count in codon_count.items():
        if count == max_count:
            most_frequent_codon.append(codon)
    return "most frequent codon is: " + ' '.join(most_frequent_codon)
print(most_frequent_nucleotide("AUGGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
print(most_frequent_nucleotide("AUGCCUCCUGUUGUUUGA"))