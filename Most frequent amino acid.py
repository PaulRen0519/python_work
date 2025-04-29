codon_table = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

def most_frequent_amino_acid(gene_sequence):
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
    results=[]
    for i in range(0, len(most_frequent_codon)):
        results.append(codon_table[most_frequent_codon[i]])
    return "most frequent amino acid(s): " + " ".join(results)
    
print(most_frequent_amino_acid("AUGGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
print(most_frequent_amino_acid("AUGCCUCCUGUUGUUUGA"))
    

