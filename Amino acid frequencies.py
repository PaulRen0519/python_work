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

import matplotlib.pyplot as plt
import numpy as np

def plot_amino_acid_frequencies(gene_sequence):
    stop_codons = {'UAA', 'UAG', 'UGA'}
    for i in range(0, len(gene_sequence), 3):
        codon = gene_sequence[i:i+3]
        if codon in stop_codons:
            actual_sequence = gene_sequence[:i]    
            break
    aa_counts = {}
    total_aas = 0
    for i in range(0, len(actual_sequence), 3):
        codon = actual_sequence[i:i+3]
        if len(codon) == 3:
            aa = codon_table[codon]
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
            total_aas += 1
    amino_acids = sorted(aa_counts.keys())
    colors = plt.cm.tab20(np.linspace(0, 1, len(amino_acids)))
    for aa in amino_acids:
        frequency = aa_counts[aa] / total_aas * 100
        plt.bar(aa, frequency, color=colors[amino_acids.index(aa)])
    plt.title('Amino acid frequencies in the gene sequence')
    plt.xlabel('Amino acid')
    plt.ylabel('Frequency (%)')
    plt.show()
plot_amino_acid_frequencies("AUGGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA") 