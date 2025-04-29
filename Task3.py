import matplotlib.pyplot as plt
import time

def plot_amino_acid_frequencies(mRNA_sequence):  # plot amino acid frequencies from mRNA sequence
    stop_codons = ['UAA', 'UAG', 'UGA']  # list of stop codons
    codon_to_amino_acid = {
        'AUG': 'Methionine',
        'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
        'UUA': 'Leucine', 'UUG': 'Leucine',
        'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
        'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine',
        'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
        'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
        'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
        'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
        'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
        'UAU': 'Tyrosine', 'UAC': 'Tyrosine',
        'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
        'CAU': 'Histidine', 'CAC': 'Histidine',
        'CAA': 'Glutamine', 'CAG': 'Glutamine',
        'AAU': 'Asparagine', 'AAC': 'Asparagine',
        'AAA': 'Lysine', 'AAG': 'Lysine',
        'GAU': 'Aspartic Acid', 'GAC': 'Aspartic Acid',
        'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
        'UGU': 'Cysteine', 'UGC': 'Cysteine',
        'UGG': 'Tryptophan',
        'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
        'AGU': 'Serine', 'AGC': 'Serine',
        'AGA': 'Arginine', 'AGG': 'Arginine',
        'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine'
    }
    
    trinucleotides = []  # store extracted codons
    for i in range(0, len(mRNA_sequence) - 2, 3):
        trinucleotide = mRNA_sequence[i:i+3]
        if trinucleotide in stop_codons:
            break
        trinucleotides.append(trinucleotide)
    
    if not trinucleotides:
        print("Error: No valid codons found in sequence.")
        return
    
    amino_acids = []  # list to store translated amino acids
    for codon in trinucleotides:
        amino_acid = codon_to_amino_acid.get(codon, None)
        if amino_acid and amino_acid != 'Stop':
            amino_acids.append(amino_acid)
    
    if not amino_acids:
        print("Error: No amino acids found after translation.")
        return
    
    freq_dict = {}  # count amino acid frequencies
    for aa in amino_acids:
        freq_dict[aa] = freq_dict.get(aa, 0) + 1
    
    # Prepare data for plotting
    amino_acid_list = list(freq_dict.keys())
    frequency_list = list(freq_dict.values())
    
    start_time = time.time()  # record start time
    
    plt.figure(figsize=(12, 6))  # set figure size
    plt.bar(amino_acid_list, frequency_list)  # plot the bar chart
    plt.title('Frequency of Amino Acids')  # set the title
    plt.xlabel('Amino Acid')  # set the x-axis label
    plt.ylabel('Frequency')  # set the y-axis label
    plt.xticks(rotation=45)  # rotate x-axis labels for better readability
    plt.tight_layout()  # adjust layout
    plt.show()  # display the plot
    
    end_time = time.time()  # record end time
    print(f"Time taken for plotting: {end_time - start_time:.4f} seconds")  # print the time taken
