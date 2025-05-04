import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# codon translation table: translate codon into amino acid
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

# stop codon list
stop_codons = ['UAA', 'UAG', 'UGA']

# transcript DNA into mRNA
def dna_to_mrna(dna):
    """
    transcript dna into mrna
    parameter: dna - input DNA sequence
    return: mrna 
    rule: A->U, T->A, C->G, G->C
    """
    dna = dna.upper()
    mrna = ''
    for base in dna:
        if base == 'A':
            mrna += 'U'
        elif base == 'T':
            mrna += 'A'
        elif base == 'C':
            mrna += 'G'
        elif base == 'G':
            mrna += 'C'
        else:
            raise ValueError("Invalid bases! Only available for A, T, C, G.")
    return mrna
        
                    
# Function 1: Find the most common codon(s)
def the_most_common_code(mRNA_sequence):
    """
    find the most frequent codon
    parameter: mRNA_sequence
    return: the most frequent codon list
    """
    # check the sequence start with AUG
    if not re.search(r"^AUG", mRNA_sequence):
        print("Warning: mRNA sequence does not start with start codon 'AUG'. Proceeding anyway.")
    
    trinucleotides = []
    start_index = mRNA_sequence.find('AUG')
    
    # begining with AUG, extract one codon every three bases
    for i in range(start_index, len(mRNA_sequence) - 2, 3):
        trinucleotide = mRNA_sequence[i:i+3]
        if trinucleotide in stop_codons:                    # stop when meeting stop codon
            break
        trinucleotides.append(trinucleotide)
    if not trinucleotides:
        print("Error: No valid codons found.")
        return None
    
    # collect the frequencies of every codon
    freq_dict = {}
    for codon in trinucleotides:
        freq_dict[codon] = freq_dict.get(codon, 0) + 1
    max_count = max(freq_dict.values())
    most_common_codons = [codon for codon, count in freq_dict.items() if count == max_count]
    
    print(f"The most common codon(s): {', '.join(most_common_codons)} with a count of {max_count}.")
    return most_common_codons

# Function 2: Translate most common codon(s)
def most_frequent_amino_acid(mRNA_sequence):
    """
    find the most frequent amino acid
    parameter: mRNA_sequence
    return: the most frequent amino acid list
    """
    mRNA_sequence = mRNA_sequence.upper().replace(" ", "")
    start_pos = mRNA_sequence.find("AUG")
    
    # check the start codon
    if start_pos == -1:
        return "Error: No start codon (AUG) found"
    
    # read each codon, translate them into amino acid, and count the frequencies of each amino acid
    aa_counts = {}
    i = start_pos
    while i <= len(mRNA_sequence) - 3:
        codon = mRNA_sequence[i:i+3]
        if codon in {'UAA', 'UAG', 'UGA'}:
            break
        if len(codon) == 3:
            aa = codon_to_amino_acid.get(codon, '?')
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        i += 3
    
    if not aa_counts:
        return "Error: No valid codons found between start and stop"
    
    # find the most frequent amino acid
    most_common = max(aa_counts.items(), key=lambda x: x[1])
    return f"Most frequent amino acid: {most_common[0]} (appeared {most_common[1]} times)"


# Function 3: Plot amino acid frequencies
def plot_amino_acid_frequencies(mRNA_sequence):
    """
    plot the amino acid frequencey figure
    parameter: mRNA_sequence
    return: plot figure
    """

    # check the start codon
    start_index = mRNA_sequence.find('AUG')
    if start_index == -1:
        print("Error: No start codon (AUG) found.")
        return
    
    # read the codon every three bases
    trinucleotides = []
    for i in range(start_index, len(mRNA_sequence) - 2, 3):
        trinucleotide = mRNA_sequence[i:i+3]
        if trinucleotide in stop_codons:
            break
        trinucleotides.append(trinucleotide)
    if not trinucleotides:
        print("Error: No valid codons found in sequence.")
        return
    
    # translate the codons into amino acids
    amino_acids = []
    for codon in trinucleotides:
        amino_acid = codon_to_amino_acid.get(codon, None)
        if amino_acid and amino_acid != 'Stop':
            amino_acids.append(amino_acid)
    if not amino_acids:
        print("Error: No amino acids found after translation.")
        return
    
    # record the frequencies of each amino acids
    freq_dict = {}
    for aa in amino_acids:
        freq_dict[aa] = freq_dict.get(aa, 0) + 1
    amino_acid_list = list(freq_dict.keys())
    frequency_list = list(freq_dict.values())
    
    # plot the figure
    plt.figure(figsize=(12, 6))
    plt.bar(amino_acid_list, frequency_list)
    plt.title('Frequency of Amino Acids')
    plt.xlabel('Amino Acid')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

# Function 4: Expression to Structure Analysis (Additional Function)
def expression_to_structure_analysis(mRNA_sequence, window_size=30, step_size=3):
    """
    plot two heatmaps for codons and amino acids
    parameters: mRNA_sequence, window_size (limit the analyzing range), step_size (read codon)
    return: two heatmaps, one for codons, another for amino acids
    """
    # check the start codon
    mRNA_sequence = mRNA_sequence.upper()
    start_index = mRNA_sequence.find('AUG')
    if start_index == -1:
        print("Invalid mRNA for analysis: No start codon.")
        return
    
    # read the codon
    coding_sequence = ''
    for i in range(start_index, len(mRNA_sequence) - 2, 3):
        codon = mRNA_sequence[i:i+3]
        if codon in stop_codons:
            break
        if len(codon) != 3 or codon not in codon_to_amino_acid:
            continue
        coding_sequence += codon
    
    # check the length of sequence, ensuring the sequence is longer than window size
    if len(coding_sequence) < window_size:
        print("Coding sequence too short for window analysis.")
        return
    codons = [coding_sequence[i:i+3] for i in range(0, len(coding_sequence), 3)]
    unique_codons = sorted(set(codons))
    num_windows = (len(coding_sequence) - window_size) // step_size + 1     # calculate the number of the window
    if num_windows <= 0:
        print("No windows available for analysis.")
        return
    
    codon_freq_matrix = np.zeros((len(unique_codons), num_windows))
    for win_idx in range(num_windows):
        start = win_idx * step_size
        window = coding_sequence[start:start + window_size]
        window_codons = [window[i:i+3] for i in range(0, len(window), 3) if len(window[i:i+3]) == 3]
        
        # statistic the frequencies of each amino acid in thee window
        for codon in window_codons:
            if codon in codon_to_amino_acid:
                codon_idx = unique_codons.index(codon)
                codon_freq_matrix[codon_idx, win_idx] += 1
    
    codon_freq_matrix = codon_freq_matrix / (window_size // 3)              # standardize the frequencies
    
    # plot the heatmap, using a spectrum of color to beautify the figure
    plt.figure(figsize=(12, 8))
    sns.heatmap(codon_freq_matrix, cmap='YlOrRd', xticklabels=range(num_windows), yticklabels=unique_codons)
    plt.xlabel('Window Position (Step Size = 3 Nucleotides)', fontsize=12)
    plt.ylabel('Codons', fontsize=12)
    plt.title('Codon Usage Bias Across mRNA Sequence', fontsize=14, pad=15)
    plt.tight_layout()
    plt.show()
    
    # translate the codon into amino acids
    amino_acids = sorted(set(codon_to_amino_acid.values()) - {'Stop'})
    aa_freq_matrix = np.zeros((len(amino_acids), num_windows))
    
    # statistic the frequencies of the amino acid
    for win_idx in range(num_windows):
        start = win_idx * step_size
        window = coding_sequence[start:start + window_size]
        window_codons = [window[i:i+3] for i in range(0, len(window), 3) if len(window[i:i+3]) == 3]
        for codon in window_codons:
            aa = codon_to_amino_acid.get(codon, None)
            if aa and aa != 'Stop':
                aa_idx = amino_acids.index(aa)
                aa_freq_matrix[aa_idx, win_idx] += 1
    aa_freq_matrix = aa_freq_matrix / (window_size // 3)
    plt.figure(figsize=(12, 8))
    sns.heatmap(aa_freq_matrix, cmap='Blues', xticklabels=range(num_windows), yticklabels=amino_acids)
    plt.xlabel('Window Position (Step Size = 3 Nucleotides)', fontsize=12)
    plt.ylabel('Amino Acids', fontsize=12)
    plt.title('Local Amino Acid Frequency Across Translated Sequence', fontsize=14, pad=15)
    plt.tight_layout()
    plt.show()

# --- Main Program ---
if __name__ == "__main__":

    #  main programe enter
    dna = input("Enter the DNA sequence: ").strip().upper()
    mRNA_sequence = dna_to_mrna(dna)
    if not mRNA_sequence:
        print("Error: Empty input.")
    else:
        # execute four analyzing function
        the_most_common_code(mRNA_sequence)                 # the most frequent codon
        most_frequent_amino_acid(mRNA_sequence)             # the most frequent amino acid
        plot_amino_acid_frequencies(mRNA_sequence)          # visualize the frequencies of amino acid
        expression_to_structure_analysis(mRNA_sequence)     # structural expressing the distribution of codons and amino acids through heatmaps
