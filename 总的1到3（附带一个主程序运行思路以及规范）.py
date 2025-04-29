import re
import matplotlib.pyplot as plt
import time
import seaborn as sns
import numpy as np

# Define the codon to amino acid translation table
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

stop_codons = ['UAA', 'UAG', 'UGA']

def the_most_common_code(mRNA_sequence):  # find the most common codon
    if not re.search(r"^AUG", mRNA_sequence):
        print("Warning: mRNA sequence does not start with start codon 'AUG'. Proceeding anyway.")
    
    stop_codons = ['UAA', 'UAG', 'UGA']
    trinucleotides = []

    start_index = mRNA_sequence.find('AUG')
    if start_index == -1:
        print("Error: No start codon (AUG) found.")
        return None

    for i in range(start_index, len(mRNA_sequence) - 2, 3):
        trinucleotide = mRNA_sequence[i:i+3]
        if trinucleotide in stop_codons:
            break
        trinucleotides.append(trinucleotide)
    
    if not trinucleotides:
        print("Error: No valid codons found.")
        return None
    
    freq_dict = {}
    for codon in trinucleotides:
        freq_dict[codon] = freq_dict.get(codon, 0) + 1
    
    max_count = max(freq_dict.values())
    most_common_codons = [codon for codon, count in freq_dict.items() if count == max_count]
    
    print(f"The most common codon(s): {', '.join(most_common_codons)} with a count of {max_count}.")
    return most_common_codons

def most_frequent_amino_acid(codons):  # translate the most common codon(s)
    if codons is None:
        print("Error: No valid codon to translate.")
        return None
    
    amino_acids = []
    for codon in codons:
        amino_acid = codon_to_amino_acid.get(codon, None)
        if amino_acid:
            amino_acids.append(amino_acid)
        else:
            print(f"Error: Codon {codon} not found in translation table.")
    
    if amino_acids:
        print(f"The most common amino acid(s): {', '.join(amino_acids)}.")
    else:
        print("Error: No valid amino acid translation found.")
    return amino_acids

def plot_amino_acid_frequencies(mRNA_sequence):  # plot amino acid frequencies
    stop_codons = ['UAA', 'UAG', 'UGA']
    trinucleotides = []

    start_index = mRNA_sequence.find('AUG')
    if start_index == -1:
        print("Error: No start codon (AUG) found.")
        return

    for i in range(start_index, len(mRNA_sequence) - 2, 3):
        trinucleotide = mRNA_sequence[i:i+3]
        if trinucleotide in stop_codons:
            break
        trinucleotides.append(trinucleotide)
    
    if not trinucleotides:
        print("Error: No valid codons found in sequence.")
        return
    
    amino_acids = []
    for codon in trinucleotides:
        amino_acid = codon_to_amino_acid.get(codon, None)
        if amino_acid and amino_acid != 'Stop':
            amino_acids.append(amino_acid)
    
    if not amino_acids:
        print("Error: No amino acids found after translation.")
        return
    
    freq_dict = {}
    for aa in amino_acids:
        freq_dict[aa] = freq_dict.get(aa, 0) + 1
    
    amino_acid_list = list(freq_dict.keys())
    frequency_list = list(freq_dict.values())
    
    start_time = time.time()
    
    plt.figure(figsize=(12, 6))
    plt.bar(amino_acid_list, frequency_list)
    plt.title('Frequency of Amino Acids')
    plt.xlabel('Amino Acid')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
    
    end_time = time.time()
    print(f"Time taken for plotting: {end_time - start_time:.4f} seconds")

# Task 4
def expression_to_structure_analysis(mRNA, window_size=30, step_size=3):
    mRNA = mRNA.upper()
    start_index = mRNA.find('AUG')
    if start_index == -1 or len(mRNA[start_index:]) % 3 != 0:
        print("Invalid mRNA for analysis: No start codon or incorrect length")
        return

    # Extract the coding sequence (from AUG to stop codon)
    coding_sequence = ''
    for i in range(start_index, len(mRNA) - 2, 3):
        codon = mRNA[i:i + 3]
        if codon in stop_codons:
            break
        if len(codon) != 3 or codon not in codon_to_amino_acid:
            print(f"Invalid codon {codon} at position {i}")
            continue
        coding_sequence += codon
    print(f"Coding sequence: {coding_sequence}, Length: {len(coding_sequence)}")

    if len(coding_sequence) < window_size:
        print(f"Coding sequence too short for window analysis (length: {len(coding_sequence)}, window_size: {window_size})")
        return

    # 1. Codon Usage Bias Heatmap (Sliding Window)
    codons = [coding_sequence[i:i + 3] for i in range(0, len(coding_sequence), 3)]
    unique_codons = sorted(set(codons))
    num_windows = (len(coding_sequence) - window_size) // step_size + 1
    if num_windows <= 0:
        print(f"No windows available for analysis (sequence length: {len(coding_sequence)}, window_size: {window_size}, step_size: {step_size})")
        return
    print(f"Number of windows: {num_windows}, Unique codons: {unique_codons}")

    # Initialize codon frequency matrix
    codon_freq_matrix = np.zeros((len(unique_codons), num_windows))

    for win_idx in range(num_windows):
        start = win_idx * step_size
        window = coding_sequence[start:start + window_size]
        window_codons = [window[i:i + 3] for i in range(0, len(window), 3) if len(window[i:i + 3]) == 3]
        for codon in window_codons:
            if codon in codon_to_amino_acid:
                codon_idx = unique_codons.index(codon)
                codon_freq_matrix[codon_idx, win_idx] += 1

    # Normalize frequencies per window
    codon_freq_matrix = codon_freq_matrix / (window_size // 3)

    # Plot Codon Usage Heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(codon_freq_matrix, cmap='YlOrRd', xticklabels=range(num_windows), yticklabels=unique_codons)
    plt.xlabel('Window Position (Step Size = 3 Nucleotides)', fontsize=12)
    plt.ylabel('Codons', fontsize=12)
    plt.title('Codon Usage Bias Across mRNA Sequence', fontsize=14, pad=15)
    plt.tight_layout()
    try:
        plt.savefig('codon_usage_heatmap.png')
        print("Saved codon_usage_heatmap.png")
    except Exception as e:
        print(f"Error saving codon_usage_heatmap.png: {e}")
    plt.close()
    # 2. Amino Acid Frequency Heatmap (Sliding Window)
    amino_acids = sorted(set(codon_to_amino_acid.values()) - {'Stop'})  # Fixed to match 'Stop'
    aa_freq_matrix = np.zeros((len(amino_acids), num_windows))

    for win_idx in range(num_windows):
        start = win_idx * step_size
        window = coding_sequence[start:start + window_size]
        window_codons = [window[i:i + 3] for i in range(0, len(window), 3) if len(window[i:i + 3]) == 3]
        for codon in window_codons:
            if codon in codon_to_amino_acid:
                aa = codon_to_amino_acid[codon]
                if aa != 'Stop':
                    aa_idx = amino_acids.index(aa)
                    aa_freq_matrix[aa_idx, win_idx] += 1

    # Normalize frequencies per window
    aa_freq_matrix = aa_freq_matrix / (window_size // 3)

    # Plot Amino Acid Frequency Heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(aa_freq_matrix, cmap='Blues', xticklabels=range(num_windows), yticklabels=amino_acids)
    plt.xlabel('Window Position (Step Size = 3 Nucleotides)', fontsize=12)
    plt.ylabel('Amino Acids', fontsize=12)
    plt.title('Local Amino Acid Frequency Across Translated Sequence', fontsize=14, pad=15)
    plt.tight_layout()
    try:
        plt.savefig('amino_acid_heatmap.png')
        print("Saved amino_acid_heatmap.png")
    except Exception as e:
        print(f"Error saving amino_acid_heatmap.png: {e}")
    plt.close()

# --- Main Program Start ---
if __name__ == "__main__":
    mRNA_sequence = input("Enter the mRNA sequence: ").strip().upper()
    
    if not mRNA_sequence:
        print("Error: Empty input.")
    else:
        final_codons = the_most_common_code(mRNA_sequence)
        most_frequent_amino_acid(final_codons)
        plot_amino_acid_frequencies(mRNA_sequence)
        expression_to_structure_analysis(mRNA_sequence, window_size=30, step_size=3)