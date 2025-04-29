import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Define codon table for translation (use U instead of T for mRNA)
codon_table = {
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
    'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W'
}

# Define stop codons
stop_codons = ['UAA', 'UAG', 'UGA']

# Task 1: Find the most frequent codon
def find_most_frequent_codon(mRNA):
    mRNA = mRNA.upper()
    codon_count = {}
    start_index = mRNA.find('AUG')
    if start_index == -1 or len(mRNA[start_index:]) % 3 != 0:
        print("Invalid mRNA: No start codon or incorrect length")
        return None
    for i in range(start_index, len(mRNA) - 2, 3):
        codon = mRNA[i:i + 3]
        if codon in stop_codons:
            break
        if len(codon) != 3 or codon not in codon_table:
            continue
        codon_count[codon] = codon_count.get(codon, 0) + 1
    if not codon_count:
        print("No valid codons found")
        return None
    return max(codon_count, key=codon_count.get)

# Task 2: Translate the most frequent codon to amino acid
def get_amino_acid_from_codon(mRNA):
    most_frequent_codon = find_most_frequent_codon(mRNA)
    if most_frequent_codon is None:
        return None
    return codon_table[most_frequent_codon]

# Task 3: Plot global amino acid frequency distribution
def plot_amino_acid_frequency(mRNA):
    mRNA = mRNA.upper()
    amino_acid_count = {}
    start_index = mRNA.find('AUG')
    if start_index == -1 or len(mRNA[start_index:]) % 3 != 0:
        print("Invalid mRNA for plotting: No start codon or incorrect length")
        return
    for i in range(start_index, len(mRNA) - 2, 3):
        codon = mRNA[i:i + 3]
        if codon in stop_codons:
            break
        if len(codon) != 3 or codon not in codon_table:
            continue
        amino_acid = codon_table[codon]
        amino_acid_count[amino_acid] = amino_acid_count.get(amino_acid, 0) + 1
    if not amino_acid_count:
        print("No amino acids to plot")
        return

    # Academic styling for the bar plot
    plt.figure(figsize=(10, 6))
    sns.barplot(x=list(amino_acid_count.values()), y=list(amino_acid_count.keys()), palette='viridis')
    plt.xlabel('Frequency', fontsize=12)
    plt.ylabel('Amino Acids', fontsize=12)
    plt.title('Global Amino Acid Frequency Distribution', fontsize=14, pad=15)
    plt.grid(True, axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('amino_acid_frequency.png')
    plt.close()

# Task 4: Expression to Structure Analysis
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
        if len(codon) != 3 or codon not in codon_table:
            continue
        coding_sequence += codon

    if len(coding_sequence) < window_size:
        print("Coding sequence too short for window analysis")
        return

    # 1. Codon Usage Bias Heatmap (Sliding Window)
    codons = [coding_sequence[i:i + 3] for i in range(0, len(coding_sequence), 3)]
    unique_codons = sorted(set(codons))
    num_windows = (len(coding_sequence) - window_size) // step_size + 1
    if num_windows <= 0:
        print("No windows available for analysis")
        return
    codon_freq_matrix = np.zeros((len(unique_codons), num_windows))

    for win_idx in range(num_windows):
        start = win_idx * step_size
        window = coding_sequence[start:start + window_size]
        window_codons = [window[i:i + 3] for i in range(0, len(window), 3) if len(window[i:i + 3]) == 3]
        for codon in window_codons:
            if codon in codon_table:
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
    plt.savefig('codon_usage_heatmap.png')
    plt.close()

    # 2. Amino Acid Frequency Heatmap (Sliding Window)
    amino_acids = sorted(set(codon_table.values()) - {'*'})
    aa_freq_matrix = np.zeros((len(amino_acids), num_windows))

    for win_idx in range(num_windows):
        start = win_idx * step_size
        window = coding_sequence[start:start + window_size]
        window_codons = [window[i:i + 3] for i in range(0, len(window), 3) if len(window[i:i + 3]) == 3]
        for codon in window_codons:
            if codon in codon_table:
                aa = codon_table[codon]
                if aa != '*':
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
    plt.savefig('amino_acid_heatmap.png')
    plt.close()

# Main execution
if __name__ == "__main__":
    # Modified mRNA sequence with stop codon
    mrna_sequence = input("Enter mRNA sequence: ")
    
    # Task 1
    most_frequent_codon = find_most_frequent_codon(mrna_sequence)
    print(f"Most frequent codon: {most_frequent_codon}")
    
    # Task 2
    most_frequent_aa = get_amino_acid_from_codon(mrna_sequence)
    print(f"Most frequent amino acid: {most_frequent_aa}")
    
    # Task 3
    plot_amino_acid_frequency(mrna_sequence)
    
    # Task 4
    expression_to_structure_analysis(mrna_sequence, window_size=30, step_size=3)