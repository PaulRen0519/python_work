import matplotlib.pyplot as plt
from collections import defaultdict

# Standard codon to amino acid mapping
CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def analyze_sequence(sequence):
    """Analyze mRNA sequence and return amino acid frequencies"""
    sequence = sequence.upper().replace(" ", "")
    start_pos = sequence.find("AUG")
    
    if start_pos == -1:
        return None, "Error: No start codon (AUG) found"
    
    aa_counts = defaultdict(int)
    i = start_pos
    while i <= len(sequence) - 3:
        codon = sequence[i:i+3]
        if codon in {'UAA', 'UAG', 'UGA'}:
            break
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, '?')
            aa_counts[aa] += 1
        i += 3
    
    if not aa_counts:
        return None, "Error: No valid codons found between start and stop"
    
    return aa_counts, None

def plot_aa_frequencies(aa_counts):
    """Generate a bar plot of amino acid frequencies"""
    # Sort amino acids by frequency (descending)
    sorted_aas = sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)
    aas, counts = zip(*sorted_aas)
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(aas, counts, color='skyblue')
    
    # Add value labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height}', ha='center', va='bottom')
    
    plt.title('Amino Acid Frequency Distribution', fontsize=14)
    plt.xlabel('Amino Acid', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

def main():
    print("mRNA Sequence Analyzer")
    print("Enter an mRNA sequence to analyze (A, U, C, G only)")
    print("Type 'exit' to quit\n")
    
    while True:
        user_input = input("Input sequence: ").strip()
        
        if user_input.lower() == 'exit':
            print("Exiting program...")
            break
        
        if not user_input:
            print("Error: Empty input\n")
            continue
            
        if not all(c.upper() in {'A', 'U', 'C', 'G'} for c in user_input):
            print("Error: Sequence contains invalid characters (only A, U, C, G allowed)\n")
            continue
            
        aa_counts, error = analyze_sequence(user_input)
        if error:
            print(error + "\n")
        else:
            print("\nAmino Acid Counts:")
            for aa, count in sorted(aa_counts.items()):
                print(f"{aa}: {count}")
            print()
            plot_aa_frequencies(aa_counts)

if __name__ == "__main__":
    main()