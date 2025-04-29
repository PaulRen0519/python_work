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

def analyze_mrna(sequence):
    """Analyze mRNA sequence and return most frequent amino acid"""
    sequence = sequence.upper().replace(" ", "")
    start_pos = sequence.find("AUG")
    
    if start_pos == -1:
        return "Error: No start codon (AUG) found"
    
    aa_counts = {}
    i = start_pos
    while i <= len(sequence) - 3:
        codon = sequence[i:i+3]
        if codon in {'UAA', 'UAG', 'UGA'}:
            break
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, '?')
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        i += 3
    
    if not aa_counts:
        return "Error: No valid codons found between start and stop"
    
    most_common = max(aa_counts.items(), key=lambda x: x[1])
    return f"Most frequent amino acid: {most_common[0]} (appeared {most_common[1]} times)"

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
            
        result = analyze_mrna(user_input)
        print(result + "\n")

if __name__ == "__main__":
    main()