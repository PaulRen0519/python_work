def is_valid_sequence(seq):
    """Check if sequence contains only A, U, C, G"""
    return all(c.upper() in {'A', 'U', 'C', 'G'} for c in seq) # Check if input is valid

def analyze_sequence(seq):
    """Analyze mRNA sequence from AUG to stop codon"""
    start = seq.find("AUG")
    if start == -1: # If AUG is found, return its starting index (starting from 0) And, return -1 if not found.
        return "Error: No start codon (AUG) found"
    
    codons = {}
    for i in range(start, len(seq)-2, 3):
        codon = seq[i:i+3]
        if codon in {'UAA', 'UAG', 'UGA'}:
            break
        codons[codon] = codons.get(codon, 0) + 1
    
    if not codons:
        return "Error: No valid codons found"
    
    most_common = max(codons.items(), key=lambda x: x[1])
    total = sum(codons.values())
    return most_common[0], most_common[1], total

def main():
    print("mRNA Sequence Analyzer")
    print("Analyzes codons from AUG to stop codon (UAA/UAG/UGA)")
    print("Enter 'quit' to exit\n")
    
    while True:
        seq = input("Enter mRNA sequence (A, U, C, G only): ").strip()
        if seq.lower() == 'quit':
            print("Exiting program")
            break
        
        if not seq:
            print("Error: Empty input\n")
            continue
            
        if not is_valid_sequence(seq):
            print("Error: Sequence must contain only A, U, C, G\n")
            continue
            
        result = analyze_sequence(seq)
        
        if isinstance(result, str):
            print(result + "\n")
        else:
            codon, count, total = result
            print(f"Most frequent codon: {codon}")
            print(f"Count: {count}/{total} ({count/total:.1%})\n")

if __name__ == "__main__":
    main()