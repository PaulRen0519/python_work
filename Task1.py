import re

def the_most_common_code(mRNA_sequence):  # find the most common codon
    if not re.search(r"^AUG", mRNA_sequence):  # check if sequence starts with start codon 'AUG'
        print("Warning: mRNA sequence does not start with start codon 'AUG'. Proceeding anyway.")
    
    stop_codons = ['UAA', 'UAG', 'UGA']  # list of stop codons
    trinucleotides = []  # store extracted codons
    
    for i in range(0, len(mRNA_sequence) - 2, 3):  # step through sequence by 3
        trinucleotide = mRNA_sequence[i:i+3]
        if trinucleotide in stop_codons:
            break
        trinucleotides.append(trinucleotide)
    
    if not trinucleotides:
        print("Error: No valid codons found in the sequence.")
        return
    
    freq_dict = {}  # count codon frequencies
    for codon in trinucleotides:
        freq_dict[codon] = freq_dict.get(codon, 0) + 1
    
    max_count = max(freq_dict.values())
    most_common_codons = [codon for codon, count in freq_dict.items() if count == max_count]
    final_codon = most_common_codons[-1]  # take the last one if multiple
    
    print(f"The most common codon is {final_codon} with a count of {max_count}.")
    return final_codon
