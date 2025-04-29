def most_frequent_amino_acid(codon):  # translate the most common codon into amino acid
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
    
    if codon is None:
        print("Error: No valid codon to translate.")
        return
    
    amino_acid = codon_to_amino_acid.get(codon, None)
    if amino_acid:
        print(f"The most common amino acid is {amino_acid}.")
    else:
        print("Error: Codon not found in translation table.")
