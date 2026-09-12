def parse_fasta(filename):
    sequence = ""

    with open(filename) as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip()
    
    return sequence


# sequence = parse_fasta("tp53.fasta")

# print(f"Sequence length: {len(sequence)}")

def translate(sequence):
    codon_table = {
        "ATG": "M",
        "GAG": "E",
        "CCG": "P",
        "CAG": "Q"
    }

    protein = ""

    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        amino_acid = codon_table[codon]
        protein += amino_acid
    return protein


print(translate("ATGGAGGAGCCGCAG"))