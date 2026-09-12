sequence = 'ATGCGTACGTTAGCGATCGATCGATGCTAGCTAGCTAACG'

a_count = sequence.count("A")
g_count = sequence.count("G")
c_count = sequence.count("C")
t_count = sequence.count("T")
length = len(sequence)

gc_total = g_count + c_count
gc_percent = (gc_total / length) * 100

print("DNA Sequence Analysis")
print("---------------------")

print(f"Length: {length}")

print(f"A: {a_count}")
print(f"C: {c_count}")
print(f"G: {g_count}")
print(f"T: {t_count}")
print(f"GC%: {gc_percent:.2f}%")



def reverse_complement(sequence):
    complement = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

    reversed_sequence = sequence[::-1]
    result = ""

    for base in reversed_sequence:
        result += complement[base]

    return result

print(reverse_complement("ATC"))