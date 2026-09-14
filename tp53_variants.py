import subprocess
from collections import Counter
from parse_variant import parse_vcf_record

# Constants

VCF_FILE = "clinvar_20260905.vcf.gz"
TP53_REGION = "17:7668421-7687490"
TP53_START = 7668421
BIN_SIZE = 1000

result = subprocess.run(
    ["tabix", "clinvar_20260905.vcf.gz", "17:7668421-7687490"],
    capture_output=True,
    text=True
)

# Functions

def get_tp53_variants():
    result = subprocess.run(
        ["tabix", VCF_FILE, TP53_REGION],
        capture_output=True,
        text=True
    )

    variants = []

    for line in result.stdout.splitlines():
        variant = parse_vcf_record(line)
        variants.append(variant)

    return variants

def calculate_genomic_bins(variants):
    bins = []

    for variant in variants:
        bin_number = (variant["position"] - TP53_START) // BIN_SIZE
        bins.append(bin_number)

    return Counter(bins)


# Main analysis

variants = get_tp53_variants()

print(f"Number of variants: {len(variants)}")

bin_counts = calculate_genomic_bins(variants)

for bin_number, num_variants in sorted(bin_counts.items()):
    genomic_start = TP53_START + (bin_number * BIN_SIZE)
    genomic_end = genomic_start + BIN_SIZE - 1

    print(
        f"Bin {bin_number}: "
        f"{genomic_start}-{genomic_end} -> "
        f"{num_variants} variants"
    )

