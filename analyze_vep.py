import pandas as pd

from tp53_variants import get_tp53_variants


TP53_GENE = "ENSG00000141510"
ALL_TRANSCRIPT_VEP_FILE = "tp53_clinvar.vep.txt"
MANE_VEP_FILE = "tp53_clinvar_mane.vep.txt"


def load_vep_file(filename):
    with open(filename) as file:
        for i, line in enumerate(file):
            if line.startswith("#Uploaded_variation"):
                header_line = i
                break
        else:
            raise ValueError(
                f"Could not find VEP header in {filename}"
            )

    df = pd.read_csv(
        filename,
        sep="\t",
        skiprows=header_line
    )

    return df.rename(
        columns={"#Uploaded_variation": "Uploaded_variation"}
    )


def classify_clinical_significance(value):
    if value == "Pathogenic":
        return "Pathogenic"
    elif value == "Likely_pathogenic":
        return "Likely pathogenic"
    elif value == "Benign":
        return "Benign"
    elif value == "Likely_benign":
        return "Likely benign"
    elif value == "Uncertain_significance":
        return "Uncertain"
    else:
        return "Other"


# -------------------------------------------------------------------
# Load ClinVar data
# -------------------------------------------------------------------

variants = get_tp53_variants()
clinvar = pd.DataFrame(variants)

clinvar["id"] = clinvar["id"].astype(int)


# -------------------------------------------------------------------
# Load all-transcript VEP data
# -------------------------------------------------------------------

all_vep = load_vep_file(ALL_TRANSCRIPT_VEP_FILE)

all_tp53 = all_vep[
    all_vep["Gene"] == TP53_GENE
].copy()


# Create one row per variant-consequence pair
variant_consequences = all_tp53[
    ["Uploaded_variation", "Consequence"]
].copy()

variant_consequences["Consequence"] = (
    variant_consequences["Consequence"].str.split(",")
)

variant_consequences = variant_consequences.explode(
    "Consequence"
)

variant_consequences = variant_consequences.drop_duplicates()


# Create variant-level consequence flags
all_transcript_flags = (
    variant_consequences
    .groupby("Uploaded_variation")["Consequence"]
    .agg(list)
    .to_frame()
)

all_transcript_flags["has_missense"] = (
    all_transcript_flags["Consequence"].apply(
        lambda x: "missense_variant" in x
    )
)

all_transcript_flags["has_frameshift"] = (
    all_transcript_flags["Consequence"].apply(
        lambda x: "frameshift_variant" in x
    )
)

all_transcript_flags["has_stop_gained"] = (
    all_transcript_flags["Consequence"].apply(
        lambda x: "stop_gained" in x
    )
)

all_transcript_flags["has_splice_site"] = (
    all_transcript_flags["Consequence"].apply(
        lambda x: (
            "splice_acceptor_variant" in x
            or "splice_donor_variant" in x
        )
    )
)

all_transcript_flags = all_transcript_flags.reset_index()


# Add ClinVar information
all_transcript_analysis = all_transcript_flags.merge(
    clinvar[["id", "clinical_significance"]],
    left_on="Uploaded_variation",
    right_on="id",
    how="left"
)

all_transcript_analysis["clinical_group"] = (
    all_transcript_analysis["clinical_significance"]
    .apply(classify_clinical_significance)
)


# -------------------------------------------------------------------
# MANE Select VEP data
# -------------------------------------------------------------------

mane_vep = load_vep_file(MANE_VEP_FILE)

mane_tp53 = mane_vep[
    mane_vep["Gene"] == TP53_GENE
].copy()

mane_tp53 = mane_tp53[
    mane_tp53["Extra"].str.contains(
        "MANE=MANE_Select",
        na=False
    )
].copy()


# Create MANE consequence flags
mane_flags = mane_tp53[
    ["Uploaded_variation", "Consequence"]
].copy()

mane_flags["has_missense"] = (
    mane_flags["Consequence"].apply(
        lambda x: "missense_variant" in x.split(",")
    )
)

mane_flags["has_frameshift"] = (
    mane_flags["Consequence"].apply(
        lambda x: "frameshift_variant" in x.split(",")
    )
)

mane_flags["has_stop_gained"] = (
    mane_flags["Consequence"].apply(
        lambda x: "stop_gained" in x.split(",")
    )
)

mane_flags["has_splice_site"] = (
    mane_flags["Consequence"].apply(
        lambda x: (
            "splice_acceptor_variant" in x.split(",")
            or "splice_donor_variant" in x.split(",")
        )
    )
)


# Add ClinVar information
mane_analysis = mane_flags.merge(
    clinvar[["id", "clinical_significance"]],
    left_on="Uploaded_variation",
    right_on="id",
    how="left"
)

mane_analysis["clinical_group"] = (
    mane_analysis["clinical_significance"]
    .apply(classify_clinical_significance)
)


# -------------------------------------------------------------------
# Descriptive comparison
# -------------------------------------------------------------------

consequence_flags = [
    "has_missense",
    "has_frameshift",
    "has_stop_gained",
    "has_splice_site"
]


main_all_transcript_groups = all_transcript_analysis[
    all_transcript_analysis["clinical_group"] != "Other"
]

all_transcript_consequence_percent = (
    main_all_transcript_groups
    .groupby("clinical_group")[consequence_flags]
    .mean()
    .mul(100)
    .round(1)
)


main_mane_groups = mane_analysis[
    mane_analysis["clinical_group"] != "Other"
]

mane_consequence_percent = (
    main_mane_groups
    .groupby("clinical_group")[consequence_flags]
    .mean()
    .mul(100)
    .round(1)
)


# -------------------------------------------------------------------
# Check that the two analyses contain the expected number of variants
# -------------------------------------------------------------------

print("All-transcript variants:",
      all_transcript_analysis["Uploaded_variation"].nunique())

print("MANE variants:",
      mane_analysis["Uploaded_variation"].nunique())

print("\nAll-transcript analysis:")
print(all_transcript_consequence_percent)

print("\nMANE analysis:")
print(mane_consequence_percent)