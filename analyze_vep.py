import pandas as pd

from tp53_variants import get_tp53_variants


# Load ClinVar data
variants = get_tp53_variants()
clinvar = pd.DataFrame(variants)
clinvar["id"] = clinvar["id"].astype(int)


# Load VEP data
df = pd.read_csv(
    "tp53_clinvar.vep.txt",
    sep="\t",
    skiprows=32
)

df = df.rename(columns={"#Uploaded_variation": "Uploaded_variation"})


# Keep only TP53 annotations
tp53 = df[df["Gene"] == "ENSG00000141510"]


# Create one row per variant-consequence pair
variant_consequences = tp53[["Uploaded_variation", "Consequence"]].copy()

variant_consequences["Consequence"] = (
    variant_consequences["Consequence"].str.split(",")
)

variant_consequences = variant_consequences.explode("Consequence")
variant_consequences = variant_consequences.drop_duplicates()


# Add ClinVar clinical significance
merged = variant_consequences.merge(
    clinvar[["id", "clinical_significance"]],
    left_on="Uploaded_variation",
    right_on="id",
    how="left"
)


# Create variant-level consequence flags
variant_flags = (
    variant_consequences
    .groupby("Uploaded_variation")["Consequence"]
    .agg(list)
    .to_frame()
)

variant_flags["has_missense"] = variant_flags["Consequence"].apply(
    lambda x: "missense_variant" in x
)

variant_flags["has_frameshift"] = variant_flags["Consequence"].apply(
    lambda x: "frameshift_variant" in x
)

variant_flags["has_stop_gained"] = variant_flags["Consequence"].apply(
    lambda x: "stop_gained" in x
)

variant_flags["has_splice_site"] = variant_flags["Consequence"].apply(
    lambda x: (
        "splice_acceptor_variant" in x
        or "splice_donor_variant" in x
    )
)

variant_flags = variant_flags.reset_index()


# Add ClinVar clinical significance to variant-level data
variant_analysis = variant_flags.merge(
    clinvar[["id", "clinical_significance"]],
    left_on="Uploaded_variation",
    right_on="id",
    how="left"
)


# Inspect the final analysis table
print(
    variant_analysis[
        [
            "Uploaded_variation",
            "clinical_significance",
            "has_missense",
            "has_frameshift",
            "has_stop_gained",
            "has_splice_site"
        ]
    ].head(20)
)