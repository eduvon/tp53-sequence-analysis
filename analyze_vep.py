import pandas as pd

df = pd.read_csv(
    "tp53_clinvar.vep.txt",
    sep="\t",
    skiprows=32
)

df = df.rename(columns={"#Uploaded_variation": "Uploaded_variation"})

tp53 = df[df["Gene"] == "ENSG00000141510"]

variant_consequences = tp53[["Uploaded_variation", "Consequence"]]

variant_consequences["Consequence"] = (
    variant_consequences["Consequence"].str.split(",")
)

variant_consequences = variant_consequences.explode("Consequence")

variant_consequences = variant_consequences.drop_duplicates()

print(
    variant_consequences
    .groupby("Consequence")["Uploaded_variation"]
    .nunique()
    .sort_values(ascending=False)
)