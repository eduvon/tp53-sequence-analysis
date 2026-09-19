import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from tp53_variants import get_tp53_variants


TP53_GENE = "ENSG00000141510"
ALL_TRANSCRIPT_VEP_FILE = "tp53_clinvar.vep.txt"
MANE_VEP_FILE = "tp53_clinvar_mane.vep.txt"

RESULTS_DIR = Path("results")
FIGURES_DIR = RESULTS_DIR / "figures"

MAIN_CLINICAL_GROUPS = [
    "Pathogenic",
    "Likely pathogenic",
    "Benign",
    "Likely benign",
    "Uncertain",
]

CONSEQUENCE_FLAGS = [
    "has_missense",
    "has_frameshift",
    "has_stop_gained",
    "has_splice_site",
]


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
        skiprows=header_line,
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


def build_consequence_flags(df):
    consequences = df[
        ["Uploaded_variation", "Consequence"]
    ].copy()

    consequences["Consequence"] = (
        consequences["Consequence"].str.split(",")
    )

    consequences = consequences.explode("Consequence")
    consequences = consequences.drop_duplicates()

    flags = (
        consequences
        .groupby("Uploaded_variation")["Consequence"]
        .agg(list)
        .to_frame()
    )

    flags["has_missense"] = flags["Consequence"].apply(
        lambda x: "missense_variant" in x
    )

    flags["has_frameshift"] = flags["Consequence"].apply(
        lambda x: "frameshift_variant" in x
    )

    flags["has_stop_gained"] = flags["Consequence"].apply(
        lambda x: "stop_gained" in x
    )

    flags["has_splice_site"] = flags["Consequence"].apply(
        lambda x: (
            "splice_acceptor_variant" in x
            or "splice_donor_variant" in x
        )
    )

    return flags.reset_index()


def analyze_consequence(table):
    table = table.reindex(
        index=["Pathogenic", "Uncertain"],
        columns=[True, False],
    )

    if table.isna().any().any():
        raise ValueError(
            "Contingency table is missing an expected row or column."
        )

    odds_ratio, p_value = fisher_exact(table)

    a = table.loc["Pathogenic", True]
    b = table.loc["Pathogenic", False]
    c = table.loc["Uncertain", True]
    d = table.loc["Uncertain", False]

    if any(count == 0 for count in [a, b, c, d]):
        raise ValueError(
            "All cells must be greater than zero to calculate "
            "the confidence interval."
        )

    standard_error = math.sqrt(
        (1 / a) +
        (1 / b) +
        (1 / c) +
        (1 / d)
    )

    log_odds_ratio = math.log(odds_ratio)

    lower = math.exp(
        log_odds_ratio - 1.96 * standard_error
    )

    upper = math.exp(
        log_odds_ratio + 1.96 * standard_error
    )

    return odds_ratio, lower, upper, p_value


def main():
    # -------------------------------------------------------------------
    # Load ClinVar data
    # -------------------------------------------------------------------

    variants = get_tp53_variants()
    clinvar = pd.DataFrame(variants)

    clinvar["id"] = clinvar["id"].astype(int)

    clinvar["clinical_group"] = (
        clinvar["clinical_significance"]
        .apply(classify_clinical_significance)
    )

    # -------------------------------------------------------------------
    # Load all-transcript VEP data
    # -------------------------------------------------------------------

    all_vep = load_vep_file(ALL_TRANSCRIPT_VEP_FILE)

    all_tp53 = all_vep[
        all_vep["Gene"] == TP53_GENE
    ].copy()

    all_transcript_flags = build_consequence_flags(
        all_tp53
    )

    all_transcript_analysis = all_transcript_flags.merge(
        clinvar[
            ["id", "clinical_significance", "clinical_group"]
        ],
        left_on="Uploaded_variation",
        right_on="id",
        how="left",
    )

    # -------------------------------------------------------------------
    # Load MANE Select VEP data
    # -------------------------------------------------------------------

    mane_vep = load_vep_file(MANE_VEP_FILE)

    mane_tp53 = mane_vep[
        mane_vep["Gene"] == TP53_GENE
    ].copy()

    mane_tp53 = mane_tp53[
        mane_tp53["Extra"].str.contains(
            "MANE=MANE_Select",
            na=False,
        )
    ].copy()

    mane_flags = build_consequence_flags(
        mane_tp53
    )

    mane_analysis = mane_flags.merge(
        clinvar[
            ["id", "clinical_significance", "clinical_group"]
        ],
        left_on="Uploaded_variation",
        right_on="id",
        how="left",
    )

    # -------------------------------------------------------------------
    # Descriptive analysis
    # -------------------------------------------------------------------

    main_all_transcript_groups = all_transcript_analysis[
        all_transcript_analysis["clinical_group"].isin(
            MAIN_CLINICAL_GROUPS
        )
    ]

    all_transcript_consequence_percent = (
        main_all_transcript_groups
        .groupby("clinical_group")[CONSEQUENCE_FLAGS]
        .mean()
        .mul(100)
        .round(1)
    )

    main_mane_groups = mane_analysis[
        mane_analysis["clinical_group"].isin(
            MAIN_CLINICAL_GROUPS
        )
    ]

    mane_consequence_percent = (
        main_mane_groups
        .groupby("clinical_group")[CONSEQUENCE_FLAGS]
        .mean()
        .mul(100)
        .round(1)
    )

    # -------------------------------------------------------------------
    # Statistical analysis: Pathogenic vs. Uncertain
    # -------------------------------------------------------------------

    stat_groups = mane_analysis[
        mane_analysis["clinical_group"].isin(
            ["Pathogenic", "Uncertain"]
        )
    ]

    results = []

    for flag in CONSEQUENCE_FLAGS:
        table = pd.crosstab(
            stat_groups["clinical_group"],
            stat_groups[flag],
        )

        odds_ratio, lower, upper, p_value = (
            analyze_consequence(table)
        )

        results.append(
            {
                "consequence": flag,
                "odds_ratio": odds_ratio,
                "ci_lower": lower,
                "ci_upper": upper,
                "p_value": p_value,
            }
        )

    results = pd.DataFrame(results)

    results["adjusted_p_value"] = multipletests(
        results["p_value"],
        method="fdr_bh",
    )[1]

    print("All-transcript variants:",
      all_transcript_analysis["Uploaded_variation"].nunique())

    print("MANE variants:",
        mane_analysis["Uploaded_variation"].nunique())

    print("\nStatistical results:")
    print(results[
        [
            "consequence",
            "odds_ratio",
            "ci_lower",
            "ci_upper",
            "adjusted_p_value",
        ]
    ])

    # -------------------------------------------------------------------
    # Save forest plot
    # -------------------------------------------------------------------

    results["label"] = (
        results["consequence"]
        .str.replace("has_", "", regex=False)
        .str.replace("_", " ", regex=False)
        .str.title()
    )

    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    y_positions = range(len(results))

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.errorbar(
        results["odds_ratio"],
        y_positions,
        xerr=[
            results["odds_ratio"] - results["ci_lower"],
            results["ci_upper"] - results["odds_ratio"],
        ],
        fmt="o",
        capsize=4,
    )

    ax.axvline(1, linestyle="--")

    ax.set_xscale("log")
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(results["label"])
    ax.set_xlabel("Odds ratio (Pathogenic vs. Uncertain)")
    ax.set_title(
        "TP53 Variant Consequences and Clinical Significance"
    )

    fig.tight_layout()

    fig.savefig(
        FIGURES_DIR / "tp53_consequence_odds_ratios.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()


if __name__ == "__main__":
    main()