from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# STEP 24 - TECHNOLOGY DECISION MATRIX
# ============================================================

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

RESULTS.mkdir(parents=True, exist_ok=True)

INPUT_FILE = RESULTS / "integrated_chemistry_summary.csv"

CHEMISTRIES = ["LFP", "NMC811"]


# ============================================================
# 1. LOAD EXPERIMENTAL DATA
# ============================================================

print("=" * 75)
print("STEP 24 - TECHNOLOGY DECISION MATRIX")
print("=" * 75)

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Required input file not found:\n{INPUT_FILE}\n"
        "Run Step 22 first."
    )

df = pd.read_csv(INPUT_FILE)

df["Chemistry"] = df["Chemistry"].astype(str).str.strip()

if df["Chemistry"].duplicated().any():
    raise ValueError("Duplicate chemistry entries detected.")

df = df.set_index("Chemistry")

if set(df.index) != set(CHEMISTRIES):
    raise ValueError(
        f"Expected chemistries: {CHEMISTRIES}. "
        f"Found: {list(df.index)}"
    )

df = df.loc[CHEMISTRIES]

print("\nExperimental data loaded successfully.")


# ============================================================
# 2. DEFINE CRITERIA
# ============================================================

# Experimental values are taken directly from Step 22.
#
# Energy-density values are taken from the reference paper.
# They describe the specific tested cells, not complete packs.

criteria = pd.DataFrame(
    index=CHEMISTRIES,
    data={
        "TR propagation [s]":
            df["Total TR propagation [s]"],

        "Cell Tmax [C]":
            df["Mean cell-associated Tmax [C]"],

        "Venting Tmax [C]":
            df["Mean venting Tmax [C]"],

        "Gravimetric energy density [Wh/kg]":
            [168.0, 248.0],

        "Volumetric energy density [Wh/L]":
            [361.0, 569.0],
    },
)

# True means higher values are preferred.
# False means lower values are preferred.

BENEFIT_CRITERIA = {
    "TR propagation [s]": True,
    "Cell Tmax [C]": False,
    "Venting Tmax [C]": False,
    "Gravimetric energy density [Wh/kg]": True,
    "Volumetric energy density [Wh/L]": True,
}

if criteria.isna().any().any():
    raise ValueError("Missing criterion values detected.")

if not np.isfinite(criteria.to_numpy(dtype=float)).all():
    raise ValueError("Non-finite criterion values detected.")

if (criteria <= 0).any().any():
    raise ValueError(
        "All criteria must be positive for this normalization."
    )

print("\nRAW CRITERIA")
print("-" * 75)
print(criteria.round(2).to_string())


# ============================================================
# 3. NORMALIZE CRITERIA
# ============================================================

normalized = pd.DataFrame(index=CHEMISTRIES)

for criterion in criteria.columns:

    values = criteria[criterion].astype(float)

    if BENEFIT_CRITERIA[criterion]:

        normalized[criterion] = (
            values / values.max()
        )

    else:

        normalized[criterion] = (
            values.min() / values
        )

print("\nNORMALIZED CRITERIA")
print("-" * 75)
print(normalized.round(4).to_string())


# ============================================================
# 4. DEFINE ENGINEERING SCENARIOS
# ============================================================

WEIGHTS = {
    "Safety-oriented": {
        "TR propagation [s]": 0.35,
        "Cell Tmax [C]": 0.20,
        "Venting Tmax [C]": 0.15,
        "Gravimetric energy density [Wh/kg]": 0.15,
        "Volumetric energy density [Wh/L]": 0.15,
    },

    "Energy-density-oriented": {
        "TR propagation [s]": 0.10,
        "Cell Tmax [C]": 0.10,
        "Venting Tmax [C]": 0.10,
        "Gravimetric energy density [Wh/kg]": 0.35,
        "Volumetric energy density [Wh/L]": 0.35,
    },
}

weights_df = pd.DataFrame(WEIGHTS)

if set(weights_df.index) != set(criteria.columns):
    raise ValueError("Weight criteria do not match input criteria.")

for scenario in weights_df.columns:

    total = weights_df[scenario].sum()

    if not np.isclose(total, 1.0):
        raise ValueError(
            f"Weights for {scenario} sum to {total}, not 1."
        )


# ============================================================
# 5. CALCULATE WEIGHTED SCORES
# ============================================================

scores = pd.DataFrame(index=CHEMISTRIES)

contributions = []

for scenario in weights_df.columns:

    weights = weights_df[scenario]

    weighted = normalized.mul(weights, axis=1)

    scores[scenario] = weighted.sum(axis=1)

    for chemistry in CHEMISTRIES:

        for criterion in criteria.columns:

            contributions.append({
                "Scenario": scenario,
                "Chemistry": chemistry,
                "Criterion": criterion,
                "Normalized score":
                    normalized.loc[chemistry, criterion],
                "Weight": weights[criterion],
                "Weighted contribution":
                    weighted.loc[chemistry, criterion],
            })

contributions_df = pd.DataFrame(contributions)

print("\nWEIGHTED SCORES")
print("-" * 75)
print(scores.round(4).to_string())

print("\nSCENARIO WINNERS")
print("-" * 75)

for scenario in scores.columns:

    winner = scores[scenario].idxmax()

    print(
        f"{scenario}: {winner} "
        f"(score = {scores.loc[winner, scenario]:.4f})"
    )


# ============================================================
# 6. SENSITIVITY ANALYSIS
# ============================================================

# Interpolate between the two weight vectors.
#
# alpha = 0: Safety-oriented
# alpha = 1: Energy-density-oriented
#
# This tests sensitivity to the assumed engineering priorities.

alphas = np.linspace(0.0, 1.0, 101)

sensitivity_rows = []

safety_weights = weights_df["Safety-oriented"]
energy_weights = weights_df["Energy-density-oriented"]

for alpha in alphas:

    weights = (
        (1.0 - alpha) * safety_weights
        + alpha * energy_weights
    )

    scenario_scores = normalized.mul(
        weights,
        axis=1
    ).sum(axis=1)

    sensitivity_rows.append({
        "Alpha": alpha,
        "LFP score": scenario_scores.loc["LFP"],
        "NMC811 score": scenario_scores.loc["NMC811"],
        "Score difference LFP-NMC811":
            scenario_scores.loc["LFP"]
            - scenario_scores.loc["NMC811"],
    })

sensitivity = pd.DataFrame(sensitivity_rows)

difference = sensitivity[
    "Score difference LFP-NMC811"
].to_numpy()

switching_alpha = None

# Exact linear interpolation of the score difference.
d0 = difference[0]
d1 = difference[-1]

if np.isclose(d0, 0.0):
    switching_alpha = 0.0

elif np.isclose(d1, 0.0):
    switching_alpha = 1.0

elif d0 * d1 < 0:

    switching_alpha = d0 / (d0 - d1)

print("\nSENSITIVITY ANALYSIS")
print("-" * 75)

if switching_alpha is not None:

    print(
        "Decision switching point: "
        f"alpha = {switching_alpha:.4f}"
    )

else:

    print(
        "No ranking reversal occurs between "
        "the two selected weight scenarios."
    )


# ============================================================
# 7. SAVE TABLES
# ============================================================

criteria.to_csv(
    RESULTS / "technology_raw_criteria.csv",
    index_label="Chemistry"
)

normalized.to_csv(
    RESULTS / "technology_normalized_criteria.csv",
    index_label="Chemistry"
)

weights_df.to_csv(
    RESULTS / "technology_scenario_weights.csv",
    index_label="Criterion"
)

scores.to_csv(
    RESULTS / "technology_decision_scores.csv",
    index_label="Chemistry"
)

contributions_df.to_csv(
    RESULTS / "technology_weighted_contributions.csv",
    index=False
)

sensitivity.to_csv(
    RESULTS / "technology_sensitivity_analysis.csv",
    index=False
)


# ============================================================
# 8. FIGURE: SCENARIO COMPARISON
# ============================================================

fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(scores.columns))
width = 0.35

lfp_values = scores.loc["LFP"].to_numpy()
nmc_values = scores.loc["NMC811"].to_numpy()

bars_lfp = ax.bar(
    x - width / 2,
    lfp_values,
    width,
    label="LFP"
)

bars_nmc = ax.bar(
    x + width / 2,
    nmc_values,
    width,
    label="NMC811"
)

ax.set_xticks(x)
ax.set_xticklabels(scores.columns)

ax.set_ylabel("Weighted normalized score [-]")
ax.set_ylim(0, 1.1)

ax.set_title(
    "Technology Decision Matrix: LFP vs NMC811\n"
    "Illustrative engineering weights"
)

ax.legend()
ax.grid(axis="y", alpha=0.25)
ax.set_axisbelow(True)

for bars in [bars_lfp, bars_nmc]:

    for bar in bars:

        height = bar.get_height()

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + 0.015,
            f"{height:.3f}",
            ha="center",
            va="bottom"
        )

fig.tight_layout()

fig.savefig(
    RESULTS / "technology_decision_matrix.png",
    dpi=300
)

plt.close(fig)


# ============================================================
# 9. FIGURE: SENSITIVITY ANALYSIS
# ============================================================

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(
    sensitivity["Alpha"],
    sensitivity["LFP score"],
    label="LFP",
    linewidth=2
)

ax.plot(
    sensitivity["Alpha"],
    sensitivity["NMC811 score"],
    label="NMC811",
    linewidth=2
)

if switching_alpha is not None:

    ax.axvline(
        switching_alpha,
        linestyle="--",
        color="gray",
        label=f"Switching point: {switching_alpha:.3f}"
    )

ax.set_xlabel(
    "Weight interpolation parameter alpha\n"
    "0 = Safety-oriented, 1 = Energy-density-oriented"
)

ax.set_ylabel("Weighted normalized score [-]")

ax.set_title(
    "Decision Sensitivity to Engineering Priorities"
)

ax.grid(alpha=0.25)
ax.legend()

fig.tight_layout()

fig.savefig(
    RESULTS / "technology_sensitivity_analysis.png",
    dpi=300
)

plt.close(fig)


# ============================================================
# 10. GENERATE ENGINEERING REPORT
# ============================================================

safety_winner = scores[
    "Safety-oriented"
].idxmax()

energy_winner = scores[
    "Energy-density-oriented"
].idxmax()

report = f"""# Technology Decision Matrix: LFP vs NMC811

## Objective

Compare LFP and NMC811 battery cells using a
multi-criteria engineering decision framework.

## Data sources

Thermal-runaway indicators:
Derived from the experimental analysis in Steps 21-23.

Energy-density indicators:
Reported characteristics of the tested cells in
Schoberl et al., eTransportation 19 (2024), 100305.

## Method

Weighted-sum multi-criteria decision analysis.

Benefit criteria:
score = value / maximum value

Cost criteria:
score = minimum value / value

The normalized scores are relative to the two
technologies compared here.

## Scenario results

Safety-oriented:
LFP = {scores.loc["LFP", "Safety-oriented"]:.4f}
NMC811 = {scores.loc["NMC811", "Safety-oriented"]:.4f}
Preferred technology = {safety_winner}

Energy-density-oriented:
LFP = {scores.loc["LFP", "Energy-density-oriented"]:.4f}
NMC811 = {scores.loc["NMC811", "Energy-density-oriented"]:.4f}
Preferred technology = {energy_winner}

## Sensitivity analysis

"""

if switching_alpha is not None:

    report += (
        f"A ranking reversal occurs at "
        f"alpha = {switching_alpha:.4f}.\n"
    )

else:

    report += (
        "No ranking reversal occurs along the "
        "selected interpolation path.\n"
    )

report += """
## Important limitations

1. TR onset estimates remain preliminary.

2. The experiments used different cell geometries
   and different heater powers.

3. Thermal-runaway propagation time is not equivalent
   to guaranteed detection or response time.

4. Temperature measurements have sensor-specific
   limitations.

5. The five criteria are not statistically independent.
   In particular, the three thermal-safety indicators
   may partially describe the same underlying behavior.

6. Energy densities refer to the tested cells,
   not complete automotive battery packs.

7. Criterion weights are illustrative engineering
   assumptions, not experimentally measured quantities.

8. Weighted scores are relative decision-support
   indicators, not absolute safety probabilities.

9. No uncertainty intervals or repeated-experiment
   variability are included in the weighted scores.

10. A real vehicle technology decision would require
    additional information about cost, cycle life,
    charging performance, operating temperatures,
    pack architecture and system-level safety.

## Engineering interpretation

The preferred battery chemistry depends on the
application requirements and chosen priorities.

This decision matrix demonstrates a transparent
trade-off between measured thermal-runaway behavior
and cell-level energy density.

It does not establish universal superiority
of either chemistry.
"""

report_path = RESULTS / "technology_decision_report.md"

report_path.write_text(
    report,
    encoding="utf-8"
)


# ============================================================
# 11. FINAL OUTPUT
# ============================================================

print("\n" + "=" * 75)
print("SAVED FILES")
print("=" * 75)

saved_files = [
    "technology_raw_criteria.csv",
    "technology_normalized_criteria.csv",
    "technology_scenario_weights.csv",
    "technology_decision_scores.csv",
    "technology_weighted_contributions.csv",
    "technology_sensitivity_analysis.csv",
    "technology_decision_matrix.png",
    "technology_sensitivity_analysis.png",
    "technology_decision_report.md",
]

for filename in saved_files:
    print(filename)

print("\nStep 24 completed successfully.")