from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# STEP 22 - INTEGRATED EXPERIMENTAL ANALYSIS
#
# Battery Thermal Safety Assessment:
# LFP vs NMC811
#
# Inputs:
#   - TR onset validation (Step 21)
#   - Temperature severity (Step 13)
#   - Venting severity
#   - Surface thermal front propagation (Step 20.3)
#
# Outputs:
#   - integrated_cell_summary.csv
#   - integrated_chemistry_summary.csv
#   - integrated_experimental_comparison.png
#
# IMPORTANT:
# TR onset estimates remain preliminary.
# Surface threshold crossing times are NOT TR onset times.
# IR measurements above 250 C have increased uncertainty.
# ============================================================


# ============================================================
# 1. CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

RESULTS = ROOT / "results"

RESULTS.mkdir(parents=True, exist_ok=True)


TR_FILE = RESULTS / "tr_onset_validation_summary.csv"

TEMP_FILE = RESULTS / "temperature_severity_summary.csv"

VENT_FILE = RESULTS / "venting_severity_summary.csv"

FRONT_FILE = RESULTS / "surface_front_summary.csv"


OUTPUT_CELL = RESULTS / "integrated_cell_summary.csv"

OUTPUT_CHEM = RESULTS / "integrated_chemistry_summary.csv"

OUTPUT_FIG = RESULTS / "integrated_experimental_comparison.png"


CHEMISTRIES = ["LFP", "NMC811"]

CELLS = range(1, 6)

SURFACE_THRESHOLDS = [100, 150, 200]


# Cell-associated thermocouples used for temperature comparison.
#
# T6 is an additional measurement position and is not included
# in the cell-associated temperature mean.
#
# NMC811 T5 is excluded because of known measurement artifacts.

VALID_TEMPERATURE_SENSORS = {
    "LFP": {1, 2, 3, 4, 5},
    "NMC811": {1, 2, 3, 4},
}


# Known measurement limitation:
# LFP T1 was damaged after approximately 185.9 s.

LFP_T1_DAMAGE_TIME = 185.9


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

def load_csv(path):

    if not path.exists():

        raise FileNotFoundError(
            f"\nRequired file not found:\n{path}\n"
            "Check whether the corresponding analysis "
            "step has been completed."
        )

    df = pd.read_csv(path)

    df.columns = df.columns.str.strip()

    print(f"\nLoaded: {path.name}")

    print(f"  Rows: {len(df)}")

    print(f"  Columns: {list(df.columns)}")

    return df


def require_columns(df, required_columns, filename):

    missing = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing:

        raise ValueError(
            f"\nMissing columns in {filename}:\n{missing}"
        )


def numeric_values(series):

    return pd.to_numeric(
        series,
        errors="coerce"
    )


def safe_mean(series):

    values = numeric_values(series).dropna()

    if values.empty:
        return np.nan

    return float(values.mean())


def safe_max(series):

    values = numeric_values(series).dropna()

    if values.empty:
        return np.nan

    return float(values.max())


def get_single_row(df, chemistry, cell, filename):

    selected = df[
        (df["Chemistry"] == chemistry)
        & (df["Cell"] == cell)
    ]

    if len(selected) != 1:

        raise ValueError(
            f"\nExpected exactly one row for "
            f"{chemistry}, Cell {cell} "
            f"in {filename}.\n"
            f"Found: {len(selected)}"
        )

    return selected.iloc[0]


# ============================================================
# 3. LOAD INPUT FILES
# ============================================================

print("\n" + "=" * 75)

print("STEP 22 - INTEGRATED EXPERIMENTAL ANALYSIS")

print("=" * 75)


tr = load_csv(TR_FILE)

temp = load_csv(TEMP_FILE)

vent = load_csv(VENT_FILE)

front = load_csv(FRONT_FILE)


# ============================================================
# 4. CHECK REQUIRED COLUMNS
# ============================================================

require_columns(
    tr,
    [
        "Chemistry",
        "Cell",
        "Candidate TR time [s]"
    ],
    TR_FILE.name
)


require_columns(
    temp,
    [
        "Chemistry",
        "Channel",
        "Tmax [C]",
        "Time at Tmax [s]"
    ],
    TEMP_FILE.name
)


require_columns(
    vent,
    [
        "Chemistry",
        "Cell",
        "Tmax [C]"
    ],
    VENT_FILE.name
)


require_columns(
    front,
    [
        "Chemistry",
        "Threshold [C]",
        "Crossing time span [s]"
    ],
    FRONT_FILE.name
)


# ============================================================
# 5. STANDARDIZE DATA
# ============================================================

for df in [tr, temp, vent, front]:

    df["Chemistry"] = (
        df["Chemistry"]
        .astype(str)
        .str.strip()
    )


tr["Cell"] = numeric_values(tr["Cell"])

vent["Cell"] = numeric_values(vent["Cell"])


tr["Candidate TR time [s]"] = numeric_values(
    tr["Candidate TR time [s]"]
)


temp["Tmax [C]"] = numeric_values(
    temp["Tmax [C]"]
)


temp["Time at Tmax [s]"] = numeric_values(
    temp["Time at Tmax [s]"]
)


vent["Tmax [C]"] = numeric_values(
    vent["Tmax [C]"]
)


front["Threshold [C]"] = numeric_values(
    front["Threshold [C]"]
)


front["Crossing time span [s]"] = numeric_values(
    front["Crossing time span [s]"]
)


# Extract sensor numbers from T1, T2, etc.

temp["Sensor number"] = (
    temp["Channel"]
    .astype(str)
    .str.extract(r"^T(\d+)$")[0]
)


temp["Sensor number"] = numeric_values(
    temp["Sensor number"]
)


# ============================================================
# 6. VALIDATE INPUT DATA
# ============================================================

print("\n" + "=" * 75)

print("INPUT DATA VALIDATION")

print("=" * 75)


expected_cells = {
    (chemistry, cell)
    for chemistry in CHEMISTRIES
    for cell in CELLS
}


actual_tr_cells = {
    (row["Chemistry"], int(row["Cell"]))
    for _, row in tr.dropna(subset=["Cell"]).iterrows()
}


if actual_tr_cells != expected_cells:

    raise ValueError(
        "\nTR onset data does not contain "
        "the expected cells.\n"
        f"Missing: {expected_cells - actual_tr_cells}\n"
        f"Unexpected: {actual_tr_cells - expected_cells}"
    )


if tr.duplicated(
    subset=["Chemistry", "Cell"]
).any():

    raise ValueError(
        "Duplicate TR onset entries detected."
    )


if vent.duplicated(
    subset=["Chemistry", "Cell"]
).any():

    raise ValueError(
        "Duplicate venting entries detected."
    )


if tr["Candidate TR time [s]"].isna().any():

    raise ValueError(
        "Missing TR onset candidate times."
    )


if not np.isfinite(
    tr["Candidate TR time [s]"]
).all():

    raise ValueError(
        "Non-finite TR onset candidate times."
    )


actual_vent_cells = {
    (row["Chemistry"], int(row["Cell"]))
    for _, row in vent.dropna(subset=["Cell"]).iterrows()
}


if actual_vent_cells != expected_cells:

    raise ValueError(
        "\nVenting data does not contain "
        "the expected cells."
    )


if temp.duplicated(
    subset=["Chemistry", "Sensor number"]
).any():

    raise ValueError(
        "Duplicate temperature sensor entries detected."
    )


# ============================================================
# 7. VALIDATE TEMPERATURE SENSOR SELECTION
# ============================================================

print("\nTemperature sensor validation:")


for chemistry in CHEMISTRIES:

    chemistry_temp = temp[
        temp["Chemistry"] == chemistry
    ].copy()


    selected_temp = chemistry_temp[
        chemistry_temp["Sensor number"].isin(
            VALID_TEMPERATURE_SENSORS[chemistry]
        )
    ]


    actual_sensors = set(
        selected_temp["Sensor number"]
        .dropna()
        .astype(int)
    )


    expected_sensors = VALID_TEMPERATURE_SENSORS[
        chemistry
    ]


    if actual_sensors != expected_sensors:

        raise ValueError(
            f"\nUnexpected temperature sensors "
            f"for {chemistry}.\n"
            f"Expected: {expected_sensors}\n"
            f"Found: {actual_sensors}"
        )


    if selected_temp["Tmax [C]"].isna().any():

        raise ValueError(
            f"Missing Tmax values for {chemistry}."
        )


    if not np.isfinite(
        selected_temp["Tmax [C]"]
    ).all():

        raise ValueError(
            f"Non-finite Tmax values for {chemistry}."
        )


    print(
        f"{chemistry}: "
        f"{sorted(actual_sensors)}"
    )


# ============================================================
# 8. VALIDATE LFP T1 MEASUREMENT
# ============================================================

lfp_t1 = temp[
    (temp["Chemistry"] == "LFP")
    & (temp["Sensor number"] == 1)
]


if len(lfp_t1) != 1:

    raise ValueError(
        "Expected exactly one LFP T1 entry."
    )


lfp_t1_max = float(
    lfp_t1["Tmax [C]"].iloc[0]
)


lfp_t1_max_time = float(
    lfp_t1["Time at Tmax [s]"].iloc[0]
)


print("\nLFP T1 measurement check:")

print(
    f"Maximum temperature: "
    f"{lfp_t1_max:.2f} C"
)

print(
    f"Time at maximum: "
    f"{lfp_t1_max_time:.3f} s"
)

print(
    f"Known sensor damage time: "
    f"{LFP_T1_DAMAGE_TIME:.3f} s"
)


if not np.isfinite(lfp_t1_max_time):

    raise ValueError(
        "Invalid time at maximum for LFP T1."
    )


if lfp_t1_max_time > LFP_T1_DAMAGE_TIME:

    raise ValueError(
        "\nLFP T1 maximum occurs after the "
        "known sensor damage time.\n"
        "Temperature severity results must "
        "be reviewed before continuing."
    )


print(
    "LFP T1 maximum occurs before the "
    "reported damage time."
)

print(
    "Note: this does not independently prove "
    "measurement accuracy near sensor failure."
)


print("\nInput validation passed.")


# ============================================================
# 9. CELL-LEVEL ANALYSIS
# ============================================================

print("\n" + "=" * 75)

print("CELL-LEVEL ANALYSIS")

print("=" * 75)


cell_rows = []


for chemistry in CHEMISTRIES:

    chemistry_tr = (
        tr[tr["Chemistry"] == chemistry]
        .sort_values("Cell")
        .copy()
    )


    previous_tr_time = None


    for cell in CELLS:

        tr_row = get_single_row(
            chemistry_tr,
            chemistry,
            cell,
            TR_FILE.name
        )


        tr_time = float(
            tr_row["Candidate TR time [s]"]
        )


        # --------------------------------------------
        # TR propagation delay
        # --------------------------------------------

        if previous_tr_time is None:

            propagation_delay = np.nan

        else:

            propagation_delay = (
                tr_time - previous_tr_time
            )


            if propagation_delay <= 0:

                raise ValueError(
                    f"Invalid TR propagation sequence: "
                    f"{chemistry}, Cell {cell}"
                )


        previous_tr_time = tr_time


        # --------------------------------------------
        # Cell-associated temperature sensor
        # --------------------------------------------

        temperature_sensor = f"T{cell}"


        if (
            cell
            not in VALID_TEMPERATURE_SENSORS[chemistry]
        ):

            cell_temp_max = np.nan

            temp_quality = (
                "EXCLUDED - known sensor artifacts"
            )

        else:

            selected_temp = temp[
                (temp["Chemistry"] == chemistry)
                & (temp["Sensor number"] == cell)
            ]


            if len(selected_temp) != 1:

                raise ValueError(
                    f"Missing temperature sensor "
                    f"{temperature_sensor} "
                    f"for {chemistry}."
                )


            cell_temp_max = float(
                selected_temp["Tmax [C]"].iloc[0]
            )


            if chemistry == "LFP" and cell == 1:

                temp_quality = (
                    "AVAILABLE - near sensor damage time"
                )

            else:

                temp_quality = "AVAILABLE"


        # --------------------------------------------
        # Venting temperature
        # --------------------------------------------

        vent_row = get_single_row(
            vent,
            chemistry,
            cell,
            VENT_FILE.name
        )


        vent_max = float(
            vent_row["Tmax [C]"]
        )


        # --------------------------------------------
        # Store results
        # --------------------------------------------

        cell_rows.append({

            "Chemistry": chemistry,

            "Cell": cell,

            "TR onset candidate [s]": tr_time,

            "Propagation delay [s]": propagation_delay,

            "Temperature sensor": temperature_sensor,

            "Cell-associated Tmax [C]": cell_temp_max,

            "Temperature data quality": temp_quality,

            "Venting sensor": f"V{cell}",

            "Venting Tmax [C]": vent_max,

            "TR onset status": (
                "PRELIMINARY - manual validation pending"
            )

        })


cell_summary = pd.DataFrame(cell_rows)


# ============================================================
# 10. CHEMISTRY-LEVEL ANALYSIS
# ============================================================

print("\n" + "=" * 75)

print("CHEMISTRY-LEVEL ANALYSIS")

print("=" * 75)


chem_rows = []


for chemistry in CHEMISTRIES:

    subset = (
        cell_summary[
            cell_summary["Chemistry"] == chemistry
        ]
        .sort_values("Cell")
        .copy()
    )


    # --------------------------------------------
    # TR propagation
    # --------------------------------------------

    tr_times = subset[
        "TR onset candidate [s]"
    ].to_numpy(dtype=float)


    propagation_delays = np.diff(tr_times)


    total_propagation = (
        tr_times[-1] - tr_times[0]
    )


    mean_propagation_delay = float(
        np.mean(propagation_delays)
    )


    sd_propagation_delay = float(
        np.std(
            propagation_delays,
            ddof=1
        )
    )


    # --------------------------------------------
    # Temperature severity
    # --------------------------------------------

    temp_subset = temp[
        (temp["Chemistry"] == chemistry)
        & (
            temp["Sensor number"].isin(
                VALID_TEMPERATURE_SENSORS[chemistry]
            )
        )
    ].copy()


    temp_subset = temp_subset.sort_values(
        "Sensor number"
    )


    number_of_temp_sensors = len(temp_subset)


    mean_temp_max = safe_mean(
        temp_subset["Tmax [C]"]
    )


    maximum_temp = safe_max(
        temp_subset["Tmax [C]"]
    )


    # --------------------------------------------
    # Venting severity
    # --------------------------------------------

    vent_subset = vent[
        vent["Chemistry"] == chemistry
    ].copy()


    mean_vent_max = safe_mean(
        vent_subset["Tmax [C]"]
    )


    maximum_vent = safe_max(
        vent_subset["Tmax [C]"]
    )


    # --------------------------------------------
    # Chemistry summary
    # --------------------------------------------

    row = {

        "Chemistry": chemistry,

        "First TR onset candidate [s]": tr_times[0],

        "Last TR onset candidate [s]": tr_times[-1],

        "Total TR propagation [s]": total_propagation,

        "Mean inter-cell TR delay [s]":
            mean_propagation_delay,

        "SD inter-cell TR delay [s]":
            sd_propagation_delay,

        "Temperature sensors used":
            number_of_temp_sensors,

        "Mean cell-associated Tmax [C]":
            mean_temp_max,

        "Maximum cell-associated Tmax [C]":
            maximum_temp,

        "Mean venting Tmax [C]":
            mean_vent_max,

        "Maximum venting Tmax [C]":
            maximum_vent,

    }


    # --------------------------------------------
    # Surface threshold crossing spans
    # --------------------------------------------

    for threshold in SURFACE_THRESHOLDS:

        selected_front = front[
            (front["Chemistry"] == chemistry)
            & np.isclose(
                front["Threshold [C]"],
                threshold
            )
        ]


        if len(selected_front) != 1:

            raise ValueError(
                f"Expected one surface-front entry "
                f"for {chemistry}, {threshold} C.\n"
                f"Found: {len(selected_front)}"
            )


        crossing_span = float(
            selected_front[
                "Crossing time span [s]"
            ].iloc[0]
        )


        if not np.isfinite(crossing_span):

            raise ValueError(
                f"Invalid surface crossing span "
                f"for {chemistry}, {threshold} C."
            )


        row[
            f"Surface {threshold}C crossing span [s]"
        ] = crossing_span


    chem_rows.append(row)


chem_summary = pd.DataFrame(chem_rows)


# ============================================================
# 11. SAVE CSV RESULTS
# ============================================================

cell_summary.to_csv(
    OUTPUT_CELL,
    index=False,
    float_format="%.3f"
)


chem_summary.to_csv(
    OUTPUT_CHEM,
    index=False,
    float_format="%.3f"
)


# ============================================================
# 12. COMPARATIVE VISUALIZATION
# ============================================================

print("\n" + "=" * 75)

print("GENERATING COMPARATIVE FIGURE")

print("=" * 75)


plot_metrics = [

    (
        "Total TR propagation [s]",
        "Total TR propagation",
        "Time [s]"
    ),

    (
        "Mean inter-cell TR delay [s]",
        "Mean inter-cell TR delay",
        "Time [s]"
    ),

    (
        "Mean cell-associated Tmax [C]",
        "Mean cell-associated Tmax\n"
        "(LFP n=5; NMC811 n=4)",
        "Temperature [°C]"
    ),

    (
        "Mean venting Tmax [C]",
        "Mean venting-sensor Tmax",
        "Temperature [°C]"
    ),

    (
        "Surface 150C crossing span [s]",
        "150°C surface crossing span",
        "Time [s]"
    )

]


fig, axes = plt.subplots(
    2,
    3,
    figsize=(17, 10)
)


axes = axes.flatten()


colors = {
    "LFP": "#377eb8",
    "NMC811": "#e6550d"
}


for ax, (column, title, ylabel) in zip(
    axes,
    plot_metrics
):

    values = []

    for chemistry in CHEMISTRIES:

        selected = chem_summary[
            chem_summary["Chemistry"] == chemistry
        ]


        value = float(
            selected[column].iloc[0]
        )


        values.append(value)


    bars = ax.bar(
        CHEMISTRIES,
        values,
        color=[
            colors[chemistry]
            for chemistry in CHEMISTRIES
        ],
        width=0.55
    )


    ax.set_title(
        title,
        fontsize=12
    )


    ax.set_ylabel(
        ylabel,
        fontsize=11
    )


    ax.grid(
        axis="y",
        alpha=0.3
    )


    ax.set_axisbelow(True)


    max_value = max(values)


    ax.set_ylim(
        0,
        max_value * 1.18
    )


    for bar, value in zip(bars, values):

        if np.isfinite(value):

            ax.annotate(
                f"{value:.1f}",
                (
                    bar.get_x()
                    + bar.get_width() / 2,
                    bar.get_height()
                ),
                ha="center",
                va="bottom",
                xytext=(0, 5),
                textcoords="offset points",
                fontsize=10
            )


# Hide unused subplot.

axes[-1].axis("off")


fig.suptitle(
    "LFP vs NMC811 - Integrated Experimental Comparison\n"
    "Specific experimental conditions; "
    "preliminary TR onset estimates\n"
    "Sensor limitations and IR measurement uncertainty apply",
    fontsize=15
)


fig.tight_layout(
    rect=[0, 0, 1, 0.90]
)


fig.savefig(
    OUTPUT_FIG,
    dpi=300,
    bbox_inches="tight"
)


plt.close(fig)


# ============================================================
# 13. PRINT CELL-LEVEL RESULTS
# ============================================================

print("\n" + "=" * 75)

print("CELL-LEVEL RESULTS")

print("=" * 75)


print(
    cell_summary.to_string(
        index=False
    )
)


# ============================================================
# 14. PRINT CHEMISTRY COMPARISON
# ============================================================

print("\n" + "=" * 75)

print("CHEMISTRY COMPARISON")

print("=" * 75)


print(
    chem_summary.round(2).to_string(
        index=False
    )
)


# ============================================================
# 15. PROPAGATION RATIO
# ============================================================

print("\n" + "=" * 75)

print("PROPAGATION RATIO")

print("=" * 75)


lfp_total = float(
    chem_summary.loc[
        chem_summary["Chemistry"] == "LFP",
        "Total TR propagation [s]"
    ].iloc[0]
)


nmc_total = float(
    chem_summary.loc[
        chem_summary["Chemistry"] == "NMC811",
        "Total TR propagation [s]"
    ].iloc[0]
)


if nmc_total <= 0:

    raise ValueError(
        "NMC811 propagation duration must be positive."
    )


propagation_ratio = (
    lfp_total / nmc_total
)


print(
    f"LFP / NMC811 propagation time ratio: "
    f"{propagation_ratio:.2f}"
)


print(
    "\nInterpretation:"
)

print(
    "This is a ratio of total TR propagation durations."
)

print(
    "It is NOT a ratio of physical thermal-front "
    "velocities in m/s."
)


# ============================================================
# 16. MEASUREMENT LIMITATIONS
# ============================================================

print("\n" + "=" * 75)

print("MEASUREMENT LIMITATIONS")

print("=" * 75)


print(
    "\n1. TR onset times are preliminary candidates."
)


print(
    "2. LFP temperature mean uses T1-T5 (n=5)."
)


print(
    "3. NMC811 temperature mean uses T1-T4 (n=4)."
)


print(
    "4. NMC811 T5 is excluded due to sensor artifacts."
)


print(
    "5. T6 is excluded from cell-associated "
    "temperature means for both chemistries."
)


print(
    "6. LFP T1 Tmax occurs shortly before "
    "the reported sensor damage time."
)


print(
    "7. Venting temperatures are sensor measurements, "
    "not direct measurements of total heat release."
)


print(
    "8. IR measurements above 250 C have "
    "increased measurement uncertainty."
)


print(
    "9. Surface threshold crossing spans "
    "are not TR propagation durations."
)


print(
    "10. The experiments used different cell "
    "geometries and heating powers."
)


# ============================================================
# 17. FINAL OUTPUT
# ============================================================

print("\n" + "=" * 75)

print("SAVED FILES")

print("=" * 75)


print(OUTPUT_CELL.name)

print(OUTPUT_CHEM.name)

print(OUTPUT_FIG.name)


print("\nStep 22 completed successfully.")