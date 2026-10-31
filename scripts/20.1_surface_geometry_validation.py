from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# STEP 20.1 - SURFACE GEOMETRY & DATA VALIDATION
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT / "data" / "surface"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]

# Thermal camera measurement limitation:
# Emissivity tape specified up to approximately 250 C.
TAPE_LIMIT_C = 250.0

# Sampling gap diagnostic
MAX_GAP_FACTOR = 3.0

# Channel number pattern
CHANNEL_PATTERN = re.compile(
    r"^Temperature\s+(\d+)\s+\[°C\]$"
)

all_channel_rows = []
all_summary_rows = []


# ============================================================
# LOAD AND VALIDATE
# ============================================================

for chemistry in CHEMISTRIES:

    file = (
        DATA_DIR
        / chemistry
        / f"T_surf_{chemistry}.txt"
    )

    if not file.exists():
        raise FileNotFoundError(file)

    df = pd.read_csv(
        file,
        encoding="utf-8-sig"
    )

    if "Time [s]" not in df.columns:
        raise ValueError(
            f"Time column missing: {file.name}"
        )

    channel_info = []

    for column in df.columns:

        match = CHANNEL_PATTERN.match(str(column))

        if match:
            channel_info.append(
                (int(match.group(1)), column)
            )

    channel_info.sort(key=lambda x: x[0])

    channel_numbers = np.array(
        [item[0] for item in channel_info],
        dtype=int
    )

    channel_columns = [
        item[1] for item in channel_info
    ]

    N = len(channel_numbers)

    if N < 2:
        raise ValueError(
            f"Insufficient channels: {file.name}"
        )

    expected_numbers = np.arange(1, N + 1)

    numbering_ok = np.array_equal(
        channel_numbers,
        expected_numbers
    )

    if not numbering_ok:
        raise ValueError(
            f"Unexpected channel numbering in {file.name}"
        )

    # --------------------------------------------------------
    # TIME VALIDATION
    # --------------------------------------------------------

    time = pd.to_numeric(
        df["Time [s]"],
        errors="coerce"
    ).to_numpy(dtype=float)

    if not np.isfinite(time).all():
        raise ValueError(
            f"Invalid time values: {file.name}"
        )

    dt = np.diff(time)

    if np.any(dt <= 0):
        raise ValueError(
            f"Non-increasing timestamps: {file.name}"
        )

    median_dt = float(np.median(dt))

    large_gaps = int(
        np.sum(dt > MAX_GAP_FACTOR * median_dt)
    )

    # --------------------------------------------------------
    # TEMPERATURE MATRIX
    # --------------------------------------------------------

    temperature_df = df[channel_columns].apply(
        pd.to_numeric,
        errors="coerce"
    )

    matrix = temperature_df.to_numpy(
        dtype=float
    )

    finite = np.isfinite(matrix)

    invalid_samples = int(
        np.sum(~finite)
    )

    negative_samples = int(
        np.sum(matrix[finite] < 0)
    )

    # Use t >= 0 for thermal exposure diagnostics.
    active_mask = time >= 0

    active_matrix = matrix[active_mask, :]

    finite_active = np.isfinite(active_matrix)

    above_limit = (
        finite_active
        & (active_matrix > TAPE_LIMIT_C)
    )

    # --------------------------------------------------------
    # NORMALIZED CHANNEL COORDINATE
    # --------------------------------------------------------

    x_normalized = (
        channel_numbers - 1
    ) / (N - 1)

    # --------------------------------------------------------
    # CHANNEL-BY-CHANNEL DIAGNOSTICS
    # --------------------------------------------------------

    for j, (number, column) in enumerate(
        channel_info
    ):

        T = matrix[:, j]

        valid = np.isfinite(T)

        T_active = active_matrix[:, j]

        valid_active = np.isfinite(T_active)

        if valid.any():

            Tmin = float(np.min(T[valid]))
            Tmax = float(np.max(T[valid]))

        else:

            Tmin = np.nan
            Tmax = np.nan

        active_count = int(
            valid_active.sum()
        )

        above_count = int(
            np.sum(
                valid_active
                & (T_active > TAPE_LIMIT_C)
            )
        )

        fraction_above = (
            above_count / active_count
            if active_count > 0
            else np.nan
        )

        row = {
            "Chemistry": chemistry,
            "Channel": column,
            "Channel number": number,
            "Normalized channel coordinate": float(
                x_normalized[j]
            ),
            "Physical x [mm]": np.nan,
            "Cell assignment": "Unverified",
            "Tmin [C]": Tmin,
            "Tmax [C]": Tmax,
            "Valid active samples": active_count,
            "Samples above 250C": above_count,
            "Fraction above 250C": fraction_above,
            "Absolute temperature reliability":
                (
                    "Exceeds tape rating"
                    if above_count > 0
                    else "Within tape rating"
                ),
        }

        all_channel_rows.append(row)

    # --------------------------------------------------------
    # CHEMISTRY SUMMARY
    # --------------------------------------------------------

    total_active_valid = int(
        finite_active.sum()
    )

    total_above_limit = int(
        above_limit.sum()
    )

    channels_above_limit = int(
        np.sum(
            np.any(above_limit, axis=0)
        )
    )

    summary = {
        "Chemistry": chemistry,
        "Channels": N,
        "Samples per channel": len(df),
        "Start time [s]": float(time[0]),
        "End time [s]": float(time[-1]),
        "Median dt [s]": median_dt,
        "Large time gaps": large_gaps,
        "Invalid samples": invalid_samples,
        "Negative samples": negative_samples,
        "Channels exceeding 250C":
            channels_above_limit,
        "Total samples above 250C":
            total_above_limit,
        "Fraction of active samples above 250C":
            (
                total_above_limit / total_active_valid
                if total_active_valid > 0
                else np.nan
            ),
        "Channel numbering verified": numbering_ok,
        "Physical coordinates verified": False,
    }

    all_summary_rows.append(summary)

    # --------------------------------------------------------
    # TERMINAL OUTPUT
    # --------------------------------------------------------

    print("\n" + "=" * 75)
    print(f"SURFACE GEOMETRY VALIDATION - {chemistry}")
    print("=" * 75)

    print(f"File: {file.name}")
    print(f"Channels: {N}")
    print(f"Samples: {len(df):,}")
    print(f"Time: {time[0]:.2f} to {time[-1]:.2f} s")
    print(f"Median dt: {median_dt:.5f} s")

    print("\nDATA QUALITY:")

    print(f"Invalid samples: {invalid_samples}")
    print(f"Negative samples: {negative_samples}")
    print(f"Large time gaps: {large_gaps}")

    print("\nGEOMETRY:")

    print(
        f"Sequential channel numbering: {numbering_ok}"
    )

    print(
        "Normalized coordinate range: "
        f"{x_normalized[0]:.2f} to "
        f"{x_normalized[-1]:.2f}"
    )

    print(
        "Physical coordinates: NOT VERIFIED"
    )

    print(
        "Cell-to-channel assignment: NOT VERIFIED"
    )

    print("\nTHERMAL CAMERA LIMITATION:")

    print(
        f"Channels exceeding {TAPE_LIMIT_C:.0f} C: "
        f"{channels_above_limit}/{N}"
    )

    print(
        f"Samples above {TAPE_LIMIT_C:.0f} C: "
        f"{total_above_limit:,}"
    )

    print(
        "Fraction of active samples above limit: "
        f"{summary['Fraction of active samples above 250C']:.2%}"
    )

    # --------------------------------------------------------
    # FIGURE 1 - NORMALIZED CHANNEL GEOMETRY
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(13, 4)
    )

    ax.scatter(
        x_normalized,
        np.zeros(N),
        s=25
    )

    ax.set_title(
        f"{chemistry} - Normalized Surface Channel Positions"
    )

    ax.set_xlabel(
        "Normalized channel index coordinate [-]"
    )

    ax.set_yticks([])

    ax.set_xlim(-0.03, 1.03)

    ax.grid(
        axis="x",
        alpha=0.3
    )

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_normalized_geometry_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    # --------------------------------------------------------
    # FIGURE 2 - TEMPERATURE LIMIT BY CHANNEL
    # --------------------------------------------------------

    channel_above_fraction = (
        np.sum(above_limit, axis=0)
        / np.maximum(
            np.sum(finite_active, axis=0),
            1
        )
    )

    fig, ax = plt.subplots(
        figsize=(13, 5)
    )

    ax.plot(
        channel_numbers,
        channel_above_fraction * 100,
        marker=".",
        linewidth=1.2
    )

    ax.set_title(
        f"{chemistry} - Fraction of Samples Above 250 C"
    )

    ax.set_xlabel("Temperature channel number")

    ax.set_ylabel(
        "Fraction of valid active samples [%]"
    )

    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_above_250C_fraction_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)


# ============================================================
# SAVE CSV FILES
# ============================================================

channel_summary = pd.DataFrame(
    all_channel_rows
)

chemistry_summary = pd.DataFrame(
    all_summary_rows
)

channel_summary.to_csv(
    RESULTS_DIR
    / "surface_geometry_channel_mapping.csv",
    index=False
)

chemistry_summary.to_csv(
    RESULTS_DIR
    / "surface_geometry_validation_summary.csv",
    index=False
)

print("\n" + "=" * 75)
print("FINAL VALIDATION SUMMARY")
print("=" * 75)

print(
    chemistry_summary.round(4).to_string(
        index=False
    )
)

print("\nSAVED FILES:")

print("surface_geometry_channel_mapping.csv")
print("surface_geometry_validation_summary.csv")

for chemistry in CHEMISTRIES:

    print(
        f"surface_normalized_geometry_{chemistry}.png"
    )

    print(
        f"surface_above_250C_fraction_{chemistry}.png"
    )

print("\nStep 20.1 completed.")