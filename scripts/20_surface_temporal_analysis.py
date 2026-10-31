
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# STEP 20 - SURFACE TEMPERATURE TEMPORAL ANALYSIS
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "surface"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]

THRESHOLDS = [100.0, 200.0, 300.0]

# Time-based smoothing window.
SMOOTHING_WINDOW_S = 1.0

# Ignore isolated numerical derivatives.
# Heating rates are calculated from smoothed data.
MAX_GAP_FACTOR = 3.0

# Plot resolution only.
MAX_PLOT_POINTS = 10000

TR_FILE = (
    RESULTS_DIR /
    "refined_tr_detection_summary.csv"
)

# ============================================================
# HELPERS
# ============================================================

def first_crossing_time(time, temperature, threshold):
    """
    First upward crossing of a threshold at t >= 0.
    Linear interpolation between adjacent samples.
    """

    valid = (
        np.isfinite(time)
        & np.isfinite(temperature)
        & (time >= 0)
    )

    t = time[valid]
    T = temperature[valid]

    if len(t) == 0:
        return np.nan

    if T[0] >= threshold:
        return float(t[0])

    crossings = np.flatnonzero(
        (T[:-1] < threshold)
        & (T[1:] >= threshold)
    )

    if len(crossings) == 0:
        return np.nan

    i = int(crossings[0])

    dT = T[i + 1] - T[i]

    if dT <= 0:
        return float(t[i + 1])

    fraction = (
        (threshold - T[i]) / dT
    )

    return float(
        t[i] + fraction * (t[i + 1] - t[i])
    )


def calculate_smoothed_rate(
    time,
    temperature,
    window_samples,
    max_gap
):
    """
    Rolling-median smoothing followed by numerical derivative.
    """

    smoothed = (
        pd.Series(temperature)
        .rolling(
            window=window_samples,
            center=True,
            min_periods=1
        )
        .median()
        .to_numpy(dtype=float)
    )

    rate = np.full(
        len(time),
        np.nan,
        dtype=float
    )

    dt = np.diff(time)

    valid_pairs = (
        np.isfinite(smoothed[:-1])
        & np.isfinite(smoothed[1:])
        & (dt > 0)
        & (dt <= max_gap)
    )

    differences = np.diff(smoothed)

    rate_values = np.full(
        len(dt),
        np.nan,
        dtype=float
    )

    rate_values[valid_pairs] = (
        differences[valid_pairs]
        / dt[valid_pairs]
    )

    rate[1:] = rate_values

    return smoothed, rate


# ============================================================
# LOAD PREVIOUS TR ESTIMATES
# ============================================================

tr_events = {}

if TR_FILE.exists():

    tr_df = pd.read_csv(TR_FILE)

    required = {
        "Chemistry",
        "Cell",
        "Matched temperature event [s]"
    }

    if required.issubset(tr_df.columns):

        for _, row in tr_df.iterrows():

            key = (
                str(row["Chemistry"]),
                int(row["Cell"])
            )

            tr_events[key] = float(
                row["Matched temperature event [s]"]
            )

        print(
            f"Loaded {len(tr_events)} "
            "previous TR event estimates."
        )

    else:
        print(
            "WARNING: TR summary has unexpected columns."
        )

else:
    print(
        "WARNING: TR summary not found. "
        "TR markers will be skipped."
    )

# ============================================================
# MAIN ANALYSIS
# ============================================================

all_channel_rows = []
all_chemistry_rows = []

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

    time = pd.to_numeric(
        df["Time [s]"],
        errors="coerce"
    ).to_numpy(dtype=float)

    temp_columns = [
        c for c in df.columns
        if str(c).startswith("Temperature ")
    ]

    if not temp_columns:
        raise ValueError(
            f"No surface temperature channels in {file.name}"
        )

    if not np.isfinite(time).all():
        raise ValueError("Invalid timestamps.")

    dt = np.diff(time)

    if np.any(dt <= 0):
        raise ValueError("Non-increasing timestamps.")

    median_dt = float(np.median(dt))
    max_gap = MAX_GAP_FACTOR * median_dt

    window_samples = max(
        3,
        int(round(
            SMOOTHING_WINDOW_S / median_dt
        ))
    )

    if window_samples % 2 == 0:
        window_samples += 1

    temperatures = df[temp_columns].apply(
        pd.to_numeric,
        errors="coerce"
    )

    print("\n" + "=" * 75)
    print(f"SURFACE TEMPORAL ANALYSIS - {chemistry}")
    print("=" * 75)

    print(f"Channels: {len(temp_columns)}")
    print(f"Samples: {len(df):,}")
    print(
        f"Smoothing window: "
        f"{window_samples} samples "
        f"(~{SMOOTHING_WINDOW_S:.1f} s)"
    )

    chemistry_rows = []

    for col in temp_columns:

        T = temperatures[col].to_numpy(
            dtype=float
        )

        smoothed, rate = calculate_smoothed_rate(
            time,
            T,
            window_samples,
            max_gap
        )

        positive_time = time >= 0

        rate_for_max = rate.copy()
        rate_for_max[~positive_time] = np.nan

        if np.isfinite(rate_for_max).any():

            idx_rate = int(
                np.nanargmax(rate_for_max)
            )

            max_heating_rate = float(
                rate_for_max[idx_rate]
            )

            time_max_heating = float(
                time[idx_rate]
            )

        else:
            max_heating_rate = np.nan
            time_max_heating = np.nan

        row = {
            "Chemistry": chemistry,
            "Channel": col,
            "Max smoothed heating rate [C/s]":
                max_heating_rate,
            "Time of max heating [s]":
                time_max_heating,
        }

        for threshold in THRESHOLDS:

            row[
                f"First crossing {int(threshold)}C [s]"
            ] = first_crossing_time(
                time,
                smoothed,
                threshold
            )

        chemistry_rows.append(row)
        all_channel_rows.append(row)

    channel_df = pd.DataFrame(
        chemistry_rows
    )

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    summary = {
        "Chemistry": chemistry,
        "Channels": len(temp_columns),
        "Median max heating rate [C/s]": float(
            channel_df[
                "Max smoothed heating rate [C/s]"
            ].median()
        ),
        "Highest max heating rate [C/s]": float(
            channel_df[
                "Max smoothed heating rate [C/s]"
            ].max()
        ),
    }

    for threshold in THRESHOLDS:

        column = (
            f"First crossing {int(threshold)}C [s]"
        )

        summary[
            f"Channels reaching {int(threshold)}C"
        ] = int(
            channel_df[column].notna().sum()
        )

        summary[
            f"Median first crossing {int(threshold)}C [s]"
        ] = float(
            channel_df[column].median()
        ) if channel_df[column].notna().any() else np.nan

    all_chemistry_rows.append(summary)

    print("\nTHRESHOLD CROSSINGS:")

    for threshold in THRESHOLDS:

        column = (
            f"First crossing {int(threshold)}C [s]"
        )

        count = channel_df[column].notna().sum()

        print(
            f"{int(threshold)} C: "
            f"{count}/{len(temp_columns)} channels"
        )

    print("\nHEATING RATE:")

    print(
        "Median maximum heating rate: "
        f"{summary['Median max heating rate [C/s]']:.2f} C/s"
    )

    print(
        "Highest maximum heating rate: "
        f"{summary['Highest max heating rate [C/s]']:.2f} C/s"
    )

    # --------------------------------------------------------
    # FIGURE 1 - TEMPERATURE ENVELOPE
    # --------------------------------------------------------

    matrix = temperatures.to_numpy(
        dtype=float
    )

    plot_step = max(
        1,
        len(time) // MAX_PLOT_POINTS
    )

    t_plot = time[::plot_step]

    matrix_plot = matrix[::plot_step, :]

    mean_T = np.nanmean(
        matrix_plot,
        axis=1
    )

    median_T = np.nanmedian(
        matrix_plot,
        axis=1
    )

    min_T = np.nanmin(
        matrix_plot,
        axis=1
    )

    max_T = np.nanmax(
        matrix_plot,
        axis=1
    )

    fig, ax = plt.subplots(
        figsize=(14, 7)
    )

    ax.fill_between(
        t_plot,
        min_T,
        max_T,
        alpha=0.18,
        label="Channel min-max range"
    )

    ax.plot(
        t_plot,
        mean_T,
        linewidth=1.6,
        label="Mean temperature"
    )

    ax.plot(
        t_plot,
        median_T,
        linewidth=1.6,
        label="Median temperature"
    )

    ax.plot(
        t_plot,
        max_T,
        linewidth=1.0,
        alpha=0.8,
        label="Maximum channel temperature"
    )

    # Previous TR estimates, not independently
    # detected from surface temperatures.
    for cell in range(1, 6):

        event_time = tr_events.get(
            (chemistry, cell)
        )

        if (
            event_time is not None
            and np.isfinite(event_time)
        ):

            ax.axvline(
                event_time,
                linestyle="--",
                linewidth=1,
                alpha=0.5
            )

            ax.text(
                event_time,
                ax.get_ylim()[1] * 0.95,
                f"C{cell}",
                rotation=90,
                ha="right",
                va="top",
                fontsize=8
            )

    ax.set_title(
        f"{chemistry} - Surface Temperature Envelope"
    )

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [°C]")
    ax.grid(True, alpha=0.25)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_temperature_envelope_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    # --------------------------------------------------------
    # FIGURE 2 - FIRST THRESHOLD CROSSING
    # --------------------------------------------------------

    channel_numbers = [
        int(
            col.split("Temperature ")[1].split(" ")[0]
        )
        for col in temp_columns
    ]

    fig, ax = plt.subplots(
        figsize=(14, 7)
    )

    for threshold in THRESHOLDS:

        column = (
            f"First crossing {int(threshold)}C [s]"
        )

        ax.scatter(
            channel_numbers,
            channel_df[column],
            s=16,
            alpha=0.8,
            label=f"{int(threshold)}°C"
        )

    ax.set_title(
        f"{chemistry} - First Surface Temperature Threshold Crossings"
    )

    ax.set_xlabel("Temperature channel number")
    ax.set_ylabel("First crossing time [s]")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_threshold_crossings_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    # --------------------------------------------------------
    # FIGURE 3 - MAXIMUM HEATING RATE
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(14, 6)
    )

    ax.plot(
        channel_numbers,
        channel_df[
            "Max smoothed heating rate [C/s]"
        ],
        marker=".",
        linewidth=1.2
    )

    ax.set_title(
        f"{chemistry} - Maximum Smoothed Surface Heating Rate"
    )

    ax.set_xlabel("Temperature channel number")
    ax.set_ylabel("Maximum heating rate [°C/s]")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_heating_rates_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

# ============================================================
# SAVE RESULTS
# ============================================================

channel_summary = pd.DataFrame(
    all_channel_rows
)

chemistry_summary = pd.DataFrame(
    all_chemistry_rows
)

channel_summary.to_csv(
    RESULTS_DIR
    / "surface_temporal_channel_summary.csv",
    index=False
)

chemistry_summary.to_csv(
    RESULTS_DIR
    / "surface_temporal_chemistry_summary.csv",
    index=False
)

print("\n" + "=" * 75)
print("FINAL COMPARISON")
print("=" * 75)

print(
    chemistry_summary.round(2).to_string(index=False)
)

print("\nSAVED FILES:")
print("surface_temporal_channel_summary.csv")
print("surface_temporal_chemistry_summary.csv")

for chemistry in CHEMISTRIES:

    print(
        f"surface_temperature_envelope_{chemistry}.png"
    )

    print(
        f"surface_threshold_crossings_{chemistry}.png"
    )

    print(
        f"surface_heating_rates_{chemistry}.png"
    )

print("\nStep 20 completed.")
