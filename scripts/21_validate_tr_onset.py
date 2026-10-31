from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# STEP 21 - THERMAL RUNAWAY ONSET VALIDATION
# Corrected version
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

TEMPERATURE_DIR = ROOT / "data_cleaned"
VOLTAGE_DIR = ROOT / "data" / "voltage"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

PREVIOUS_RESULTS = RESULTS_DIR / "refined_tr_detection_summary.csv"

CHEMISTRIES = ["LFP", "NMC811"]
CELLS = range(1, 6)

# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

# Search window around the previous TR estimate
SEARCH_BEFORE = 45.0
SEARCH_AFTER = 45.0

# Common interpolation grid
DT = 0.1

# Smoothing windows in seconds
TEMP_SMOOTHING = 0.5
VOLTAGE_SMOOTHING = 0.5

# Original paper temperature-rate criterion
TEMP_RATE_THRESHOLD = 10.0  # K/s

# Additional robustness criteria
# These are NOT claimed to be the paper's original criteria.
MIN_VOLTAGE_DROP = 0.005     # V over confirmation window
CONFIRMATION_WINDOW = 1.0    # s
MIN_TEMP_INCREASE = 5.0      # K over confirmation window

# Reject isolated voltage derivative artifacts
MIN_NEGATIVE_DUDT = 1e-4     # V/s

# For checking suspiciously large heating rates
SUSPICIOUS_TEMP_RATE = 100.0  # K/s

# Temperature sensors adjacent to each cell
CELL_SENSORS = {
    1: [1, 2],
    2: [2, 3],
    3: [3, 4],
    4: [4, 5],
    5: [5, 6],
}

# Known damaged / unreliable sensor regions
SENSOR_END_TIMES = {
    ("LFP", 1): 185.9,
}

# NMC811 T5 has known negative artifacts.
# Negative temperatures are removed below; remaining spikes
# are assessed by the diagnostics.

# ------------------------------------------------------------
# DATA LOADING
# ------------------------------------------------------------

def load_signal(path, value_keyword):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    df = pd.read_csv(path)

    df.columns = [str(c).strip() for c in df.columns]

    time_col = next(
        (c for c in df.columns if "Time" in c),
        None
    )

    value_col = next(
        (c for c in df.columns if value_keyword in c),
        None
    )

    if time_col is None or value_col is None:
        raise ValueError(
            f"Unexpected columns in {path.name}: "
            f"{list(df.columns)}"
        )

    df = df[[time_col, value_col]].copy()
    df.columns = ["time", "value"]

    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    df = df.dropna(subset=["time"])
    df = df.sort_values("time")
    df = df.drop_duplicates(subset="time")

    return df


def load_temperature(chemistry, sensor):
    path = (
        TEMPERATURE_DIR
        / chemistry
        / f"T{sensor}_{chemistry}.txt"
    )

    df = load_signal(path, "Temperature")

    # Negative readings are not physically meaningful here
    df.loc[df["value"] < 0, "value"] = np.nan

    end_time = SENSOR_END_TIMES.get((chemistry, sensor))

    if end_time is not None:
        df.loc[df["time"] > end_time, "value"] = np.nan

    return df


def load_voltage(chemistry, cell):
    path = (
        VOLTAGE_DIR
        / chemistry
        / f"U{cell}_{chemistry}.txt"
    )

    return load_signal(path, "Voltage")


# ------------------------------------------------------------
# PREVIOUS TR ESTIMATES
# ------------------------------------------------------------

def load_previous_estimates():
    if not PREVIOUS_RESULTS.exists():
        raise FileNotFoundError(
            f"Previous TR results not found:\n{PREVIOUS_RESULTS}"
        )

    df = pd.read_csv(PREVIOUS_RESULTS)

    required = ["Chemistry", "Cell"]

    for column in required:
        if column not in df.columns:
            raise ValueError(
                f"Missing column '{column}' in previous results."
            )

    possible_time_columns = [
        "Matched temperature event [s]",
        "Refined TR time [s]",
        "TR time [s]",
        "Estimated TR time [s]",
    ]

    time_column = next(
        (c for c in possible_time_columns if c in df.columns),
        None
    )

    if time_column is None:
        raise ValueError(
            "Cannot identify previous TR time column.\n"
            f"Available columns: {list(df.columns)}"
        )

    estimates = {}

    for _, row in df.iterrows():
        key = (str(row["Chemistry"]).strip(), int(row["Cell"]))
        estimates[key] = float(row[time_column])

    return estimates


# ------------------------------------------------------------
# SIGNAL PROCESSING
# ------------------------------------------------------------

def interpolate_signal(df, grid):
    x = df["time"].to_numpy(dtype=float)
    y = df["value"].to_numpy(dtype=float)

    valid = np.isfinite(x) & np.isfinite(y)

    x = x[valid]
    y = y[valid]

    result = np.full(len(grid), np.nan)

    if len(x) < 2:
        return result

    inside = (grid >= x[0]) & (grid <= x[-1])

    result[inside] = np.interp(
        grid[inside],
        x,
        y
    )

    # Do not interpolate across substantial missing-data gaps
    gaps = np.diff(x)
    large_gaps = np.where(gaps > 1.0)[0]

    for idx in large_gaps:
        result[
            (grid > x[idx]) & (grid < x[idx + 1])
        ] = np.nan

    return result


def smooth_signal(values, seconds):
    window = max(3, int(round(seconds / DT)))

    if window % 2 == 0:
        window += 1

    series = pd.Series(values)

    # Median suppresses isolated spikes
    median = series.rolling(
        window,
        center=True,
        min_periods=max(2, window // 2)
    ).median()

    # Mean smooths remaining small fluctuations
    smoothed = median.rolling(
        window,
        center=True,
        min_periods=max(2, window // 2)
    ).mean()

    return smoothed.to_numpy()


def calculate_derivative(values):
    derivative = np.full_like(values, np.nan)

    valid = np.isfinite(values)

    # Process contiguous valid sections separately
    edges = np.diff(
        np.concatenate(([False], valid, [False])).astype(int)
    )

    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]

    for start, end in zip(starts, ends):
        if end - start >= 3:
            derivative[start:end] = np.gradient(
                values[start:end],
                DT
            )

    return derivative


def future_change(values, seconds):
    shift = max(1, int(round(seconds / DT)))

    result = np.full(len(values), np.nan)

    result[:-shift] = (
        values[shift:] - values[:-shift]
    )

    return result


# ------------------------------------------------------------
# CANDIDATE DETECTION
# ------------------------------------------------------------

def analyze_cell(chemistry, cell, previous_time):
    start = previous_time - SEARCH_BEFORE
    end = previous_time + SEARCH_AFTER

    grid = np.arange(start, end + DT / 2, DT)

    voltage_df = load_voltage(chemistry, cell)

    voltage = interpolate_signal(voltage_df, grid)
    voltage_smooth = smooth_signal(
        voltage,
        VOLTAGE_SMOOTHING
    )

    dudt = calculate_derivative(voltage_smooth)

    voltage_change = future_change(
        voltage_smooth,
        CONFIRMATION_WINDOW
    )

    voltage_valid = (
        np.isfinite(dudt)
        & np.isfinite(voltage_change)
        & (dudt < -MIN_NEGATIVE_DUDT)
        & (voltage_change <= -MIN_VOLTAGE_DROP)
    )

    all_candidates = []
    sensor_results = {}

    for sensor in CELL_SENSORS[cell]:
        temp_df = load_temperature(chemistry, sensor)

        temperature = interpolate_signal(temp_df, grid)

        temperature_smooth = smooth_signal(
            temperature,
            TEMP_SMOOTHING
        )

        dtdt = calculate_derivative(temperature_smooth)

        temperature_change = future_change(
            temperature_smooth,
            CONFIRMATION_WINDOW
        )

        temperature_valid = (
            np.isfinite(dtdt)
            & np.isfinite(temperature_change)
            & (dtdt > TEMP_RATE_THRESHOLD)
            & (temperature_change >= MIN_TEMP_INCREASE)
        )

        simultaneous = temperature_valid & voltage_valid

        candidate_indices = np.where(simultaneous)[0]

        # Group consecutive detections into events
        groups = []

        if len(candidate_indices) > 0:
            current_group = [candidate_indices[0]]

            for idx in candidate_indices[1:]:
                if idx - current_group[-1] <= 2:
                    current_group.append(idx)
                else:
                    groups.append(current_group)
                    current_group = [idx]

            groups.append(current_group)

        for group in groups:
            idx = group[0]

            candidate = {
                "Chemistry": chemistry,
                "Cell": cell,
                "Temperature sensor": sensor,
                "Candidate TR time [s]": float(grid[idx]),
                "Previous TR estimate [s]": previous_time,
                "Difference [s]": float(
                    grid[idx] - previous_time
                ),
                "dT/dt [K/s]": float(dtdt[idx]),
                "dU/dt [V/s]": float(dudt[idx]),
                "Temperature increase [K]": float(
                    temperature_change[idx]
                ),
                "Voltage change [V]": float(
                    voltage_change[idx]
                ),
                "Event duration [s]": float(
                    len(group) * DT
                ),
            }

            all_candidates.append(candidate)

        sensor_results[sensor] = {
            "temperature": temperature_smooth,
            "dtdt": dtdt,
        }

    # --------------------------------------------------------
    # Candidate selection
    # --------------------------------------------------------

    # Select the candidate closest to the previous estimate.
    # Previous estimate is used only as a search reference.
    # This does NOT independently validate the event.

    if all_candidates:
        all_candidates.sort(
            key=lambda c: (
                abs(c["Difference [s]"]),
                -c["Event duration [s]"]
            )
        )

        selected = all_candidates[0].copy()

        if (
            selected["dT/dt [K/s]"]
            > SUSPICIOUS_TEMP_RATE
        ):
            selected["Status"] = (
                "CANDIDATE - HIGH HEATING RATE; REVIEW"
            )
        else:
            selected["Status"] = (
                "CANDIDATE - MANUAL REVIEW REQUIRED"
            )

    else:
        selected = {
            "Chemistry": chemistry,
            "Cell": cell,
            "Previous TR estimate [s]": previous_time,
            "Candidate TR time [s]": np.nan,
            "Difference [s]": np.nan,
            "Temperature sensor": np.nan,
            "dT/dt [K/s]": np.nan,
            "dU/dt [V/s]": np.nan,
            "Temperature increase [K]": np.nan,
            "Voltage change [V]": np.nan,
            "Event duration [s]": np.nan,
            "Status": "NO ROBUST CANDIDATE",
        }

    # --------------------------------------------------------
    # Diagnostic plot
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        3, 1,
        figsize=(13, 10),
        sharex=True
    )

    axes[0].plot(
        grid,
        voltage_smooth,
        color="black",
        label=f"U{cell}"
    )

    axes[0].set_ylabel("Voltage [V]")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    for sensor, signals in sensor_results.items():
        axes[1].plot(
            grid,
            signals["temperature"],
            label=f"T{sensor}"
        )

        axes[2].plot(
            grid,
            signals["dtdt"],
            label=f"dT{sensor}/dt"
        )

    axes[1].set_ylabel("Temperature [°C]")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    axes[2].axhline(
        TEMP_RATE_THRESHOLD,
        color="red",
        linestyle=":",
        label="10 K/s threshold"
    )

    axes[2].set_ylabel("Heating rate [K/s]")
    axes[2].set_xlabel("Time [s]")
    axes[2].legend()
    axes[2].grid(alpha=0.3)

    for ax in axes:
        ax.axvline(
            previous_time,
            color="gray",
            linestyle="--",
            label="Previous estimate"
        )

        if np.isfinite(
            selected["Candidate TR time [s]"]
        ):
            ax.axvline(
                selected["Candidate TR time [s]"],
                color="green",
                linestyle="-.",
                label="Selected candidate"
            )

    fig.suptitle(
        f"{chemistry} - Cell {cell} - TR Onset Validation"
    )

    plt.tight_layout()

    output = (
        RESULTS_DIR
        / f"tr_validation_{chemistry}_cell{cell}.png"
    )

    plt.savefig(output, dpi=180)
    plt.close(fig)

    return selected, all_candidates


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():
    previous_estimates = load_previous_estimates()

    summary_rows = []
    candidate_rows = []

    for chemistry in CHEMISTRIES:
        print("\n" + "=" * 72)
        print(f"TR ONSET VALIDATION - {chemistry}")
        print("=" * 72)

        for cell in CELLS:
            key = (chemistry, cell)

            if key not in previous_estimates:
                print(
                    f"Cell {cell}: previous estimate missing"
                )
                continue

            previous_time = previous_estimates[key]

            selected, candidates = analyze_cell(
                chemistry,
                cell,
                previous_time
            )

            summary_rows.append(selected)
            candidate_rows.extend(candidates)

            candidate_time = selected[
                "Candidate TR time [s]"
            ]

            if np.isfinite(candidate_time):
                candidate_text = f"{candidate_time:.2f} s"
            else:
                candidate_text = "None"

            print(
                f"Cell {cell}: "
                f"previous={previous_time:.2f} s | "
                f"candidate={candidate_text} | "
                f"status={selected['Status']}"
            )

    summary = pd.DataFrame(summary_rows)
    all_candidates = pd.DataFrame(candidate_rows)

    summary_path = (
        RESULTS_DIR / "tr_onset_validation_summary.csv"
    )

    candidates_path = (
        RESULTS_DIR / "tr_onset_all_candidates.csv"
    )

    summary.to_csv(summary_path, index=False)
    all_candidates.to_csv(
        candidates_path,
        index=False
    )

    print("\n" + "=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)

    columns = [
        "Chemistry",
        "Cell",
        "Previous TR estimate [s]",
        "Candidate TR time [s]",
        "Difference [s]",
        "Temperature sensor",
        "dT/dt [K/s]",
        "dU/dt [V/s]",
        "Status",
    ]

    print(summary[columns].to_string(index=False))

    print("\nSaved files:")
    print(summary_path.name)
    print(candidates_path.name)
    print("tr_validation_<chemistry>_cell<1-5>.png")

    print("\nStep 21 diagnostic analysis completed.")


if __name__ == "__main__":
    main()