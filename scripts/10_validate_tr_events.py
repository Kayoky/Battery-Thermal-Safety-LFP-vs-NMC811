
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# 1. Directories and settings
# --------------------------------------------------

project_dir = Path(__file__).resolve().parent.parent
temp_dir = project_dir / "data_cleaned"
voltage_dir = project_dir / "data" / "voltage"
results_dir = project_dir / "results"

results_dir.mkdir(exist_ok=True)

dt = 0.1
smoothing_window = 11
temperature_threshold = 10.0

sensor_mapping = {
    1: [1, 2],
    2: [2, 3],
    3: [3, 4],
    4: [4, 5],
    5: [5, 6]
}

# Diagnostic windows around candidate events
window_before = 120  # seconds
window_after = 180   # seconds

# --------------------------------------------------
# 2. Helper functions
# --------------------------------------------------

def load_signal(path, column):
    df = pd.read_csv(path)

    return (
        df["Time [s]"].to_numpy(dtype=float, copy=True),
        df[column].to_numpy(dtype=float, copy=True)
    )

def interpolate_without_gaps(t, values, grid):

    result = np.full(grid.shape, np.nan, dtype=float)

    valid = np.isfinite(t) & np.isfinite(values)

    # Find contiguous valid segments
    starts = np.flatnonzero(valid & ~np.r_[False, valid[:-1]])
    ends = np.flatnonzero(valid & ~np.r_[valid[1:], False])

    for start, end in zip(starts, ends):
        if end <= start:
            continue

        segment = (
            (grid >= t[start])
            & (grid <= t[end])
        )

        result[segment] = np.interp(
            grid[segment],
            t[start:end + 1],
            values[start:end + 1]
        )

    return result


def smooth_signal(values):
    return (
        pd.Series(values)
        .rolling(
            window=smoothing_window,
            center=True,
            min_periods=smoothing_window
        )
        .median()
        .to_numpy()
    )


def derivative(values, time):
    result = np.full_like(values, np.nan, dtype=float)
    valid = np.isfinite(values)

    starts = np.flatnonzero(valid & ~np.r_[False, valid[:-1]])
    ends = np.flatnonzero(valid & ~np.r_[valid[1:], False])

    for start, end in zip(starts, ends):
        if end - start < 2:
            continue

        result[start:end + 1] = np.gradient(
            values[start:end + 1],
            time[start:end + 1]
        )

    return result


def voltage_collapse(t, voltage):

    # Diagnostic threshold, not the paper's TR criterion
    baseline = np.nanmedian(
        voltage[(t >= 0) & (t <= 30)]
    )

    if not np.isfinite(baseline):
        return np.nan

    threshold = 0.5

    # First downward crossing below 0.5 V
    crossing = (
        (voltage[1:] < threshold)
        & (voltage[:-1] >= threshold)
        & (t[1:] >= 0)
    )

    indices = np.flatnonzero(crossing)

    return (
        t[indices[0] + 1]
        if len(indices) else np.nan
    )


# --------------------------------------------------
# 3. Load previous diagnostic candidates
# --------------------------------------------------

summary_file = (
    results_dir / "detection_diagnostics_summary.csv"
)

summary = pd.read_csv(summary_file)

validation_rows = []

# --------------------------------------------------
# 4. Main analysis
# --------------------------------------------------

for chemistry in ["LFP", "NMC811"]:

    print(f"\n{'=' * 60}")
    print(f"VALIDATION - {chemistry}")
    print(f"{'=' * 60}")

    for cell in range(1, 6):

        # Previous candidate from script 09
        previous = summary[
            (summary["Chemistry"] == chemistry)
            & (summary["Cell"] == cell)
        ]

        if previous.empty:
            print(f"Missing diagnostic result: {chemistry} Cell {cell}")
            continue

        candidate = float(
            previous["Combined detection [s]"].iloc[0]
        )

        # Load voltage on its full time range
        voltage_file = (
            voltage_dir / chemistry /
            f"U{cell}_{chemistry}.txt"
        )

        t_u, u = load_signal(
            voltage_file, "Voltage [V]"
        )

        time = np.arange(
            max(0, t_u[0]),
            t_u[-1],
            dt
        )

        voltage = interpolate_without_gaps(
            t_u, u, time
        )

        voltage_smooth = smooth_signal(voltage)
        dU_dt = derivative(voltage_smooth, time)

        collapse_time = voltage_collapse(
            time, voltage_smooth
        )

        # Load adjacent temperatures
        temperatures = {}
        temperature_rates = {}

        for sensor in sensor_mapping[cell]:

            temp_file = (
                temp_dir / chemistry /
                f"T{sensor}_{chemistry}.txt"
            )

            t_t, temp = load_signal(
                temp_file, "Temperature [°C]"
            )

            # Documented sensor limitation
            if chemistry == "LFP" and sensor == 1:
                temp[t_t >= 185.9] = np.nan

            # Unreliable channel excluded from TR detection
            if chemistry == "NMC811" and sensor == 5:
                temp[:] = np.nan

            interpolated = interpolate_without_gaps(
                t_t, temp, time
            )

            smoothed = smooth_signal(interpolated)

            temperatures[sensor] = smoothed
            temperature_rates[sensor] = derivative(
                smoothed, time
            )

        # Select a window around the previous candidate
        center = (
            candidate
            if np.isfinite(candidate)
            else collapse_time
        )

        if not np.isfinite(center):
            print(f"No candidate: {chemistry} Cell {cell}")
            continue

        left = max(0, center - window_before)
        right = center + window_after

        mask = (time >= left) & (time <= right)

        # ------------------------------------------
        # Plot: voltage, temperatures, derivatives
        # ------------------------------------------

        fig, axes = plt.subplots(
            3, 1,
            figsize=(13, 10),
            sharex=True
        )

        axes[0].plot(
            time[mask],
            voltage_smooth[mask],
            label=f"U{cell}"
        )

        axes[0].set_ylabel("Voltage [V]")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        for sensor, values in temperatures.items():
            axes[1].plot(
                time[mask],
                values[mask],
                label=f"T{sensor}"
            )

        axes[1].set_ylabel("Temperature [°C]")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        for sensor, rate in temperature_rates.items():
            axes[2].plot(
                time[mask],
                rate[mask],
                label=f"dT{sensor}/dt"
            )

        axes[2].axhline(
            temperature_threshold,
            color="red",
            linestyle="--",
            label="10 °C/s threshold"
        )

        axes[2].set_ylabel("Heating rate [°C/s]")
        axes[2].set_xlabel("Time [s]")
        axes[2].legend()
        axes[2].grid(True, alpha=0.3)

        for ax in axes:

            if np.isfinite(candidate):
                ax.axvline(
                    candidate,
                    color="purple",
                    linestyle="--",
                    linewidth=1.5,
                    label="Combined candidate"
                )

            if np.isfinite(collapse_time):
                ax.axvline(
                    collapse_time,
                    color="black",
                    linestyle=":",
                    linewidth=1.5,
                    label="Voltage collapse"
                )

        fig.suptitle(
            f"TR Event Validation - {chemistry} Cell {cell}"
        )

        fig.tight_layout()

        output_file = (
            results_dir /
            f"validation_{chemistry}_cell{cell}.png"
        )

        fig.savefig(output_file, dpi=200)
        plt.close(fig)

        validation_rows.append({
            "Chemistry": chemistry,
            "Cell": cell,
            "Combined candidate [s]": candidate,
            "Voltage collapse [s]": collapse_time,
            "Difference [s]": (
                collapse_time - candidate
                if np.isfinite(candidate)
                and np.isfinite(collapse_time)
                else np.nan
            ),
            "Status": "Requires visual validation"
        })

        print(
            f"Cell {cell}: "
            f"Candidate = {candidate:.1f} s | "
            f"Collapse = {collapse_time:.1f} s"
        )

# --------------------------------------------------
# 5. Export validation summary
# --------------------------------------------------

validation_df = pd.DataFrame(validation_rows)

output_file = results_dir / "tr_validation_summary.csv"

validation_df.to_csv(output_file, index=False)

print("\n" + "=" * 60)
print("VALIDATION SUMMARY")
print("=" * 60)

print(validation_df.to_string(index=False))

print(f"\nSaved: {output_file.name}")
print("Saved: 10 validation figures")
