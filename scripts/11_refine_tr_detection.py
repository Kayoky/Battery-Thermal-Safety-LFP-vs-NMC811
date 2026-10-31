
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 1. PROJECT SETTINGS
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

TEMP_DIR = ROOT / "data_cleaned"
VOLT_DIR = ROOT / "data" / "voltage"
RESULTS = ROOT / "results"

RESULTS.mkdir(exist_ok=True)

DT = 0.1
SMOOTHING_SECONDS = 1.1
TEMP_THRESHOLD = 10.0

# Diagnostic voltage thresholds.
# These are NOT the paper's TR definition.
VOLTAGE_LEVELS = [0.95, 0.50, 0.10]

SENSOR_MAPPING = {
    1: [1, 2],
    2: [2, 3],
    3: [3, 4],
    4: [4, 5],
    5: [5, 6],
}

CHEMISTRIES = ["LFP", "NMC811"]

# ============================================================
# 2. SIGNAL PROCESSING
# ============================================================

def load_signal(path, column):
    df = pd.read_csv(path)

    t = df["Time [s]"].to_numpy(
        dtype=float, copy=True
    )

    y = df[column].to_numpy(
        dtype=float, copy=True
    )

    return t, y


def interpolate_segments(t, y, grid):
    """
    Interpolate only inside contiguous valid segments.
    Missing measurements remain missing.
    """

    output = np.full(grid.shape, np.nan)

    valid = np.isfinite(t) & np.isfinite(y)

    starts = np.flatnonzero(
        valid & ~np.r_[False, valid[:-1]]
    )

    ends = np.flatnonzero(
        valid & ~np.r_[valid[1:], False]
    )

    for start, end in zip(starts, ends):

        if end <= start:
            continue

        mask = (
            (grid >= t[start])
            & (grid <= t[end])
        )

        output[mask] = np.interp(
            grid[mask],
            t[start:end + 1],
            y[start:end + 1]
        )

    return output


def smooth_segments(values, window):
    """
    Smooth valid segments independently.
    Do not bridge missing-data gaps.
    """

    result = np.full_like(values, np.nan)

    valid = np.isfinite(values)

    starts = np.flatnonzero(
        valid & ~np.r_[False, valid[:-1]]
    )

    ends = np.flatnonzero(
        valid & ~np.r_[valid[1:], False]
    )

    for start, end in zip(starts, ends):

        segment = values[start:end + 1]

        result[start:end + 1] = (
            pd.Series(segment)
            .rolling(
                window=window,
                center=True,
                min_periods=1
            )
            .median()
            .to_numpy()
        )

    return result


def derivative_segments(values, time):

    result = np.full_like(values, np.nan)

    valid = np.isfinite(values)

    starts = np.flatnonzero(
        valid & ~np.r_[False, valid[:-1]]
    )

    ends = np.flatnonzero(
        valid & ~np.r_[valid[1:], False]
    )

    for start, end in zip(starts, ends):

        if end - start < 2:
            continue

        result[start:end + 1] = np.gradient(
            values[start:end + 1],
            time[start:end + 1]
        )

    return result


def threshold_events(time, rate, threshold):
    """
    Return the start time of each separate interval
    where the heating rate exceeds the threshold.
    """

    active = (
        np.isfinite(rate)
        & (rate > threshold)
    )

    starts = np.flatnonzero(
        active & ~np.r_[False, active[:-1]]
    )

    return time[starts]


def first_crossing(time, voltage, threshold):

    valid = (
        np.isfinite(voltage[1:])
        & np.isfinite(voltage[:-1])
    )

    crossings = (
        valid
        & (voltage[:-1] >= threshold)
        & (voltage[1:] < threshold)
        & (time[1:] >= 0)
    )

    indices = np.flatnonzero(crossings)

    if len(indices) == 0:
        return np.nan

    return float(time[indices[0] + 1])


def significant_voltage_drop(time, voltage):
    """
    Find first sustained drop below 95% of the
    initial voltage, not necessarily final collapse.

    Diagnostic criterion, not TR definition.
    """

    baseline_mask = (
        (time >= 0)
        & (time <= 30)
        & np.isfinite(voltage)
    )

    if baseline_mask.sum() == 0:
        return np.nan, np.nan

    baseline = float(
        np.median(voltage[baseline_mask])
    )

    threshold = baseline * 0.95

    below = (
        np.isfinite(voltage)
        & (voltage < threshold)
    )

    required = int(np.ceil(1.0 / DT))

    count = 0

    for i, condition in enumerate(below):

        count = count + 1 if condition else 0

        if count >= required:
            return (
                float(time[i - required + 1]),
                threshold
            )

    return np.nan, threshold


# ============================================================
# 3. ANALYSIS
# ============================================================

summary_rows = []
event_rows = []

WINDOW = max(
    3,
    int(round(SMOOTHING_SECONDS / DT))
)

if WINDOW % 2 == 0:
    WINDOW += 1

for chemistry in CHEMISTRIES:

    print("\n" + "=" * 65)
    print(f"REFINED DETECTION - {chemistry}")
    print("=" * 65)

    for cell in range(1, 6):

        # ----------------------------------------------------
        # Voltage
        # ----------------------------------------------------

        voltage_path = (
            VOLT_DIR
            / chemistry
            / f"U{cell}_{chemistry}.txt"
        )

        t_u, u = load_signal(
            voltage_path,
            "Voltage [V]"
        )

        grid = np.arange(
            max(0.0, t_u[0]),
            t_u[-1],
            DT
        )

        voltage = interpolate_segments(
            t_u, u, grid
        )

        voltage_smooth = smooth_segments(
            voltage, WINDOW
        )

        dU_dt = derivative_segments(
            voltage_smooth, grid
        )

        # Significant voltage decrease
        drop_time, drop_threshold = (
            significant_voltage_drop(
                grid, voltage_smooth
            )
        )

        # Voltage collapse thresholds
        collapse_times = {}

        for level in VOLTAGE_LEVELS:

            collapse_times[level] = first_crossing(
                grid,
                voltage_smooth,
                level
            )

        # ----------------------------------------------------
        # Temperatures
        # ----------------------------------------------------

        temperatures = {}
        heating_rates = {}
        sensor_events = {}

        for sensor in SENSOR_MAPPING[cell]:

            temp_path = (
                TEMP_DIR
                / chemistry
                / f"T{sensor}_{chemistry}.txt"
            )

            t_t, temp = load_signal(
                temp_path,
                "Temperature [°C]"
            )

            # Damaged heater-side sensor
            if chemistry == "LFP" and sensor == 1:
                temp[t_t >= 185.9] = np.nan

            # Unreliable thermocouple:
            # retain exclusion from automated detection
            if chemistry == "NMC811" and sensor == 5:
                temp[:] = np.nan

            temp_interp = interpolate_segments(
                t_t, temp, grid
            )

            temp_smooth = smooth_segments(
                temp_interp, WINDOW
            )

            rate = derivative_segments(
                temp_smooth, grid
            )

            temperatures[sensor] = temp_smooth
            heating_rates[sensor] = rate

            events = threshold_events(
                grid,
                rate,
                TEMP_THRESHOLD
            )

            sensor_events[sensor] = events

            for event_time in events:

                event_rows.append({
                    "Chemistry": chemistry,
                    "Cell": cell,
                    "Sensor": f"T{sensor}",
                    "Temperature event [s]": event_time
                })

        # ----------------------------------------------------
        # Candidate selection
        # ----------------------------------------------------

        all_events = np.concatenate(
            list(sensor_events.values())
        )

        all_events = np.sort(all_events)

        first_temperature_event = (
            float(all_events[0])
            if len(all_events)
            else np.nan
        )

        # Temperature event closest to significant
        # voltage drop, within a diagnostic window.
        # This is a candidate, not confirmed TR.
        matched_event = np.nan

        if (
            np.isfinite(drop_time)
            and len(all_events) > 0
        ):

            distances = np.abs(
                all_events - drop_time
            )

            nearest_index = np.argmin(distances)

            if distances[nearest_index] <= 30.0:
                matched_event = float(
                    all_events[nearest_index]
                )

        # ----------------------------------------------------
        # Save summary
        # ----------------------------------------------------

        summary_rows.append({
            "Chemistry": chemistry,
            "Cell": cell,
            "First temperature event [s]":
                first_temperature_event,
            "Significant voltage drop [s]":
                drop_time,
            "Matched temperature event [s]":
                matched_event,
            "Voltage below 0.95 V [s]":
                collapse_times[0.95],
            "Voltage below 0.50 V [s]":
                collapse_times[0.50],
            "Voltage below 0.10 V [s]":
                collapse_times[0.10],
            "Status": "Candidate - requires validation"
        })

        print(
            f"Cell {cell} | "
            f"First T event: {first_temperature_event:.1f} s | "
            f"Voltage drop: {drop_time:.1f} s | "
            f"Matched T event: {matched_event:.1f} s | "
            f"U < 0.5 V: {collapse_times[0.50]:.1f} s"
        )

        # ----------------------------------------------------
        # Plot
        # ----------------------------------------------------

        fig, axes = plt.subplots(
            3, 1,
            figsize=(13, 10),
            sharex=True
        )

        axes[0].plot(
            grid, voltage_smooth,
            label=f"U{cell}"
        )

        axes[0].set_ylabel("Voltage [V]")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        for sensor, values in temperatures.items():

            axes[1].plot(
                grid, values,
                label=f"T{sensor}"
            )

        axes[1].set_ylabel("Temperature [°C]")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        for sensor, values in heating_rates.items():

            axes[2].plot(
                grid, values,
                label=f"dT{sensor}/dt"
            )

        axes[2].axhline(
            TEMP_THRESHOLD,
            color="red",
            linestyle="--",
            label="10 °C/s"
        )

        axes[2].set_ylabel("Heating rate [°C/s]")
        axes[2].set_xlabel("Time [s]")
        axes[2].legend()
        axes[2].grid(True, alpha=0.3)

        markers = [
            (
                first_temperature_event,
                "orange",
                "First T event"
            ),
            (
                drop_time,
                "purple",
                "Significant voltage drop"
            ),
            (
                collapse_times[0.50],
                "black",
                "Voltage < 0.5 V"
            )
        ]

        for ax in axes:

            for event_time, color, label in markers:

                if np.isfinite(event_time):

                    ax.axvline(
                        event_time,
                        color=color,
                        linestyle="--",
                        linewidth=1.2,
                        label=label
                    )

        fig.suptitle(
            f"Refined TR Diagnostics - {chemistry} Cell {cell}"
        )

        fig.tight_layout()

        output = (
            RESULTS
            / f"refined_detection_{chemistry}_cell{cell}.png"
        )

        fig.savefig(output, dpi=200)
        plt.close(fig)

# ============================================================
# 4. EXPORT
# ============================================================

summary_df = pd.DataFrame(summary_rows)
events_df = pd.DataFrame(event_rows)

summary_path = (
    RESULTS / "refined_tr_detection_summary.csv"
)

events_path = (
    RESULTS / "temperature_event_candidates.csv"
)

summary_df.to_csv(
    summary_path,
    index=False
)

events_df.to_csv(
    events_path,
    index=False
)

print("\n" + "=" * 65)
print("REFINED DETECTION SUMMARY")
print("=" * 65)

print(summary_df.to_string(index=False))

print("\nSaved:")
print(summary_path.name)
print(events_path.name)
print("10 refined diagnostic figures")
