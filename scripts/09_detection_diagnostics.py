
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# 1. Project directories
# --------------------------------------------------

project_dir = Path(__file__).resolve().parent.parent

temperature_dir = project_dir / "data_cleaned"
voltage_dir = project_dir / "data" / "voltage"
results_dir = project_dir / "results"

results_dir.mkdir(exist_ok=True)

# --------------------------------------------------
# 2. Settings
# --------------------------------------------------

dt = 0.1
smoothing_window = 11
temperature_threshold = 10.0
minimum_duration = 0.3

sensor_mapping = {
    1: [1, 2],
    2: [2, 3],
    3: [3, 4],
    4: [4, 5],
    5: [5, 6]
}

# --------------------------------------------------
# 3. Functions
# --------------------------------------------------

def load_signal(file, column):
    df = pd.read_csv(file)

    return (
        df["Time [s]"].to_numpy(dtype=float),
        df[column].to_numpy(dtype=float)
    )


def interpolate_signal(time, values, common_time):

    valid = np.isfinite(time) & np.isfinite(values)

    if valid.sum() < 2:
        return np.full_like(common_time, np.nan)

    result = np.interp(
        common_time,
        time[valid],
        values[valid],
        left=np.nan,
        right=np.nan
    )

    # Avoid interpolation across missing-data gaps
    invalid_indices = np.where(~np.isfinite(values))[0]

    if len(invalid_indices) > 0:
        for index in invalid_indices:
            if 0 < index < len(time) - 1:
                gap = (
                    (common_time >= time[index - 1])
                    & (common_time <= time[index + 1])
                )
                result[gap] = np.nan

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


def first_sustained_event(condition, time):

    required_samples = max(
        1,
        int(np.ceil(minimum_duration / dt))
    )

    count = 0

    for i, active in enumerate(condition):

        count = count + 1 if active else 0

        if count >= required_samples:
            return time[i - required_samples + 1]

    return np.nan


def first_voltage_collapse(time, voltage):

    # Diagnostic definition:
    # First crossing below 0.5 V after initially
    # being above 2 V.
    if not np.any(voltage > 2.0):
        return np.nan

    below = np.isfinite(voltage) & (voltage < 0.5)

    crossings = (
        below
        & np.roll(~below, 1)
        & (time >= 0)
    )
    crossings[0] = False

    indices = np.where(crossings)[0]

    if len(indices) == 0:
        return np.nan

    return time[indices[0]]


# --------------------------------------------------
# 4. Main analysis
# --------------------------------------------------

summary = []

for chemistry in ["LFP", "NMC811"]:

    print(f"\n{'=' * 65}")
    print(f"DETECTION DIAGNOSTICS - {chemistry}")
    print(f"{'=' * 65}")

    for cell in range(1, 6):

        # Load voltage
        voltage_file = (
            voltage_dir / chemistry /
            f"U{cell}_{chemistry}.txt"
        )

        t_u, u = load_signal(
            voltage_file, "Voltage [V]"
        )

        # Load adjacent temperature sensors
        sensors = sensor_mapping[cell]
        temperature_data = {}

        for sensor in sensors:

            temp_file = (
                temperature_dir / chemistry /
                f"T{sensor}_{chemistry}.txt"
            )

            temperature_data[sensor] = load_signal(
                temp_file, "Temperature [°C]"
            )

        # Shared time axis
        start = max(
            0,
            t_u[0],
            *(temperature_data[s][0][0] for s in sensors)
        )

        end = min(
            t_u[-1],
            *(temperature_data[s][0][-1] for s in sensors)
        )

        time = np.arange(start, end, dt)

        voltage = interpolate_signal(t_u, u, time)
        voltage_smooth = smooth_signal(voltage)
        dU_dt = np.gradient(voltage_smooth, time)

        # Temperature derivatives
        temperature_rates = []

        for sensor in sensors:

            t_t, temperature = temperature_data[sensor]

            temp = interpolate_signal(
                t_t, temperature, time
            )

            if chemistry == "LFP" and sensor == 1:
                temp[time >= 185.9] = np.nan

            if chemistry == "NMC811" and sensor == 5:
                temp[:] = np.nan

            temp_smooth = smooth_signal(temp)

            dT_dt = np.gradient(temp_smooth, time)
            temperature_rates.append(dT_dt)

        # Use the largest valid heating rate
        rates = np.vstack(temperature_rates)

        max_dT_dt = np.full(len(time), np.nan)

        for i in range(len(time)):
            valid_rates = rates[:, i]
            valid_rates = valid_rates[
                np.isfinite(valid_rates)
            ]

            if len(valid_rates) > 0:
                max_dT_dt[i] = valid_rates.max()

        # Independent conditions
        temperature_condition = (
            max_dT_dt > temperature_threshold
        )

        voltage_condition = dU_dt < 0

        combined_condition = (
            temperature_condition & voltage_condition
        )

        # Event times
        t_temperature = first_sustained_event(
            temperature_condition, time
        )

        t_voltage_negative = first_sustained_event(
            voltage_condition, time
        )

        t_combined = first_sustained_event(
            combined_condition, time
        )

        t_voltage_collapse = first_voltage_collapse(
            time, voltage_smooth
        )

        # Maximum heating rate
        finite_rates = max_dT_dt[np.isfinite(max_dT_dt)]

        max_heating_rate = (
            finite_rates.max()
            if len(finite_rates) > 0
            else np.nan
        )

        summary.append({
            "Chemistry": chemistry,
            "Cell": cell,
            "First dT/dt > 10 [s]": t_temperature,
            "First dU/dt < 0 [s]": t_voltage_negative,
            "Combined detection [s]": t_combined,
            "Voltage collapse [s]": t_voltage_collapse,
            "Max heating rate [C/s]": max_heating_rate
        })

        print(
            f"Cell {cell} | "
            f"Temperature: {t_temperature:.1f} s | "
            f"Voltage decrease: {t_voltage_negative:.1f} s | "
            f"Combined: {t_combined:.1f} s | "
            f"Collapse: {t_voltage_collapse:.1f} s | "
            f"Max dT/dt: {max_heating_rate:.1f} °C/s"
        )

        # ------------------------------------------
        # Diagnostic plot
        # ------------------------------------------

        fig, axes = plt.subplots(
            3, 1,
            figsize=(12, 9),
            sharex=True
        )

        axes[0].plot(time, voltage_smooth)
        axes[0].set_ylabel("Voltage [V]")
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(time, dU_dt)
        axes[1].axhline(
            0, color="red", linestyle="--"
        )
        axes[1].set_ylabel("dU/dt [V/s]")
        axes[1].grid(True, alpha=0.3)

        axes[2].plot(time, max_dT_dt)
        axes[2].axhline(
            temperature_threshold,
            color="red",
            linestyle="--",
            label="10 °C/s threshold"
        )
        axes[2].set_ylabel("Max dT/dt [°C/s]")
        axes[2].set_xlabel("Time [s]")
        axes[2].legend()
        axes[2].grid(True, alpha=0.3)

        fig.suptitle(
            f"TR Detection Diagnostics - {chemistry} Cell {cell}"
        )

        fig.tight_layout()

        output_file = (
            results_dir /
            f"detection_diagnostics_{chemistry}_cell{cell}.png"
        )

        fig.savefig(output_file, dpi=200)
        plt.close(fig)

# --------------------------------------------------
# 5. Export summary
# --------------------------------------------------

summary_df = pd.DataFrame(summary)

output_file = (
    results_dir / "detection_diagnostics_summary.csv"
)

summary_df.to_csv(output_file, index=False)

print("\n" + "=" * 65)
print("DIAGNOSTIC SUMMARY")
print("=" * 65)

print(summary_df.to_string(index=False))

print(f"\nSaved: {output_file.name}")
print("Saved: 10 diagnostic figures")
