
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
# 2. Analysis settings
# --------------------------------------------------

sampling_interval = 0.1  # seconds

# Rolling median reduces isolated measurement spikes
smoothing_window = 11  # approximately 1.1 seconds

temperature_threshold = 10.0  # °C/s

# Require the conditions to persist for 0.3 seconds
minimum_duration = 0.3

# Mapping: cell -> adjacent temperature sensors
sensor_mapping = {
    1: [1, 2],
    2: [2, 3],
    3: [3, 4],
    4: [4, 5],
    5: [5, 6]
}

# --------------------------------------------------
# 3. Helper functions
# --------------------------------------------------

def load_signal(file, column):

    df = pd.read_csv(file)

    time = df["Time [s]"].to_numpy(dtype=float)
    values = df[column].to_numpy(dtype=float)

    return time, values


def interpolate_signal(time, values, common_time):

    valid = np.isfinite(time) & np.isfinite(values)

    if valid.sum() < 2:
        return np.full_like(common_time, np.nan)

    return np.interp(
        common_time,
        time[valid],
        values[valid],
        left=np.nan,
        right=np.nan
    )


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


def calculate_derivative(values, time):

    return np.gradient(values, time)


def first_sustained_event(condition, time):

    required_samples = max(
        1,
        int(np.ceil(minimum_duration / sampling_interval))
    )

    count = 0

    for i, active in enumerate(condition):

        if active:
            count += 1
        else:
            count = 0

        if count >= required_samples:
            return time[i - required_samples + 1]

    return np.nan


# --------------------------------------------------
# 4. Main analysis
# --------------------------------------------------

results = []

for chemistry in ["LFP", "NMC811"]:

    print(f"\n{'=' * 60}")
    print(f"THERMAL RUNAWAY DETECTION - {chemistry}")
    print(f"{'=' * 60}")

    # Load all signals
    voltage_signals = {}
    temperature_signals = {}

    for cell in range(1, 6):

        file = (
            voltage_dir
            / chemistry
            / f"U{cell}_{chemistry}.txt"
        )

        voltage_signals[cell] = load_signal(
            file,
            "Voltage [V]"
        )

    for sensor in range(1, 7):

        file = (
            temperature_dir
            / chemistry
            / f"T{sensor}_{chemistry}.txt"
        )

        temperature_signals[sensor] = load_signal(
            file,
            "Temperature [°C]"
        )

    # Shared time range
    all_signals = (
        list(voltage_signals.values())
        + list(temperature_signals.values())
    )

    start_time = max(signal[0][0] for signal in all_signals)
    end_time = min(signal[0][-1] for signal in all_signals)

    common_time = np.arange(
        max(0, start_time),
        end_time,
        sampling_interval
    )

    # Interpolate and smooth voltages
    voltage_rates = {}

    for cell, (time, values) in voltage_signals.items():

        interpolated = interpolate_signal(
            time, values, common_time
        )

        smoothed = smooth_signal(interpolated)

        voltage_rates[cell] = calculate_derivative(
            smoothed, common_time
        )

    # Interpolate and smooth temperatures
    temperature_rates = {}

    for sensor, (time, values) in temperature_signals.items():

        interpolated = interpolate_signal(
            time, values, common_time
        )

        # Exclude reported damaged LFP T1 sensor
        if chemistry == "LFP" and sensor == 1:
            interpolated[common_time >= 185.9] = np.nan

        # Do not use NMC811 T5 for automatic detection
        if chemistry == "NMC811" and sensor == 5:
            interpolated[:] = np.nan

        smoothed = smooth_signal(interpolated)

        temperature_rates[sensor] = calculate_derivative(
            smoothed, common_time
        )

    # Detect TR for each cell
    detection_times = {}

    for cell in range(1, 6):

        dU_dt = voltage_rates[cell]

        adjacent_sensors = sensor_mapping[cell]

        temperature_conditions = []

        for sensor in adjacent_sensors:

            dT_dt = temperature_rates[sensor]

            temperature_conditions.append(
                dT_dt > temperature_threshold
            )

        # At least one adjacent temperature sensor
        temperature_event = np.logical_or.reduce(
            temperature_conditions
        )

        voltage_event = dU_dt < 0

        combined_condition = (
            temperature_event
            & voltage_event
            & (common_time >= 0)
        )

        tr_time = first_sustained_event(
            combined_condition,
            common_time
        )

        detection_times[cell] = tr_time

        print(
            f"Cell {cell}: "
            f"TR onset = {tr_time:.2f} s"
        )

    # Calculate propagation intervals
    for cell in range(1, 6):

        tr_time = detection_times[cell]

        if cell == 1:
            propagation_delay = np.nan
        else:
            previous_time = detection_times[cell - 1]

            if np.isfinite(tr_time) and np.isfinite(previous_time):
                propagation_delay = tr_time - previous_time
            else:
                propagation_delay = np.nan

        results.append({
            "Chemistry": chemistry,
            "Cell": cell,
            "TR onset [s]": tr_time,
            "Propagation delay [s]": propagation_delay
        })

# --------------------------------------------------
# 5. Export results
# --------------------------------------------------

results_df = pd.DataFrame(results)

output_file = (
    results_dir / "thermal_runaway_detection.csv"
)

results_df.to_csv(output_file, index=False)

print("\n" + "=" * 60)
print("THERMAL RUNAWAY RESULTS")
print("=" * 60)

print(results_df.to_string(index=False))

print(f"\nSaved: {output_file.name}")
