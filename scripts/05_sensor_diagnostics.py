
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Project directories
project_dir = Path(__file__).resolve().parent.parent
raw_dir = project_dir / "data"
cleaned_dir = project_dir / "data_cleaned"
results_dir = project_dir / "results"

results_dir.mkdir(exist_ok=True)

# Diagnostic threshold (not a TR criterion)
threshold = 100.0  # °C/s

# Sensor configurations
sensors = [
    {
        "chemistry": "NMC811",
        "sensor": "T5",
        "time_limit": 650,
        "reference_time": None
    },
    {
        "chemistry": "LFP",
        "sensor": "T1",
        "time_limit": 500,
        "reference_time": 185.9
    }
]

for config in sensors:

    chemistry = config["chemistry"]
    sensor = config["sensor"]
    time_limit = config["time_limit"]

    filename = f"{sensor}_{chemistry}.txt"

    # Load original and cleaned measurements
    raw = pd.read_csv(raw_dir / chemistry / filename)
    cleaned = pd.read_csv(cleaned_dir / chemistry / filename)

    time = raw["Time [s]"]
    raw_temp = raw["Temperature [°C]"]
    clean_temp = cleaned["Temperature [°C]"]

    # Calculate temperature change rate
    dt = time.diff()
    heating_rate = clean_temp.diff() / dt

    valid = (
        dt.gt(0)
        & clean_temp.notna()
        & clean_temp.shift(1).notna()
    )

    suspicious = valid & (heating_rate.abs() > threshold)

    # Identify values removed during cleaning
    invalid = raw_temp.notna() & clean_temp.isna()

    # Restrict the displayed time interval
    window = (time >= 0) & (time <= time_limit)

    fig, ax = plt.subplots(figsize=(13, 6))

    # Original temperature signal
    ax.plot(
        time[window],
        raw_temp[window],
        label="Raw temperature",
        color="gray",
        alpha=0.55,
        linewidth=1
    )

    # Cleaned temperature signal
    ax.plot(
        time[window],
        clean_temp[window],
        label="Cleaned temperature",
        color="blue",
        linewidth=1.3
    )

    # Mark invalid measurements
    invalid_window = invalid & window

    ax.scatter(
        time[invalid_window],
        raw_temp[invalid_window],
        color="red",
        s=10,
        label="Removed invalid values",
        zorder=4
    )

    # Mark rapid temperature changes
    suspicious_window = suspicious & window

    ax.scatter(
        time[suspicious_window],
        clean_temp[suspicious_window],
        color="orange",
        s=24,
        label="|dT/dt| > 100 °C/s",
        zorder=5
    )

    # Reference from scientific publication
    reference_time = config["reference_time"]

    if reference_time is not None:
        ax.axvline(
            reference_time,
            color="purple",
            linestyle="--",
            linewidth=1.5,
            label=f"Reported sensor issue: {reference_time} s"
        )

    ax.set_title(
        f"Sensor Diagnostics - {chemistry} {sensor}"
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [°C]")

    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    # Save figure
    output_file = (
        results_dir / f"diagnostics_{chemistry}_{sensor}.png"
    )

    fig.savefig(output_file, dpi=300)

    print(f"\n--- {chemistry} {sensor} ---")
    print(f"Invalid values: {invalid_window.sum()}")
    print(f"Suspicious jumps: {suspicious_window.sum()}")
    print(f"Saved: {output_file.name}")

# Show both figures
plt.show()
