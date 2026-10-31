
from pathlib import Path
import pandas as pd
import numpy as np

# Project directories
project_dir = Path(__file__).resolve().parent.parent
data_dir = project_dir / "data_cleaned"
results_dir = project_dir / "results"

# Preliminary threshold for unusually rapid changes
threshold = 100.0  # °C/s

summary = []

for chemistry in ["LFP", "NMC811"]:

    print(f"\n{'=' * 55}")
    print(f"Battery chemistry: {chemistry}")
    print(f"{'=' * 55}")

    folder = data_dir / chemistry

    for file in sorted(folder.glob("T[1-6]_*.txt")):

        df = pd.read_csv(file)

        time = df["Time [s]"]
        temperature = df["Temperature [°C]"]

        # Calculate time and temperature differences
        dt = time.diff()
        dT = temperature.diff()

        # Temperature change rate
        heating_rate = dT / dt

        # Only compare consecutive valid measurements
        valid = (
            dt.gt(0)
            & temperature.notna()
            & temperature.shift(1).notna()
        )

        # Identify unusually large changes
        suspicious = valid & (heating_rate.abs() > threshold)

        # Separate heating and cooling events
        rapid_heating = valid & (heating_rate > threshold)
        rapid_cooling = valid & (heating_rate < -threshold)

        sensor = file.stem.split("_")[0]

        summary.append({
            "Chemistry": chemistry,
            "Sensor": sensor,
            "Suspicious jumps": int(suspicious.sum()),
            "Rapid heating": int(rapid_heating.sum()),
            "Rapid cooling": int(rapid_cooling.sum())
        })

        print(
            f"{sensor}: "
            f"{suspicious.sum()} suspicious jumps | "
            f"Heating: {rapid_heating.sum()} | "
            f"Cooling: {rapid_cooling.sum()}"
        )

        if suspicious.any():

            examples = pd.DataFrame({
                "Time [s]": time[suspicious],
                "Temperature [°C]": temperature[suspicious],
                "Rate [°C/s]": heating_rate[suspicious]
            })

            print(examples.head(5).to_string(index=False))

# Export summary
summary_df = pd.DataFrame(summary)

output_file = results_dir / "temperature_anomaly_summary.csv"
summary_df.to_csv(output_file, index=False)

print(f"\nSummary saved: {output_file}")
