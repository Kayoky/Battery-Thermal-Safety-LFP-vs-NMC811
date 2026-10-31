
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# PATHS
# --------------------------------------------------

project_dir = Path(__file__).resolve().parent.parent
data_dir = project_dir / "data_cleaned"
results_dir = project_dir / "results"
results_dir.mkdir(exist_ok=True)

chemistries = ["LFP", "NMC811"]
thresholds = [200, 300, 500]

# Sensor exclusions based on documented data issues.
# LFP T1: damaged sensor after 185.9 s.
# NMC811 T5: known invalid readings and spikes.
# Keep NMC811 T5 out of quantitative severity metrics
# until its valid intervals have been established.
excluded_channels = {
    "LFP": set(),
    "NMC811": {"T5"}
}

rows = []

# --------------------------------------------------
# ANALYSIS
# --------------------------------------------------

for chemistry in chemistries:

    print("\n" + "=" * 60)
    print(f"TEMPERATURE SEVERITY - {chemistry}")
    print("=" * 60)

    folder = data_dir / chemistry

    for file in sorted(folder.glob("T[1-6]_*.txt")):

        channel = file.stem.split("_")[0]
        df = pd.read_csv(file)

        t = df["Time [s]"].to_numpy(
            dtype=float, copy=True
        )
        temp = df["Temperature [°C]"].to_numpy(
            dtype=float, copy=True
        )

        # Apply known sensor restrictions
        if chemistry == "LFP" and channel == "T1":
            temp[t >= 185.9] = np.nan

        if channel in excluded_channels[chemistry]:
            print(f"{channel}: EXCLUDED - sensor anomaly")
            continue

        # Analyze from heater activation (t = 0)
        valid = (
            np.isfinite(t)
            & np.isfinite(temp)
            & (t >= 0)
        )

        t = t[valid]
        temp = temp[valid]

        if len(t) < 2:
            print(f"{channel}: insufficient valid data")
            continue

        # Maximum measured temperature
        idx_max = np.argmax(temp)
        t_max = float(t[idx_max])
        Tmax = float(temp[idx_max])

        # Heating rate with 1-second rolling smoothing
        series = pd.Series(temp)
        dt_median = float(np.median(np.diff(t)))

        window = max(3, round(1.0 / dt_median))
        smooth = series.rolling(
            window=window,
            center=True,
            min_periods=1
        ).mean().to_numpy()

        heating_rate = np.gradient(smooth, t)
        max_rate = float(np.nanmax(heating_rate))

        # Time above thresholds
        # Linear interpolation at threshold crossings
        durations = {}

        for threshold in thresholds:
            duration = 0.0

            for i in range(len(t) - 1):
                a, b = temp[i], temp[i + 1]
                dt = t[i + 1] - t[i]

                if dt <= 0:
                    continue

                if a >= threshold and b >= threshold:
                    duration += dt

                elif a < threshold and b < threshold:
                    continue

                else:
                    fraction = (threshold - a) / (b - a)
                    fraction = np.clip(fraction, 0, 1)

                    if a < threshold <= b:
                        duration += dt * (1 - fraction)
                    elif a >= threshold > b:
                        duration += dt * fraction

            durations[threshold] = duration

        row = {
            "Chemistry": chemistry,
            "Channel": channel,
            "Tmax [C]": Tmax,
            "Time at Tmax [s]": t_max,
            "Max heating rate [C/s]": max_rate,
            "Time above 200C [s]": durations[200],
            "Time above 300C [s]": durations[300],
            "Time above 500C [s]": durations[500],
        }

        rows.append(row)

        print(
            f"{channel} | Tmax = {Tmax:.1f} °C"
            f" | Max dT/dt = {max_rate:.1f} °C/s"
            f" | >200°C = {durations[200]:.1f} s"
            f" | >300°C = {durations[300]:.1f} s"
            f" | >500°C = {durations[500]:.1f} s"
        )

# --------------------------------------------------
# SAVE RESULTS
# --------------------------------------------------

results = pd.DataFrame(rows)

output = results_dir / "temperature_severity_summary.csv"
results.to_csv(output, index=False)

print("\n" + "=" * 60)
print("TEMPERATURE SEVERITY SUMMARY")
print("=" * 60)
print(results.round(2).to_string(index=False))

# --------------------------------------------------
# PLOT 1: MAXIMUM TEMPERATURE
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(11, 6))

for chemistry in chemistries:
    subset = results[results["Chemistry"] == chemistry]
    ax.plot(
        subset["Channel"],
        subset["Tmax [C]"],
        marker="o",
        linewidth=2,
        label=chemistry
    )

ax.set_title("Maximum Measured Temperature by Sensor Channel")
ax.set_xlabel("Temperature sensor channel")
ax.set_ylabel("Maximum temperature [°C]")
ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
fig.savefig(
    results_dir / "temperature_max_comparison.png",
    dpi=300
)
plt.close(fig)

# --------------------------------------------------
# PLOT 2: DURATION ABOVE THRESHOLDS
# --------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, chemistry in zip(axes, chemistries):

    subset = results[results["Chemistry"] == chemistry]
    x = np.arange(len(subset))
    width = 0.25

    for j, threshold in enumerate(thresholds):
        ax.bar(
            x + (j - 1) * width,
            subset[f"Time above {threshold}C [s]"],
            width,
            label=f">{threshold}°C"
        )

    ax.set_title(chemistry)
    ax.set_xticks(x)
    ax.set_xticklabels(subset["Channel"])
    ax.set_xlabel("Temperature sensor channel")
    ax.set_ylabel("Cumulative time above threshold [s]")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

fig.suptitle("Temperature Exposure Duration")
fig.tight_layout()
fig.savefig(
    results_dir / "temperature_exposure_comparison.png",
    dpi=300
)
plt.close(fig)

print("\nSaved:")
print("temperature_severity_summary.csv")
print("temperature_max_comparison.png")
print("temperature_exposure_comparison.png")
print("\nTemperature severity analysis completed.")
