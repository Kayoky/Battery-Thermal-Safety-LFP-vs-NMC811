
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 17 - VENTING SEVERITY QUANTIFICATION
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "venting"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]
CELLS = range(1, 6)

THRESHOLDS = [200.0, 300.0, 500.0]
INTEGRAL_THRESHOLDS = [200.0, 300.0]

# Gaps larger than this multiple of median dt
# are not integrated across.
MAX_GAP_FACTOR = 3.0

rows = []


# ============================================================
# EXACT LINEAR-SEGMENT METRICS
# ============================================================

def exposure_metrics(time, temperature, threshold):
    """
    Calculate:
      1. Observed duration above threshold [s]
      2. Integrated temperature exceedance [C*s]

    Assumes linear temperature variation between
    consecutive valid samples.

    Does not integrate across large sampling gaps.
    """

    if len(time) < 2:
        return 0.0, 0.0

    dt_all = np.diff(time)
    positive_dt = dt_all[dt_all > 0]

    if len(positive_dt) == 0:
        return 0.0, 0.0

    median_dt = np.median(positive_dt)
    max_gap = MAX_GAP_FACTOR * median_dt

    duration = 0.0
    integral = 0.0

    for i in range(len(time) - 1):

        t0 = time[i]
        t1 = time[i + 1]

        y0 = temperature[i]
        y1 = temperature[i + 1]

        if not np.isfinite([t0, t1, y0, y1]).all():
            continue

        dt = t1 - t0

        if dt <= 0 or dt > max_gap:
            continue

        z0 = y0 - threshold
        z1 = y1 - threshold

        # Both endpoints above threshold
        if z0 >= 0 and z1 >= 0:

            duration += dt

            integral += (
                (z0 + z1) / 2.0
            ) * dt

        # Both endpoints below threshold
        elif z0 <= 0 and z1 <= 0:
            continue

        # Upward crossing
        elif z0 < 0 and z1 > 0:

            fraction = -z0 / (z1 - z0)
            fraction = np.clip(fraction, 0, 1)

            active_dt = dt * (1 - fraction)

            duration += active_dt
            integral += 0.5 * z1 * active_dt

        # Downward crossing
        elif z0 > 0 and z1 < 0:

            fraction = z0 / (z0 - z1)
            fraction = np.clip(fraction, 0, 1)

            active_dt = dt * fraction

            duration += active_dt
            integral += 0.5 * z0 * active_dt

    return float(duration), float(integral)


# ============================================================
# MAIN ANALYSIS
# ============================================================

for chemistry in CHEMISTRIES:

    print("\n" + "=" * 70)
    print(f"VENTING SEVERITY - {chemistry}")
    print("=" * 70)

    for cell in CELLS:

        file = (
            DATA_DIR
            / chemistry
            / f"V{cell}_{chemistry}.txt"
        )

        if not file.exists():
            raise FileNotFoundError(file)

        df = pd.read_csv(file)

        time = pd.to_numeric(
            df["Time [s]"],
            errors="coerce"
        ).to_numpy(dtype=float)

        temperature = pd.to_numeric(
            df["Temperature [°C]"],
            errors="coerce"
        ).to_numpy(dtype=float)

        # Keep t >= 0, preserve NaNs as gaps.
        # Do not bridge missing data.
        valid_time = np.isfinite(time) & (time >= 0)

        time = time[valid_time]
        temperature = temperature[valid_time]

        if len(time) < 2:
            raise ValueError(
                f"Insufficient data in {file.name}"
            )

        if np.any(np.diff(time) <= 0):
            raise ValueError(
                f"Invalid time order in {file.name}"
            )

        finite = np.isfinite(temperature)

        if not finite.any():
            raise ValueError(
                f"No valid temperatures in {file.name}"
            )

        valid_indices = np.flatnonzero(finite)
        idx_max = valid_indices[
            np.argmax(temperature[finite])
        ]

        Tmax = float(temperature[idx_max])
        t_Tmax = float(time[idx_max])

        durations = {}
        integrals = {}

        for threshold in THRESHOLDS:

            duration, integral = exposure_metrics(
                time,
                temperature,
                threshold
            )

            durations[threshold] = duration

            if threshold in INTEGRAL_THRESHOLDS:
                integrals[threshold] = integral

        row = {
            "Chemistry": chemistry,
            "Cell": cell,
            "Tmax [C]": Tmax,
            "Time at Tmax [s]": t_Tmax,
            "Time above 200C [s]": durations[200.0],
            "Time above 300C [s]": durations[300.0],
            "Time above 500C [s]": durations[500.0],
            "Exposure integral above 200C [C*s]":
                integrals[200.0],
            "Exposure integral above 300C [C*s]":
                integrals[300.0],
        }

        rows.append(row)

        print(
            f"Cell {cell}: "
            f"Tmax={Tmax:.1f} C | "
            f">200C={durations[200.0]:.1f} s | "
            f">300C={durations[300.0]:.1f} s | "
            f">500C={durations[500.0]:.1f} s | "
            f"I200={integrals[200.0]:.1f} C*s"
        )


# ============================================================
# SAVE DETAILED RESULTS
# ============================================================

results = pd.DataFrame(rows)

output_csv = (
    RESULTS_DIR / "venting_severity_summary.csv"
)

results.to_csv(output_csv, index=False)

print("\n" + "=" * 70)
print("VENTING SEVERITY SUMMARY")
print("=" * 70)

print(results.round(2).to_string(index=False))


# ============================================================
# CHEMISTRY SUMMARY
# ============================================================

metrics = [
    "Tmax [C]",
    "Time above 200C [s]",
    "Time above 300C [s]",
    "Time above 500C [s]",
    "Exposure integral above 200C [C*s]",
    "Exposure integral above 300C [C*s]",
]

chemistry_summary = (
    results.groupby("Chemistry")[metrics]
    .agg(["mean", "median", "max"])
)

chemistry_summary.columns = [
    f"{metric} - {stat}"
    for metric, stat in chemistry_summary.columns
]

chemistry_summary.to_csv(
    RESULTS_DIR / "venting_severity_chemistry_summary.csv"
)

print("\nCHEMISTRY SUMMARY:")
print(chemistry_summary.round(2).to_string())


# ============================================================
# FIGURE 1: MAXIMUM TEMPERATURES
# ============================================================

fig, ax = plt.subplots(figsize=(10, 6))

for chemistry in CHEMISTRIES:

    subset = results[
        results["Chemistry"] == chemistry
    ].sort_values("Cell")

    ax.plot(
        subset["Cell"],
        subset["Tmax [C]"],
        marker="o",
        linewidth=2,
        label=chemistry
    )

ax.set_title(
    "Maximum Measured Safety Vent Temperatures"
)
ax.set_xlabel("Cell number")
ax.set_ylabel("Maximum temperature [°C]")
ax.set_xticks(list(CELLS))
ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()

fig.savefig(
    RESULTS_DIR / "venting_severity_max_temperatures.png",
    dpi=300
)

plt.show()
plt.close(fig)


# ============================================================
# FIGURE 2: TIME ABOVE TEMPERATURE THRESHOLDS
# ============================================================

fig, axes = plt.subplots(
    1, 2,
    figsize=(14, 5),
    sharey=True
)

for ax, chemistry in zip(axes, CHEMISTRIES):

    subset = results[
        results["Chemistry"] == chemistry
    ].sort_values("Cell")

    x = np.arange(len(subset))
    width = 0.25

    for j, threshold in enumerate(THRESHOLDS):

        column = f"Time above {int(threshold)}C [s]"

        ax.bar(
            x + (j - 1) * width,
            subset[column],
            width,
            label=f">{int(threshold)}°C"
        )

    ax.set_title(chemistry)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"Cell {i}" for i in subset["Cell"]]
    )
    ax.set_xlabel("Cell")
    ax.set_ylabel("Observed duration [s]")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

fig.suptitle(
    "Observed Vent Temperature Exposure Duration"
)
fig.tight_layout()

fig.savefig(
    RESULTS_DIR / "venting_severity_threshold_durations.png",
    dpi=300
)

plt.show()
plt.close(fig)


# ============================================================
# FIGURE 3: THERMAL EXPOSURE INTEGRALS
# ============================================================

fig, axes = plt.subplots(
    1, 2,
    figsize=(14, 5),
    sharey=True
)

for ax, chemistry in zip(axes, CHEMISTRIES):

    subset = results[
        results["Chemistry"] == chemistry
    ].sort_values("Cell")

    x = np.arange(len(subset))
    width = 0.35

    for j, threshold in enumerate(INTEGRAL_THRESHOLDS):

        column = (
            f"Exposure integral above "
            f"{int(threshold)}C [C*s]"
        )

        ax.bar(
            x + (j - 0.5) * width,
            subset[column],
            width,
            label=f"Above {int(threshold)}°C"
        )

    ax.set_title(chemistry)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"Cell {i}" for i in subset["Cell"]]
    )
    ax.set_xlabel("Cell")
    ax.set_ylabel("Temperature exposure [°C·s]")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

fig.suptitle(
    "Integrated Vent Temperature Exceedance"
)
fig.tight_layout()

fig.savefig(
    RESULTS_DIR / "venting_severity_exposure_integrals.png",
    dpi=300
)

plt.show()
plt.close(fig)


print("\n" + "=" * 70)
print("SAVED FILES")
print("=" * 70)

print("venting_severity_summary.csv")
print("venting_severity_chemistry_summary.csv")
print("venting_severity_max_temperatures.png")
print("venting_severity_threshold_durations.png")
print("venting_severity_exposure_integrals.png")

print("\nStep 17 completed.")
