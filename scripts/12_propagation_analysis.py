
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 12 - THERMAL RUNAWAY PROPAGATION ANALYSIS
# ============================================================

project_dir = Path(__file__).resolve().parent.parent
results_dir = project_dir / "results"
results_dir.mkdir(parents=True, exist_ok=True)

input_file = results_dir / "refined_tr_detection_summary.csv"

# Expected output from script 11
time_column = "Matched temperature event [s]"

if not input_file.exists():
    raise FileNotFoundError(
        f"Missing input file: {input_file}"
    )

df = pd.read_csv(input_file)

required_columns = ["Chemistry", "Cell", time_column]

for column in required_columns:
    if column not in df.columns:
        raise ValueError(
            f"Missing column: {column}\n"
            f"Available columns: {list(df.columns)}"
        )

df["Cell"] = pd.to_numeric(df["Cell"], errors="raise")
df[time_column] = pd.to_numeric(
    df[time_column], errors="coerce"
)

chemistries = ["LFP", "NMC811"]

propagation_rows = []
summary_rows = []
cell_times = {}

# ============================================================
# 1. PROPAGATION CALCULATIONS
# ============================================================

for chemistry in chemistries:

    subset = df[df["Chemistry"] == chemistry].copy()
    subset = subset.sort_values("Cell")

    if subset["Cell"].tolist() != [1, 2, 3, 4, 5]:
        raise ValueError(
            f"{chemistry}: Expected exactly cells 1-5."
        )

    if subset[time_column].isna().any():
        raise ValueError(
            f"{chemistry}: Missing TR onset estimates."
        )

    times = subset[time_column].to_numpy(dtype=float)
    cell_times[chemistry] = times

    delays = np.diff(times)

    if np.any(delays <= 0):
        raise ValueError(
            f"{chemistry}: TR onset times are not increasing."
        )

    total_time = times[-1] - times[0]
    mean_delay = np.mean(delays)
    median_delay = np.median(delays)
    std_delay = np.std(delays, ddof=1)

    print("\n" + "=" * 60)
    print(f"PROPAGATION ANALYSIS - {chemistry}")
    print("=" * 60)

    for i, delay in enumerate(delays, start=1):
        print(
            f"Cell {i} -> Cell {i+1}: "
            f"{delay:.1f} s"
        )

        propagation_rows.append({
            "Chemistry": chemistry,
            "From cell": i,
            "To cell": i + 1,
            "Propagation delay [s]": delay,
            "Method": "Refined matched TR onset estimate"
        })

    print(f"\nTotal propagation time: {total_time:.1f} s")
    print(f"Mean propagation delay: {mean_delay:.1f} s")
    print(f"Median propagation delay: {median_delay:.1f} s")
    print(f"Standard deviation: {std_delay:.1f} s")

    summary_rows.append({
        "Chemistry": chemistry,
        "First TR onset [s]": times[0],
        "Last TR onset [s]": times[-1],
        "Total propagation time [s]": total_time,
        "Mean delay [s]": mean_delay,
        "Median delay [s]": median_delay,
        "Std delay [s]": std_delay,
        "Number of propagation intervals": len(delays),
        "Method": "Refined matched TR onset estimate"
    })

propagation_df = pd.DataFrame(propagation_rows)
summary_df = pd.DataFrame(summary_rows)

# ============================================================
# 2. COMPARISON
# ============================================================

lfp_total = summary_df.loc[
    summary_df["Chemistry"] == "LFP",
    "Total propagation time [s]"
].iloc[0]

nmc_total = summary_df.loc[
    summary_df["Chemistry"] == "NMC811",
    "Total propagation time [s]"
].iloc[0]

ratio = lfp_total / nmc_total

print("\n" + "=" * 60)
print("LFP VS NMC811 COMPARISON")
print("=" * 60)

print(f"LFP total time: {lfp_total:.1f} s")
print(f"NMC811 total time: {nmc_total:.1f} s")
print(f"Time ratio LFP/NMC811: {ratio:.2f}")

print(
    f"\nNMC811 propagation occurred over a "
    f"{ratio:.2f} times shorter total time interval."
)

# ============================================================
# 3. SAVE TABLES
# ============================================================

propagation_file = results_dir / "propagation_delays.csv"
summary_file = results_dir / "propagation_summary.csv"

propagation_df.to_csv(propagation_file, index=False)
summary_df.to_csv(summary_file, index=False)

print("\nSaved:")
print(propagation_file.name)
print(summary_file.name)

# ============================================================
# 4. FIGURE 1 - TR ONSET TIMELINE
# ============================================================

fig1, ax1 = plt.subplots(figsize=(11, 6))

for chemistry in chemistries:
    times = cell_times[chemistry]

    relative_times = times - times[0]

    ax1.plot(
        range(1, 6),
        relative_times,
        marker="o",
        linewidth=2,
        markersize=7,
        label=chemistry
    )

ax1.set_title("Thermal Runaway Propagation Timeline")
ax1.set_xlabel("Cell number")
ax1.set_ylabel("Time since first TR onset [s]")
ax1.set_xticks(range(1, 6))
ax1.grid(True, alpha=0.3)
ax1.legend()

fig1.tight_layout()

fig1.savefig(
    results_dir / "tr_propagation_timeline.png",
    dpi=300
)

# ============================================================
# 5. FIGURE 2 - PROPAGATION DELAYS
# ============================================================

fig2, ax2 = plt.subplots(figsize=(11, 6))

x = np.arange(4)
width = 0.35

for index, chemistry in enumerate(chemistries):

    delays = propagation_df.loc[
        propagation_df["Chemistry"] == chemistry,
        "Propagation delay [s]"
    ].to_numpy()

    positions = x + (index - 0.5) * width

    bars = ax2.bar(
        positions,
        delays,
        width,
        label=chemistry
    )

    ax2.bar_label(
        bars,
        fmt="%.1f",
        padding=3
    )

ax2.set_title("Propagation Delay Between Adjacent Cells")
ax2.set_xlabel("Cell transition")
ax2.set_ylabel("Propagation delay [s]")

ax2.set_xticks(x)
ax2.set_xticklabels([
    "1 → 2",
    "2 → 3",
    "3 → 4",
    "4 → 5"
])

ax2.grid(axis="y", alpha=0.3)
ax2.legend()

fig2.tight_layout()

fig2.savefig(
    results_dir / "tr_propagation_delays.png",
    dpi=300
)

print("tr_propagation_timeline.png")
print("tr_propagation_delays.png")

# Open both figures
plt.show()

print("\nPropagation analysis completed.")
