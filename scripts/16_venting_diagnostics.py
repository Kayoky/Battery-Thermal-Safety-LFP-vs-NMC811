
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 16 - VENTING TEMPERATURE DIAGNOSTICS
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "venting"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]
CELLS = range(1, 6)

# Diagnostic settings
SMOOTHING_WINDOW_S = 1.0
PEAK_WINDOW_S = 30.0
RATE_THRESHOLD = 500.0  # °C/s, diagnostic flag only

# No automatic removal of flagged measurements.
summary = []

# ============================================================
# LOAD PREVIOUS TR EVENT ESTIMATES
# ============================================================

tr_file = RESULTS_DIR / "refined_tr_detection_summary.csv"

if tr_file.exists():
    tr_df = pd.read_csv(tr_file)
    print("Loaded refined TR event estimates.")
else:
    tr_df = pd.DataFrame()
    print("WARNING: TR summary not found.")

def get_tr_time(chemistry, cell):
    if tr_df.empty:
        return np.nan

    required = {
        "Chemistry",
        "Cell",
        "Matched temperature event [s]"
    }

    if not required.issubset(tr_df.columns):
        return np.nan

    selection = tr_df[
        (tr_df["Chemistry"] == chemistry)
        & (tr_df["Cell"] == cell)
    ]

    if selection.empty:
        return np.nan

    return float(
        selection.iloc[0]["Matched temperature event [s]"]
    )

# ============================================================
# MAIN ANALYSIS
# ============================================================

for chemistry in CHEMISTRIES:

    print("\n" + "=" * 70)
    print(f"VENTING DIAGNOSTICS - {chemistry}")
    print("=" * 70)

    for cell in CELLS:

        file = (
            DATA_DIR
            / chemistry
            / f"V{cell}_{chemistry}.txt"
        )

        if not file.exists():
            print(f"Missing: {file.name}")
            continue

        df = pd.read_csv(file)

        time = pd.to_numeric(
            df["Time [s]"], errors="coerce"
        ).to_numpy(dtype=float)

        temp = pd.to_numeric(
            df["Temperature [°C]"], errors="coerce"
        ).to_numpy(dtype=float)

        valid = (
            np.isfinite(time)
            & np.isfinite(temp)
            & (time >= 0)
        )

        time = time[valid]
        temp = temp[valid]

        if len(time) < 3:
            print(f"Insufficient data: {file.name}")
            continue

        if np.any(np.diff(time) <= 0):
            raise ValueError(
                f"{file.name}: non-increasing time."
            )

        # ----------------------------------------------------
        # SMOOTHING
        # ----------------------------------------------------

        dt = float(np.median(np.diff(time)))

        window = max(
            3,
            int(round(SMOOTHING_WINDOW_S / dt))
        )

        smooth = (
            pd.Series(temp)
            .rolling(
                window=window,
                center=True,
                min_periods=1
            )
            .mean()
            .to_numpy()
        )

        # ----------------------------------------------------
        # TEMPERATURE DERIVATIVES
        # ----------------------------------------------------

        raw_rate = np.gradient(temp, time)
        smooth_rate = np.gradient(smooth, time)

        max_heating_rate = float(np.max(smooth_rate))
        max_cooling_rate = float(np.min(smooth_rate))

        flagged = np.abs(raw_rate) > RATE_THRESHOLD
        flagged_count = int(np.sum(flagged))

        # ----------------------------------------------------
        # MAXIMUM TEMPERATURE
        # ----------------------------------------------------

        idx_max = int(np.argmax(temp))

        Tmax = float(temp[idx_max])
        t_Tmax = float(time[idx_max])

        # Peak duration above 90% of Tmax.
        # Find the contiguous interval around the raw maximum.
        peak_threshold = 0.9 * Tmax

        left = idx_max
        right = idx_max

        while (
            left > 0
            and temp[left - 1] >= peak_threshold
        ):
            left -= 1

        while (
            right < len(temp) - 1
            and temp[right + 1] >= peak_threshold
        ):
            right += 1

        peak_duration = float(
            time[right] - time[left]
        )

        # ----------------------------------------------------
        # TR ALIGNMENT
        # ----------------------------------------------------

        tr_time = get_tr_time(chemistry, cell)

        if np.isfinite(tr_time):
            peak_delay = t_Tmax - tr_time
        else:
            peak_delay = np.nan

        # ----------------------------------------------------
        # SUMMARY
        # ----------------------------------------------------

        summary.append({
            "Chemistry": chemistry,
            "Cell": cell,
            "Tmax [C]": Tmax,
            "Time at Tmax [s]": t_Tmax,
            "Max smoothed heating rate [C/s]":
                max_heating_rate,
            "Max smoothed cooling rate [C/s]":
                max_cooling_rate,
            "Samples with abs raw dTdt > 500 C/s":
                flagged_count,
            "Peak duration above 90pct Tmax [s]":
                peak_duration,
            "Estimated TR event [s]": tr_time,
            "Peak minus TR event [s]": peak_delay
        })

        # ----------------------------------------------------
        # FIGURE: FULL SIGNAL + LOCAL PEAK
        # ----------------------------------------------------

        fig, axes = plt.subplots(
            3, 1,
            figsize=(12, 10),
            gridspec_kw={
                "height_ratios": [1.1, 1.1, 0.9]
            }
        )

        # Full temperature history
        axes[0].plot(
            time,
            temp,
            linewidth=1,
            label="Raw vent temperature"
        )

        axes[0].plot(
            time,
            smooth,
            linewidth=1.5,
            alpha=0.8,
            label="1 s smoothed"
        )

        axes[0].axvline(
            t_Tmax,
            color="red",
            linestyle="--",
            label="Maximum temperature"
        )

        if np.isfinite(tr_time):
            axes[0].axvline(
                tr_time,
                color="green",
                linestyle=":",
                label="Estimated TR event"
            )

        axes[0].set_title(
            f"{chemistry} - Cell {cell}: Vent Temperature"
        )
        axes[0].set_ylabel("Temperature [°C]")
        axes[0].grid(True, alpha=0.3)
        axes[0].legend()

        # Local peak window
        local = (
            (time >= t_Tmax - PEAK_WINDOW_S)
            & (time <= t_Tmax + PEAK_WINDOW_S)
        )

        axes[1].plot(
            time[local],
            temp[local],
            linewidth=1,
            label="Raw"
        )

        axes[1].plot(
            time[local],
            smooth[local],
            linewidth=2,
            label="Smoothed"
        )

        axes[1].axhline(
            peak_threshold,
            color="gray",
            linestyle="--",
            label="90% of raw Tmax"
        )

        axes[1].axvline(
            t_Tmax,
            color="red",
            linestyle=":"
        )

        axes[1].set_title(
            "Local Temperature Peak (±30 s)"
        )
        axes[1].set_ylabel("Temperature [°C]")
        axes[1].grid(True, alpha=0.3)
        axes[1].legend()

        # Heating/cooling rates
        axes[2].plot(
            time[local],
            raw_rate[local],
            alpha=0.4,
            linewidth=0.8,
            label="Raw dT/dt"
        )

        axes[2].plot(
            time[local],
            smooth_rate[local],
            linewidth=1.5,
            label="Smoothed dT/dt"
        )

        axes[2].axhline(
            RATE_THRESHOLD,
            color="red",
            linestyle="--",
            alpha=0.6
        )

        axes[2].axhline(
            -RATE_THRESHOLD,
            color="red",
            linestyle="--",
            alpha=0.6
        )

        axes[2].set_title(
            "Temperature Rate Around Maximum"
        )
        axes[2].set_xlabel("Time [s]")
        axes[2].set_ylabel("dT/dt [°C/s]")
        axes[2].grid(True, alpha=0.3)
        axes[2].legend()

        fig.tight_layout()

        output = (
            RESULTS_DIR
            / f"venting_diagnostics_{chemistry}_cell{cell}.png"
        )

        fig.savefig(output, dpi=300)
        plt.close(fig)

        print(
            f"Cell {cell}: "
            f"Tmax={Tmax:.1f} °C, "
            f"t={t_Tmax:.1f} s, "
            f"max heating={max_heating_rate:.1f} °C/s, "
            f"flags={flagged_count}, "
            f"peak duration={peak_duration:.2f} s"
        )

# ============================================================
# SAVE SUMMARY
# ============================================================

summary_df = pd.DataFrame(summary)

if summary_df.empty:
    raise RuntimeError("No venting files were processed.")

summary_df.to_csv(
    RESULTS_DIR / "venting_diagnostics_summary.csv",
    index=False
)

print("\n" + "=" * 70)
print("VENTING DIAGNOSTICS SUMMARY")
print("=" * 70)

print(summary_df.round(2).to_string(index=False))

print("\nSaved:")
print("venting_diagnostics_summary.csv")
print("10 diagnostic PNG figures (5 per chemistry)")
print("\nStep 16 completed.")
