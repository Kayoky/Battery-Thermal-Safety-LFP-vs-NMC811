
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# STEP 19 - SURFACE TEMPERATURE DIAGNOSTICS
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "surface"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]

# Diagnostic threshold only.
# Exceeding it does NOT automatically mean invalid data.
RATE_THRESHOLD = 100.0  # °C/s

# Do not calculate rates across unusually large gaps.
MAX_GAP_FACTOR = 3.0

# Number of hottest channels shown individually.
TOP_N = 10

all_rows = []
chemistry_rows = []

# ============================================================
# MAIN ANALYSIS
# ============================================================

for chemistry in CHEMISTRIES:

    file = (
        DATA_DIR
        / chemistry
        / f"T_surf_{chemistry}.txt"
    )

    if not file.exists():
        raise FileNotFoundError(file)

    df = pd.read_csv(file, encoding="utf-8-sig")

    if "Time [s]" not in df.columns:
        raise ValueError(
            f"Time column not found in {file.name}"
        )

    temp_columns = [
        c for c in df.columns
        if str(c).startswith("Temperature ")
    ]

    if not temp_columns:
        raise ValueError(
            f"No temperature columns in {file.name}"
        )

    time = pd.to_numeric(
        df["Time [s]"],
        errors="coerce"
    ).to_numpy(dtype=float)

    temperatures = df[temp_columns].apply(
        pd.to_numeric,
        errors="coerce"
    )

    if not np.isfinite(time).all():
        raise ValueError(
            f"Invalid timestamps in {file.name}"
        )

    dt = np.diff(time)

    if np.any(dt <= 0):
        raise ValueError(
            f"Non-increasing timestamps in {file.name}"
        )

    median_dt = float(np.median(dt))
    max_gap = MAX_GAP_FACTOR * median_dt

    print("\n" + "=" * 75)
    print(f"SURFACE DIAGNOSTICS - {chemistry}")
    print("=" * 75)

    print(f"File: {file.name}")
    print(f"Samples: {len(df):,}")
    print(f"Temperature channels: {len(temp_columns)}")
    print(f"Time range: {time[0]:.2f} to {time[-1]:.2f} s")
    print(f"Median dt: {median_dt:.5f} s")

    # --------------------------------------------------------
    # CHANNEL-BY-CHANNEL ANALYSIS
    # --------------------------------------------------------

    for col in temp_columns:

        T = temperatures[col].to_numpy(
            dtype=float
        )

        finite = np.isfinite(T)
        valid_indices = np.flatnonzero(finite)

        invalid_count = int((~finite).sum())
        negative_count = int(
            np.sum(T[finite] < 0)
        )

        if len(valid_indices) == 0:
            row = {
                "Chemistry": chemistry,
                "Channel": col,
                "Valid samples": 0,
                "Invalid samples": invalid_count,
                "Negative samples": 0,
                "Tmin [C]": np.nan,
                "Tmax [C]": np.nan,
                "Time at Tmax [s]": np.nan,
                "Max abs rate [C/s]": np.nan,
                "Rate flags >100C/s": 0,
            }

            all_rows.append(row)
            continue

        Tmin = float(np.min(T[finite]))
        Tmax = float(np.max(T[finite]))

        idx_max = valid_indices[
            np.argmax(T[finite])
        ]

        t_max = float(time[idx_max])

        # Adjacent finite samples only.
        # Never bridge NaN values.
        pair_valid = finite[:-1] & finite[1:]

        pair_valid &= dt <= max_gap

        dT = np.diff(T)

        rates = np.full(
            len(dt),
            np.nan,
            dtype=float
        )

        rates[pair_valid] = (
            dT[pair_valid] / dt[pair_valid]
        )

        finite_rates = rates[np.isfinite(rates)]

        if len(finite_rates):
            max_abs_rate = float(
                np.max(np.abs(finite_rates))
            )

            rate_flags = int(
                np.sum(
                    np.abs(finite_rates)
                    > RATE_THRESHOLD
                )
            )
        else:
            max_abs_rate = np.nan
            rate_flags = 0

        row = {
            "Chemistry": chemistry,
            "Channel": col,
            "Valid samples": int(finite.sum()),
            "Invalid samples": invalid_count,
            "Negative samples": negative_count,
            "Tmin [C]": Tmin,
            "Tmax [C]": Tmax,
            "Time at Tmax [s]": t_max,
            "Max abs rate [C/s]": max_abs_rate,
            "Rate flags >100C/s": rate_flags,
        }

        all_rows.append(row)

    # --------------------------------------------------------
    # CHEMISTRY SUMMARY
    # --------------------------------------------------------

    chemistry_data = pd.DataFrame(
        [
            r for r in all_rows
            if r["Chemistry"] == chemistry
        ]
    )

    max_values = chemistry_data["Tmax [C]"]

    hottest = chemistry_data.loc[
        max_values.idxmax()
    ]

    summary = {
        "Chemistry": chemistry,
        "Channels": len(temp_columns),
        "Samples per channel": len(df),
        "Total invalid samples": int(
            chemistry_data["Invalid samples"].sum()
        ),
        "Total negative samples": int(
            chemistry_data["Negative samples"].sum()
        ),
        "Total rate flags": int(
            chemistry_data["Rate flags >100C/s"].sum()
        ),
        "Global Tmax [C]": float(
            max_values.max()
        ),
        "Hottest channel": hottest["Channel"],
        "Median channel Tmax [C]": float(
            max_values.median()
        ),
    }

    chemistry_rows.append(summary)

    print("\nDATA QUALITY:")
    print(
        f"Invalid samples: "
        f"{summary['Total invalid samples']}"
    )
    print(
        f"Negative samples: "
        f"{summary['Total negative samples']}"
    )
    print(
        f"Rate flags >100 C/s: "
        f"{summary['Total rate flags']}"
    )

    print("\nTEMPERATURE EXTREMES:")
    print(
        f"Global Tmax: "
        f"{summary['Global Tmax [C]']:.2f} C"
    )
    print(
        f"Hottest channel: "
        f"{summary['Hottest channel']}"
    )
    print(
        f"Median channel Tmax: "
        f"{summary['Median channel Tmax [C]']:.2f} C"
    )

    top = chemistry_data.nlargest(
        TOP_N,
        "Tmax [C]"
    )

    print(f"\nTOP {TOP_N} HOTTEST CHANNELS:")
    print(
        top[
            [
                "Channel",
                "Tmax [C]",
                "Time at Tmax [s]",
                "Rate flags >100C/s",
            ]
        ].round(2).to_string(index=False)
    )

    # --------------------------------------------------------
    # FIGURE 1 - ALL TEMPERATURE CHANNELS
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(14, 7))

    # Plot every channel.
    # Downsample only for visualization.
    plot_step = max(1, len(df) // 10000)

    plot_time = time[::plot_step]

    for col in temp_columns:
        ax.plot(
            plot_time,
            temperatures[col].to_numpy()[::plot_step],
            linewidth=0.55,
            alpha=0.35
        )

    ax.set_title(
        f"{chemistry} - All Surface Temperature Channels"
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [°C]")
    ax.grid(True, alpha=0.25)

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_all_channels_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    # --------------------------------------------------------
    # FIGURE 2 - CHANNEL MAXIMUM TEMPERATURES
    # --------------------------------------------------------

    channel_numbers = [
        int(col.split("Temperature ")[1].split(" ")[0])
        for col in temp_columns
    ]

    fig, ax = plt.subplots(figsize=(14, 6))

    ax.plot(
        channel_numbers,
        chemistry_data["Tmax [C]"],
        marker=".",
        markersize=5,
        linewidth=1.2
    )

    ax.set_title(
        f"{chemistry} - Maximum Temperature by Channel"
    )
    ax.set_xlabel("Temperature channel number")
    ax.set_ylabel("Maximum temperature [°C]")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_channel_maxima_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    # --------------------------------------------------------
    # FIGURE 3 - TOP 10 HOTTEST CHANNELS
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(14, 7))

    for col in top["Channel"]:

        ax.plot(
            plot_time,
            temperatures[col].to_numpy()[::plot_step],
            linewidth=1.2,
            label=col
        )

    ax.set_title(
        f"{chemistry} - Top {TOP_N} Hottest Channels"
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [°C]")
    ax.grid(True, alpha=0.3)
    ax.legend(
        fontsize=8,
        ncol=2
    )

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_top_channels_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)


# ============================================================
# SAVE CSV RESULTS
# ============================================================

channel_summary = pd.DataFrame(all_rows)
chemistry_summary = pd.DataFrame(chemistry_rows)

channel_summary.to_csv(
    RESULTS_DIR / "surface_channel_diagnostics.csv",
    index=False
)

chemistry_summary.to_csv(
    RESULTS_DIR / "surface_chemistry_diagnostics.csv",
    index=False
)

print("\n" + "=" * 75)
print("FINAL COMPARISON")
print("=" * 75)

print(
    chemistry_summary.round(2).to_string(index=False)
)

print("\nSAVED FILES:")
print("surface_channel_diagnostics.csv")
print("surface_chemistry_diagnostics.csv")

for chemistry in CHEMISTRIES:
    print(f"surface_all_channels_{chemistry}.png")
    print(f"surface_channel_maxima_{chemistry}.png")
    print(f"surface_top_channels_{chemistry}.png")

print("\nStep 19 completed.")
