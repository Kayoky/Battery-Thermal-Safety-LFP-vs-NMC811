from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# STEP 20.3 - THERMAL FRONT PROPAGATION ANALYSIS
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT / "data" / "surface"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]

THRESHOLDS_C = [100.0, 150.0, 200.0]

# Require the threshold to remain reached for at least
# this duration to suppress isolated measurement spikes.
MIN_PERSISTENCE_S = 0.5

# Previously estimated TR events from Step 11.
# These are algorithmic estimates, not definitive TR onset.
TR_EVENTS = {
    "LFP": [183.9, 649.1, 1154.0, 1646.3, 2022.0],
    "NMC811": [156.9, 211.0, 287.5, 379.4, 551.5],
}

CHANNEL_PATTERN = re.compile(
    r"^Temperature\s+(\d+)\s+\[°C\]$"
)

all_crossing_rows = []
all_summary_rows = []
all_neighbor_rows = []


# ============================================================
# THRESHOLD DETECTION
# ============================================================

def find_persistent_crossing(
    time,
    temperature,
    threshold,
    persistence_s
):
    """
    Return the first time at which temperature reaches
    threshold and remains at or above it for at least
    persistence_s.

    The first crossing is linearly interpolated between
    neighboring samples.

    Returns NaN if no persistent crossing is found.
    """

    valid = np.isfinite(temperature)

    above = (
        valid
        & (temperature >= threshold)
    )

    if not np.any(above):
        return np.nan

    # Find contiguous True intervals.
    padded = np.concatenate(
        ([False], above, [False])
    )

    changes = np.diff(
        padded.astype(int)
    )

    starts = np.where(changes == 1)[0]
    ends = np.where(changes == -1)[0] - 1

    for start, end in zip(starts, ends):

        duration = time[end] - time[start]

        if duration < persistence_s:
            continue

        # First point already above threshold.
        if start == 0:
            return float(time[start])

        t0 = time[start - 1]
        t1 = time[start]

        T0 = temperature[start - 1]
        T1 = temperature[start]

        # Interpolate if both neighboring samples
        # are finite and form a valid crossing.
        if (
            np.isfinite(T0)
            and np.isfinite(T1)
            and T0 < threshold <= T1
            and T1 != T0
        ):
            fraction = (
                (threshold - T0)
                / (T1 - T0)
            )

            return float(
                t0 + fraction * (t1 - t0)
            )

        return float(t1)

    return np.nan


# ============================================================
# PROCESS EACH CHEMISTRY
# ============================================================

for chemistry in CHEMISTRIES:

    file = (
        DATA_DIR
        / chemistry
        / f"T_surf_{chemistry}.txt"
    )

    if not file.exists():
        raise FileNotFoundError(file)

    df = pd.read_csv(
        file,
        encoding="utf-8-sig"
    )

    if "Time [s]" not in df.columns:
        raise ValueError(
            f"Missing time column in {file.name}"
        )

    channel_info = []

    for column in df.columns:

        match = CHANNEL_PATTERN.match(
            str(column)
        )

        if match:
            channel_info.append(
                (int(match.group(1)), column)
            )

    channel_info.sort(
        key=lambda item: item[0]
    )

    channel_numbers = np.array(
        [item[0] for item in channel_info],
        dtype=int
    )

    channel_columns = [
        item[1] for item in channel_info
    ]

    n_channels = len(channel_columns)

    if n_channels < 2:
        raise ValueError(
            f"Not enough channels in {file.name}"
        )

    expected = np.arange(
        1,
        n_channels + 1
    )

    if not np.array_equal(
        channel_numbers,
        expected
    ):
        raise ValueError(
            f"Unexpected channel numbering: {file.name}"
        )

    time = pd.to_numeric(
        df["Time [s]"],
        errors="coerce"
    ).to_numpy(dtype=float)

    temperature = df[
        channel_columns
    ].apply(
        pd.to_numeric,
        errors="coerce"
    ).to_numpy(dtype=float)

    # Analyze experimental time only.
    active_mask = (
        np.isfinite(time)
        & (time >= 0)
    )

    time = time[active_mask]

    temperature = temperature[
        active_mask, :
    ]

    if len(time) < 2:
        raise ValueError(
            f"Insufficient time samples: {file.name}"
        )

    if np.any(np.diff(time) <= 0):
        raise ValueError(
            f"Time is not strictly increasing: {file.name}"
        )

    # Normalized channel-index coordinate.
    # This is NOT a verified physical coordinate.
    x_normalized = (
        channel_numbers - 1
    ) / (n_channels - 1)

    crossing_matrix = np.full(
        (
            n_channels,
            len(THRESHOLDS_C)
        ),
        np.nan
    )

    print("\n" + "=" * 75)
    print(
        f"THERMAL FRONT PROPAGATION - {chemistry}"
    )
    print("=" * 75)

    print(f"Channels: {n_channels}")
    print(f"Time samples: {len(time):,}")
    print(
        f"Time range: "
        f"{time[0]:.2f} to "
        f"{time[-1]:.2f} s"
    )

    # ========================================================
    # COMPUTE THRESHOLD CROSSINGS
    # ========================================================

    for threshold_index, threshold in enumerate(
        THRESHOLDS_C
    ):

        for j, channel_number in enumerate(
            channel_numbers
        ):

            crossing = find_persistent_crossing(
                time=time,
                temperature=temperature[:, j],
                threshold=threshold,
                persistence_s=MIN_PERSISTENCE_S
            )

            crossing_matrix[
                j,
                threshold_index
            ] = crossing

            all_crossing_rows.append({
                "Chemistry": chemistry,
                "Channel number": int(
                    channel_number
                ),
                "Normalized channel coordinate":
                    float(x_normalized[j]),
                "Threshold [C]": threshold,
                "First persistent crossing [s]":
                    crossing,
                "Physical x [mm]": np.nan,
            })

        crossings = crossing_matrix[
            :,
            threshold_index
        ]

        valid = np.isfinite(
            crossings
        )

        reached_count = int(
            valid.sum()
        )

        if reached_count > 0:

            first_crossing = float(
                np.min(crossings[valid])
            )

            last_crossing = float(
                np.max(crossings[valid])
            )

            median_crossing = float(
                np.median(crossings[valid])
            )

            mean_crossing = float(
                np.mean(crossings[valid])
            )

            crossing_span = (
                last_crossing
                - first_crossing
            )

        else:

            first_crossing = np.nan
            last_crossing = np.nan
            median_crossing = np.nan
            mean_crossing = np.nan
            crossing_span = np.nan

        # Neighboring channel delays.
        neighbor_delays = []

        for j in range(
            n_channels - 1
        ):

            t_a = crossings[j]
            t_b = crossings[j + 1]

            if (
                np.isfinite(t_a)
                and np.isfinite(t_b)
            ):

                delay = float(
                    t_b - t_a
                )

                neighbor_delays.append(
                    delay
                )

                all_neighbor_rows.append({
                    "Chemistry": chemistry,
                    "Threshold [C]": threshold,
                    "Channel A": int(
                        channel_numbers[j]
                    ),
                    "Channel B": int(
                        channel_numbers[j + 1]
                    ),
                    "Delta crossing time [s]":
                        delay,
                })

        # Diagnostic only:
        # Neighboring delays may be negative because
        # a higher-index channel can heat earlier.
        median_neighbor_delay = (
            float(
                np.median(
                    neighbor_delays
                )
            )
            if neighbor_delays
            else np.nan
        )

        all_summary_rows.append({
            "Chemistry": chemistry,
            "Threshold [C]": threshold,
            "Total channels": n_channels,
            "Channels reaching threshold":
                reached_count,
            "First crossing [s]":
                first_crossing,
            "Last crossing [s]":
                last_crossing,
            "Mean crossing [s]":
                mean_crossing,
            "Median crossing [s]":
                median_crossing,
            "Crossing time span [s]":
                crossing_span,
            "Median neighboring delay [s]":
                median_neighbor_delay,
        })

        print(
            f"\nThreshold: {threshold:.0f} °C"
        )

        print(
            f"Channels reaching threshold: "
            f"{reached_count}/{n_channels}"
        )

        print(
            f"First crossing: "
            f"{first_crossing:.2f} s"
        )

        print(
            f"Last crossing: "
            f"{last_crossing:.2f} s"
        )

        print(
            f"Median crossing: "
            f"{median_crossing:.2f} s"
        )

        print(
            f"Crossing time span: "
            f"{crossing_span:.2f} s"
        )

        print(
            "Median neighboring delay: "
            f"{median_neighbor_delay:.2f} s"
        )

    # ========================================================
    # FIGURE 1 - THRESHOLD CROSSING TIMES
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(13, 6)
    )

    for threshold_index, threshold in enumerate(
        THRESHOLDS_C
    ):

        ax.plot(
            channel_numbers,
            crossing_matrix[
                :,
                threshold_index
            ],
            marker=".",
            markersize=4,
            linewidth=1.2,
            label=f"{threshold:.0f}°C"
        )

    ax.set_title(
        f"{chemistry} - Temperature Threshold Crossing Times"
    )

    ax.set_xlabel(
        "Temperature Channel Number"
    )

    ax.set_ylabel(
        "First Persistent Crossing Time [s]"
    )

    ax.grid(
        True,
        alpha=0.3
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_front_crossings_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # FIGURE 2 - NEIGHBORING CHANNEL DELAYS
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(13, 6)
    )

    for threshold_index, threshold in enumerate(
        THRESHOLDS_C
    ):

        crossings = crossing_matrix[
            :,
            threshold_index
        ]

        delays = np.diff(
            crossings
        )

        ax.plot(
            channel_numbers[:-1],
            delays,
            marker=".",
            markersize=4,
            linewidth=1,
            label=f"{threshold:.0f}°C"
        )

    ax.axhline(
        0,
        linestyle="--",
        linewidth=1
    )

    ax.set_title(
        f"{chemistry} - Neighboring Channel Crossing Delays"
    )

    ax.set_xlabel(
        "Starting Channel Number"
    )

    ax.set_ylabel(
        "Delta Crossing Time [s]"
    )

    ax.grid(
        True,
        alpha=0.3
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_front_neighbor_delays_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # FIGURE 3 - CROSSING TIMES WITH TR EVENT MARKERS
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(13, 6)
    )

    threshold_index = THRESHOLDS_C.index(
        150.0
    )

    ax.plot(
        channel_numbers,
        crossing_matrix[
            :,
            threshold_index
        ],
        marker=".",
        markersize=4,
        linewidth=1.2,
        label="150°C crossing"
    )

    for cell_index, event_time in enumerate(
        TR_EVENTS[chemistry],
        start=1
    ):

        ax.axhline(
            event_time,
            linestyle="--",
            linewidth=1,
            alpha=0.7,
            label=(
                f"Estimated TR Cell {cell_index}"
            )
        )

    ax.set_title(
        f"{chemistry} - Surface Crossings and Estimated TR Events"
    )

    ax.set_xlabel(
        "Temperature Channel Number"
    )

    ax.set_ylabel(
        "Time [s]"
    )

    ax.grid(
        True,
        alpha=0.3
    )

    ax.legend(
        fontsize=8,
        ncol=2
    )

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_front_vs_TR_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)


# ============================================================
# SAVE RESULTS
# ============================================================

crossing_df = pd.DataFrame(
    all_crossing_rows
)

summary_df = pd.DataFrame(
    all_summary_rows
)

neighbor_df = pd.DataFrame(
    all_neighbor_rows
)

crossing_df.to_csv(
    RESULTS_DIR
    / "surface_front_crossing_times.csv",
    index=False
)

summary_df.to_csv(
    RESULTS_DIR
    / "surface_front_summary.csv",
    index=False
)

neighbor_df.to_csv(
    RESULTS_DIR
    / "surface_front_neighbor_delays.csv",
    index=False
)

print("\n" + "=" * 75)
print("FINAL SUMMARY")
print("=" * 75)

print(
    summary_df.round(2).to_string(
        index=False
    )
)

print("\nSAVED CSV FILES:")
print("surface_front_crossing_times.csv")
print("surface_front_summary.csv")
print("surface_front_neighbor_delays.csv")

print("\nSAVED FIGURES:")

for chemistry in CHEMISTRIES:

    print(
        f"surface_front_crossings_{chemistry}.png"
    )

    print(
        f"surface_front_neighbor_delays_{chemistry}.png"
    )

    print(
        f"surface_front_vs_TR_{chemistry}.png"
    )

print("\nStep 20.3 completed.")