from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


# ============================================================
# STEP 20.2 - SURFACE TEMPERATURE HEATMAPS
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT / "data" / "surface"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]

TAPE_LIMIT_C = 250.0

# Plot resolution
MAX_TIME_POINTS = 3000


# ============================================================
# LOAD AND PROCESS
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

    time = df["Time [s]"].to_numpy(dtype=float)

    channel_columns = [
        col for col in df.columns
        if col.startswith("Temperature ")
    ]

    channel_columns.sort(
        key=lambda x: int(x.split(" ")[1])
    )

    temperature = df[channel_columns].to_numpy(
        dtype=float
    )

    # Only experimental time >= 0
    mask = time >= 0

    time = time[mask]
    temperature = temperature[mask, :]

    # Reduce resolution for plotting
    step = max(
        1,
        len(time) // MAX_TIME_POINTS
    )

    time_plot = time[::step]

    T_plot = temperature[::step, :]

    # Matrix: channels x time
    T_map = T_plot.T

    n_channels = T_map.shape[0]

    print("\n" + "=" * 70)
    print(f"SURFACE HEATMAP ANALYSIS - {chemistry}")
    print("=" * 70)

    print(f"Channels: {n_channels}")
    print(f"Time points used: {len(time_plot)}")
    print(
        f"Temperature range: "
        f"{np.nanmin(T_map):.2f} to "
        f"{np.nanmax(T_map):.2f} °C"
    )

    # ========================================================
    # 1. HEATMAP WITH COLOR SCALE LIMITED TO 250 C
    # ========================================================

    fig, ax = plt.subplots(figsize=(15, 6))

    im = ax.imshow(
        T_map,
        aspect="auto",
        origin="lower",
        extent=[
            time_plot[0],
            time_plot[-1],
            0.5,
            n_channels + 0.5
        ],
        cmap="inferno",
        norm=Normalize(vmin=20, vmax=TAPE_LIMIT_C)
    )

    ax.set_title(
        f"{chemistry} - Surface Temperature Heatmap "
        "(Color Scale Limited to 250°C)"
    )

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature Channel Number")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Temperature [°C]")

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_heatmap_250C_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # 2. FULL TEMPERATURE RANGE HEATMAP
    # ========================================================

    fig, ax = plt.subplots(figsize=(15, 6))

    im = ax.imshow(
        T_map,
        aspect="auto",
        origin="lower",
        extent=[
            time_plot[0],
            time_plot[-1],
            0.5,
            n_channels + 0.5
        ],
        cmap="inferno"
    )

    ax.set_title(
        f"{chemistry} - Full-Range Surface Temperature Heatmap "
        "(>250°C: Increased Measurement Uncertainty)"
    )

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature Channel Number")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Measured Temperature [°C]")

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_heatmap_full_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # 3. MASK OF VALUES ABOVE 250 C
    # ========================================================

    above_limit = T_map > TAPE_LIMIT_C

    fig, ax = plt.subplots(figsize=(15, 6))

    im = ax.imshow(
        above_limit.astype(int),
        aspect="auto",
        origin="lower",
        extent=[
            time_plot[0],
            time_plot[-1],
            0.5,
            n_channels + 0.5
        ],
        cmap="gray_r",
        vmin=0,
        vmax=1
    )

    ax.set_title(
        f"{chemistry} - Surface Measurements Above 250°C"
    )

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature Channel Number")

    cbar = fig.colorbar(
        im,
        ax=ax,
        ticks=[0, 1]
    )

    cbar.ax.set_yticklabels(
        ["≤250°C", ">250°C"]
    )

    fig.tight_layout()

    fig.savefig(
        RESULTS_DIR
        / f"surface_heatmap_above_250C_{chemistry}.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)


print("\nStep 20.2 completed.")