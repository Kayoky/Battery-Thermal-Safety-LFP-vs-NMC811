
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 15 - VENTING TEMPERATURE VISUALIZATION
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "venting"
RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]

for chemistry in CHEMISTRIES:

    fig, ax = plt.subplots(figsize=(13, 7))

    print(f"\nProcessing {chemistry}...")

    for cell in range(1, 6):

        file = (
            DATA_DIR
            / chemistry
            / f"V{cell}_{chemistry}.txt"
        )

        df = pd.read_csv(file)

        time = df["Time [s]"]
        temp = df["Temperature [°C]"]

        ax.plot(
            time,
            temp,
            linewidth=1.2,
            label=f"Cell {cell}"
        )

        idx = temp.idxmax()
        t_max = time.loc[idx]
        T_max = temp.loc[idx]

        ax.scatter(
            t_max,
            T_max,
            s=35,
            zorder=5
        )

        print(
            f"Cell {cell}: "
            f"Tmax = {T_max:.2f} °C "
            f"at t = {t_max:.2f} s"
        )

    ax.axvline(
        x=0,
        color="black",
        linestyle="--",
        linewidth=1,
        label="Heater ON"
    )

    ax.set_title(
        f"Raw Safety Vent Temperature Measurements - {chemistry}"
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Measured temperature [°C]")

    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()

    output = (
        RESULTS_DIR
        / f"venting_temperatures_{chemistry}_raw.png"
    )

    fig.savefig(output, dpi=300)
    plt.show()
    plt.close(fig)

    print(f"Saved: {output.name}")

print("\nVenting temperature visualization completed.")
