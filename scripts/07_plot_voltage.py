
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# Project directories
project_dir = Path(__file__).resolve().parent.parent
voltage_dir = project_dir / "data" / "voltage"
results_dir = project_dir / "results"

results_dir.mkdir(exist_ok=True)

# Plot order
chemistries = ["LFP", "NMC811"]

for chemistry in chemistries:

    folder = voltage_dir / chemistry

    fig, ax = plt.subplots(figsize=(12, 6))

    for cell in range(1, 6):

        file = folder / f"U{cell}_{chemistry}.txt"

        df = pd.read_csv(file)

        time = df["Time [s]"]
        voltage = df["Voltage [V]"]

        ax.plot(
            time,
            voltage,
            label=f"Cell {cell}",
            linewidth=1.2
        )

        print(
            f"{chemistry} | Cell {cell} | "
            f"U_min = {voltage.min():.3f} V | "
            f"U_max = {voltage.max():.3f} V"
        )

    ax.set_title(
        f"Cell Voltage Measurements - {chemistry}"
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Voltage [V]")

    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    # Save figure
    output_file = (
        results_dir / f"voltage_{chemistry}_raw.png"
    )

    fig.savefig(output_file, dpi=300)

    print(f"Saved: {output_file.name}\n")

# Display both figures
plt.show()
