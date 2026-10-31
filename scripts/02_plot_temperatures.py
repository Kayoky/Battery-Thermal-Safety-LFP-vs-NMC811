
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# Project directories
project_dir = Path(__file__).resolve().parent.parent
results_dir = project_dir / "results"
results_dir.mkdir(exist_ok=True)

# Define the exact order of the plots
plot_config = [
    ("LFP", "raw", project_dir / "data" / "LFP"),
    ("LFP", "cleaned", project_dir / "data_cleaned" / "LFP"),
    ("NMC811", "raw", project_dir / "data" / "NMC811"),
    ("NMC811", "cleaned", project_dir / "data_cleaned" / "NMC811")
]

for chemistry, version, folder in plot_config:

    fig, ax = plt.subplots(figsize=(12, 6))

    files = sorted(folder.glob("T[1-6]_*.txt"))

    if not files:
        print(f"No files found in: {folder}")
        plt.close(fig)
        continue

    for file in files:

        df = pd.read_csv(file)

        time = df["Time [s]"]
        temperature = df["Temperature [°C]"]

        sensor = file.stem.split("_")[0]

        ax.plot(
            time,
            temperature,
            label=sensor,
            linewidth=1.2
        )

    ax.set_title(
        f"Cell Temperature Measurements - {chemistry} ({version})"
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [°C]")

    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()

    # Save each plot
    output_file = (
        results_dir / f"temperatures_{chemistry}_{version}.png"
    )

    fig.savefig(output_file, dpi=300)

    print(f"Saved: {output_file.name}")

# Display all four figures
plt.show()
