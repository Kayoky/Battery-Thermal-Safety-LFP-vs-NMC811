
from pathlib import Path
import pandas as pd

# Main project directory
project_dir = Path(__file__).resolve().parent.parent

# Experimental data directories
data_dir = project_dir / "data"

for chemistry in ["LFP", "NMC811"]:
    folder = data_dir / chemistry

    print(f"\n{'=' * 50}")
    print(f"Battery chemistry: {chemistry}")
    print(f"{'=' * 50}")

    files = sorted(folder.glob("*.txt"))

    if not files:
        print("No TXT files found.")
        continue

    print(f"Number of files: {len(files)}")

    for file in files:
        print(f"\nFile: {file.name}")

        # Display first lines to identify the file format
        with open(file, "r", encoding="utf-8-sig") as f:
            for _ in range(5):
                line = f.readline()
                if not line:
                    break
                print(line.rstrip())

# Inspect suspicious NMC811 T5 measurements
file = data_dir / "NMC811" / "T5_NMC811.txt"
df = pd.read_csv(file)

invalid = df[df["Temperature [°C]"] < 0]

print("\nNMC811 - T5 sensor inspection")
print(f"Total measurements: {len(df)}")
print(f"Negative measurements: {len(invalid)}")

if not invalid.empty:
    print(f"Minimum temperature: {invalid['Temperature [°C]'].min():.2f} °C")
    print(f"First negative value at: {invalid['Time [s]'].iloc[0]:.2f} s")
    print(f"Last negative value at: {invalid['Time [s]'].iloc[-1]:.2f} s")

    print("\nFirst 10 negative measurements:")
    print(invalid.head(10).to_string(index=False))
