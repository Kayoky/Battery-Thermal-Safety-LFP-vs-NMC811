
from pathlib import Path
import pandas as pd

# Project directories
project_dir = Path(__file__).resolve().parent.parent
voltage_dir = project_dir / "data" / "voltage"

for chemistry in ["LFP", "NMC811"]:

    print(f"\n{'=' * 60}")
    print(f"Battery chemistry: {chemistry}")
    print(f"{'=' * 60}")

    folder = voltage_dir / chemistry

    if not folder.exists():
        print(f"Folder not found: {folder}")
        continue

    files = sorted(folder.glob("*.txt"))

    if not files:
        print("No .txt files found.")
        continue

    for file in files:

        print(f"\nFile: {file.name}")

        try:
            df = pd.read_csv(file)

            print(f"Number of measurements: {len(df)}")
            print(f"Columns: {list(df.columns)}")

            print("\nFirst 5 measurements:")
            print(df.head().to_string(index=False))

            print("\nData quality:")
            print(f"Missing values: {df.isna().sum().sum()}")
            print(f"Duplicate rows: {df.duplicated().sum()}")

            # Check time column if available
            time_columns = [
                col for col in df.columns
                if "time" in col.lower()
            ]

            if time_columns:
                time = pd.to_numeric(
                    df[time_columns[0]],
                    errors="coerce"
                )

                print(f"Start time: {time.min():.2f} s")
                print(f"End time: {time.max():.2f} s")

                dt = time.diff()

                print(
                    f"Median sampling interval: "
                    f"{dt.median():.4f} s"
                )

                print(
                    f"Non-increasing time steps: "
                    f"{(dt.dropna() <= 0).sum()}"
                )

        except Exception as error:
            print(f"Error reading file: {error}")

print("\nVoltage data inspection completed.")
