
from pathlib import Path
import numpy as np
import pandas as pd

# ============================================================
# 14 - VENTING TEMPERATURE DATA INSPECTION
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "venting"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]
summary = []

for chemistry in CHEMISTRIES:

    print("\n" + "=" * 65)
    print(f"VENTING DATA INSPECTION - {chemistry}")
    print("=" * 65)

    for cell in range(1, 6):

        file = (
            DATA_DIR
            / chemistry
            / f"V{cell}_{chemistry}.txt"
        )

        if not file.exists():
            print(f"\nMISSING FILE: {file}")
            continue

        df = pd.read_csv(file)

        print(f"\nFile: {file.name}")
        print(f"Number of measurements: {len(df)}")
        print(f"Columns: {list(df.columns)}")
        print("\nFirst 5 measurements:")
        print(df.head().to_string(index=False))

        if df.shape[1] != 2:
            print("WARNING: Expected two columns.")
            continue

        time = pd.to_numeric(
            df.iloc[:, 0], errors="coerce"
        ).to_numpy(dtype=float)

        temp = pd.to_numeric(
            df.iloc[:, 1], errors="coerce"
        ).to_numpy(dtype=float)

        valid = np.isfinite(time) & np.isfinite(temp)

        missing = int((~valid).sum())
        negative = int(np.sum(temp < 0))
        duplicate_rows = int(df.duplicated().sum())

        if valid.sum() >= 2:

            valid_time = time[valid]
            valid_temp = temp[valid]

            dt = np.diff(valid_time)

            median_dt = float(np.median(dt))
            non_increasing = int(np.sum(dt <= 0))

            t_min = float(np.min(valid_time))
            t_max = float(np.max(valid_time))

            T_min = float(np.min(valid_temp))
            T_max = float(np.max(valid_temp))

        else:
            median_dt = np.nan
            non_increasing = 0
            t_min = np.nan
            t_max = np.nan
            T_min = np.nan
            T_max = np.nan

        print("\nData quality:")
        print(f"Missing/invalid rows: {missing}")
        print(f"Duplicate rows: {duplicate_rows}")
        print(f"Negative temperatures: {negative}")
        print(f"Start time: {t_min:.2f} s")
        print(f"End time: {t_max:.2f} s")
        print(f"Median sampling interval: {median_dt:.4f} s")
        print(f"Non-increasing time steps: {non_increasing}")
        print(f"Minimum temperature: {T_min:.2f} °C")
        print(f"Maximum temperature: {T_max:.2f} °C")

        summary.append({
            "Chemistry": chemistry,
            "Cell": cell,
            "File": file.name,
            "Measurements": len(df),
            "Invalid rows": missing,
            "Duplicate rows": duplicate_rows,
            "Negative temperatures": negative,
            "Start time [s]": t_min,
            "End time [s]": t_max,
            "Median dt [s]": median_dt,
            "Non-increasing time steps": non_increasing,
            "Tmin [C]": T_min,
            "Tmax [C]": T_max
        })

# ============================================================
# SAVE SUMMARY
# ============================================================

summary_df = pd.DataFrame(summary)

output = RESULTS_DIR / "venting_data_inspection.csv"
summary_df.to_csv(output, index=False)

print("\n" + "=" * 65)
print("VENTING DATA INSPECTION SUMMARY")
print("=" * 65)

if not summary_df.empty:
    print(summary_df.round(3).to_string(index=False))

print(f"\nSaved: {output.name}")
print("Venting data inspection completed.")
