
from pathlib import Path
import pandas as pd
import numpy as np

# ============================================================
# STEP 18 - SURFACE TEMPERATURE DATA INSPECTION
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "surface"
RESULTS_DIR = ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CHEMISTRIES = ["LFP", "NMC811"]

summary = []

print("=" * 75)
print("STEP 18 - SURFACE TEMPERATURE DATA INSPECTION")
print("=" * 75)

for chemistry in CHEMISTRIES:

    folder = DATA_DIR / chemistry

    print(f"\n{'=' * 75}")
    print(f"CHEMISTRY: {chemistry}")
    print(f"{'=' * 75}")

    if not folder.exists():
        print(f"ERROR: Folder not found: {folder}")
        continue

    files = sorted(
        f for f in folder.iterdir()
        if f.is_file() and f.suffix.lower()
        in [".txt", ".csv"]
    )

    print(f"Files found: {len(files)}")

    if not files:
        print("WARNING: No TXT or CSV files found.")
        continue

    for file in files:

        print(f"\n--- {file.name} ---")

        try:
            # Automatically detect common separators.
            df = pd.read_csv(
                file,
                sep=None,
                engine="python",
                encoding="utf-8-sig"
            )

        except Exception as e:
            print(f"READ ERROR: {e}")
            continue

        print(f"Shape: {df.shape}")
        print(f"Columns: {list(df.columns)}")
        print("First 3 rows:")
        print(df.head(3).to_string(index=False))

        numeric = df.apply(
            pd.to_numeric,
            errors="coerce"
        )

        # Try to identify time and temperature columns.
        time_columns = [
            c for c in df.columns
            if "time" in str(c).lower()
            or "zeit" in str(c).lower()
        ]

        temp_columns = [
            c for c in df.columns
            if "temp" in str(c).lower()
            or "°c" in str(c).lower()
        ]

        time_col = (
            time_columns[0]
            if time_columns else None
        )

        temp_col = (
            temp_columns[0]
            if temp_columns else None
        )

        row = {
            "Chemistry": chemistry,
            "File": file.name,
            "Rows": len(df),
            "Columns": len(df.columns),
            "Column names": " | ".join(map(str, df.columns)),
            "Missing values": int(df.isna().sum().sum()),
            "Time column": time_col,
            "Temperature column": temp_col,
        }

        if time_col is not None:

            t = numeric[time_col].to_numpy(
                dtype=float
            )

            valid_t = t[np.isfinite(t)]

            if len(valid_t) > 0:

                row["Start time [s]"] = valid_t[0]
                row["End time [s]"] = valid_t[-1]

            dt = np.diff(t)
            finite_dt = dt[np.isfinite(dt)]

            row["Non-increasing time steps"] = int(
                np.sum(finite_dt <= 0)
            )

            row["Median dt [s]"] = (
                float(np.median(finite_dt))
                if len(finite_dt) else np.nan
            )

        if temp_col is not None:

            T = numeric[temp_col].to_numpy(
                dtype=float
            )

            valid_T = T[np.isfinite(T)]

            row["Invalid temperature values"] = int(
                np.sum(~np.isfinite(T))
            )

            row["Negative temperatures"] = int(
                np.sum(valid_T < 0)
            )

            if len(valid_T) > 0:

                row["Tmin [C]"] = float(
                    np.min(valid_T)
                )

                row["Tmax [C]"] = float(
                    np.max(valid_T)
                )

        summary.append(row)

        print(
            f"Missing values: "
            f"{row['Missing values']}"
        )

        if time_col is not None:
            print(
                f"Time: "
                f"{row.get('Start time [s]', np.nan):.2f}"
                f" to "
                f"{row.get('End time [s]', np.nan):.2f} s"
            )
            print(
                f"Median dt: "
                f"{row.get('Median dt [s]', np.nan):.4f} s"
            )

        if temp_col is not None:
            print(
                f"Temperature range: "
                f"{row.get('Tmin [C]', np.nan):.2f}"
                f" to "
                f"{row.get('Tmax [C]', np.nan):.2f} C"
            )
            print(
                f"Negative temperatures: "
                f"{row['Negative temperatures']}"
            )

# ============================================================
# SAVE SUMMARY
# ============================================================

if summary:

    summary_df = pd.DataFrame(summary)

    output = (
        RESULTS_DIR /
        "surface_data_inspection_summary.csv"
    )

    summary_df.to_csv(
        output,
        index=False
    )

    print("\n" + "=" * 75)
    print("OVERALL SUMMARY")
    print("=" * 75)

    print(
        summary_df.to_string(index=False)
    )

    print(f"\nSaved: {output}")

else:
    print("\nNo files were successfully inspected.")

print("\nStep 18 inspection finished.")
