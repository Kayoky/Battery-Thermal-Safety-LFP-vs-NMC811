
from pathlib import Path
import pandas as pd
import numpy as np

project_dir = Path(__file__).resolve().parent.parent
data_dir = project_dir / "data"
output_dir = project_dir / "data_cleaned"

for chemistry in ["LFP", "NMC811"]:

    input_folder = data_dir / chemistry
    output_folder = output_dir / chemistry
    output_folder.mkdir(parents=True, exist_ok=True)

    for file in sorted(input_folder.glob("T[1-6]_*.txt")):

        df = pd.read_csv(file)

        temperature = df["Temperature [°C]"]

        # Flag physically implausible negative temperatures
        invalid_mask = temperature < 0
        invalid_count = invalid_mask.sum()

        # Replace invalid values with NaN
        df.loc[invalid_mask, "Temperature [°C]"] = np.nan

        # Save cleaned data
        output_file = output_folder / file.name
        df.to_csv(output_file, index=False)

        print(
            f"{chemistry} | {file.name} | "
            f"Invalid values: {invalid_count}"
        )
