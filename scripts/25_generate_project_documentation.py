#!/usr/bin/env python3
"""Step 25 v2: reproducible documentation generator for LFP vs NMC811.
Run from the project root: python scripts/25_generate_project_documentation.py
Requires pandas and numpy. Existing scientific results are never modified.
"""
from __future__ import annotations

import ast
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
import re
import shutil
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
SCRIPTS = ROOT / "scripts"
BACKUPS = ROOT / "documentation_backups"
STAMP = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
SOURCE = ("Schöberl et al., ‘Thermal runaway propagation in automotive lithium-ion "
          "batteries with NMC-811 and LFP cathodes: Safety requirements and "
          "impact on system integration’, eTransportation 19 (2024), 100305.")
ARTICLE_DOI = "10.1016/j.etran.2023.100305"
DATA_DOI = "10.14459/2023mp1717758"
INPUTS = {
    "integrated": "integrated_chemistry_summary.csv",
    "cells": "integrated_cell_summary.csv",
    "criteria": "technology_raw_criteria.csv",
    "normalized": "technology_normalized_criteria.csv",
    "weights": "technology_scenario_weights.csv",
    "scores": "technology_decision_scores.csv",
    "sensitivity": "technology_sensitivity_analysis.csv",
}
CHEMS = ["LFP", "NMC811"]


def load_csv(name: str) -> pd.DataFrame:
    path = RESULTS / INPUTS[name]
    if not path.is_file():
        raise FileNotFoundError(f"Missing required input: {path}")
    return pd.read_csv(path)


def indexed(df: pd.DataFrame, name: str, key: str) -> pd.DataFrame:
    if key not in df.columns:
        raise ValueError(f"{name}: expected column '{key}'; found {list(df.columns)}")
    if df[key].duplicated().any():
        raise ValueError(f"{name}: duplicate {key} values")
    return df.set_index(key)


def number(frame: pd.DataFrame, row: str, col: str) -> float:
    if row not in frame.index or col not in frame.columns:
        raise ValueError(f"Missing required metric: {row} / {col}")
    val = float(frame.loc[row, col])
    if not np.isfinite(val):
        raise ValueError(f"Non-finite metric: {row} / {col}")
    return val


def fmt(x: float, digits: int = 2) -> str:
    return f"{x:.{digits}f}"


def md_table(frame: pd.DataFrame, index_name: str | None = None, digits: int = 2) -> str:
    cols = [index_name or str(frame.index.name or "Metric")] + list(map(str, frame.columns))
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for idx, row in frame.iterrows():
        items = [str(idx)]
        for v in row:
            if pd.isna(v):
                items.append("N/A")
            elif isinstance(v, (int, float, np.integer, np.floating)):
                items.append(fmt(float(v), digits))
            else:
                items.append(str(v).replace("|", "\\|"))
        lines.append("| " + " | ".join(items) + " |")
    return "\n".join(lines)


def optional_image(filename: str, caption: str) -> str:
    if (RESULTS / filename).is_file():
        return f"![{caption}](results/{filename})\n"
    return f"*Figure unavailable locally: `results/{filename}`.*\n"


def derive_switch(sensitivity: pd.DataFrame) -> float | None:
    needed = ["Alpha", "Score difference LFP-NMC811"]
    if not set(needed).issubset(sensitivity.columns):
        raise ValueError(f"Sensitivity data must include: {needed}")
    data = sensitivity[needed].copy().apply(pd.to_numeric, errors="raise").sort_values("Alpha")
    a = data["Alpha"].to_numpy(dtype=float)
    d = data["Score difference LFP-NMC811"].to_numpy(dtype=float)
    if len(a) < 2 or not np.isfinite(a).all() or not np.isfinite(d).all():
        raise ValueError("Invalid sensitivity data")
    for i in range(len(a) - 1):
        if abs(d[i]) < 1e-12:
            return float(a[i])
        if d[i] * d[i + 1] < 0:
            return float(a[i] - d[i] * (a[i + 1] - a[i]) / (d[i + 1] - d[i]))
    if abs(d[-1]) < 1e-12:
        return float(a[-1])
    return None


def inspect_dependencies() -> tuple[str, list[str]]:
    """Static best-effort import inventory; avoid deleting unknown dependencies."""
    standard = set(getattr(sys, "stdlib_module_names", set()))
    aliases = {
        "numpy": "numpy", "pandas": "pandas", "matplotlib": "matplotlib",
        "scipy": "scipy", "seaborn": "seaborn", "openpyxl": "openpyxl",
        "sklearn": "scikit-learn", "PIL": "Pillow", "cv2": "opencv-python",
        "plotly": "plotly", "tqdm": "tqdm", "statsmodels": "statsmodels",
    }
    imports: set[str] = set()
    parse_errors: list[str] = []
    for path in sorted(SCRIPTS.glob("*.py")):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        except (SyntaxError, UnicodeError) as exc:
            parse_errors.append(f"Could not inspect {path.name}: {exc}")
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                imports.add(node.module.split(".")[0])
    unknown = sorted(i for i in imports if i not in aliases and i not in standard and i not in {"__future__"})
    required = {aliases[i] for i in imports if i in aliases} | {"numpy", "pandas", "matplotlib"}
    lines = []
    for pkg in sorted(required):
        try:
            lines.append(f"{pkg}=={metadata.version(pkg)}")
        except metadata.PackageNotFoundError:
            lines.append(pkg)
            parse_errors.append(f"Version unavailable for {pkg}; emitted unpinned dependency")
    if unknown:
        parse_errors.append("Unrecognized third-party/local imports (review manually): " + ", ".join(unknown))
    return "\n".join(lines) + "\n", parse_errors


def backup_and_write(path: Path, content: str) -> None:
    if path.exists():
        BACKUPS.mkdir(parents=True, exist_ok=True)
        tag = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
        shutil.copy2(path, BACKUPS / f"{path.stem}_{tag}{path.suffix}")
    path.write_text(content, encoding="utf-8")


def main() -> None:
    print("=" * 72)
    print("STEP 25 v2 - AUTOMATED PROJECT DOCUMENTATION")
    print("=" * 72)
    data = {k: load_csv(k) for k in INPUTS}
    integrated = indexed(data["integrated"], "integrated", "Chemistry")
    criteria = indexed(data["criteria"], "criteria", "Chemistry")
    normalized = indexed(data["normalized"], "normalized", "Chemistry")
    scores = indexed(data["scores"], "scores", "Chemistry")
    weights = indexed(data["weights"], "weights", "Criterion")
    cells = data["cells"].copy()
    if not set(CHEMS).issubset(integrated.index):
        raise ValueError("Integrated summary must contain LFP and NMC811")
    for frame, name in [(criteria, "criteria"), (normalized, "normalized"), (scores, "scores")]:
        if not set(CHEMS).issubset(frame.index):
            raise ValueError(f"{name} must contain both chemistries")
    if not {"Chemistry", "Cell", "TR onset candidate [s]", "Propagation delay [s]"}.issubset(cells.columns):
        raise ValueError("Cell summary is missing required columns")
    for chem in CHEMS:
        if len(cells[cells["Chemistry"] == chem]) != 5:
            raise ValueError(f"Expected exactly five cell records for {chem}")

    key_cols = ["Total TR propagation [s]", "Mean inter-cell TR delay [s]",
                "Mean cell-associated Tmax [C]", "Mean venting Tmax [C]"]
    values = pd.DataFrame(index=key_cols, columns=CHEMS, dtype=float)
    for col in key_cols:
        for chem in CHEMS:
            values.loc[col, chem] = number(integrated, chem, col)
    values.index = ["Total TR propagation [s]", "Mean inter-cell delay [s]",
                    "Mean cell-associated Tmax [°C]", "Mean venting-sensor Tmax [°C]"]
    prop_lfp = number(integrated, "LFP", "Total TR propagation [s]")
    prop_nmc = number(integrated, "NMC811", "Total TR propagation [s]")
    if prop_lfp <= 0 or prop_nmc <= 0:
        raise ValueError("Propagation durations must be positive")
    ratio = prop_lfp / prop_nmc

    for chem in CHEMS:
        row = cells[cells["Chemistry"] == chem].sort_values("Cell")
        onsets = pd.to_numeric(row["TR onset candidate [s]"], errors="raise").to_numpy(dtype=float)
        delays = pd.to_numeric(row["Propagation delay [s]"], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(onsets).all() or not np.all(np.diff(onsets) > 0):
            raise ValueError(f"Onset times must be finite and increasing for {chem}")
        if not np.isclose(onsets[-1] - onsets[0], number(integrated, chem, "Total TR propagation [s]"), atol=0.15):
            raise ValueError(f"Cell onset times inconsistent with total propagation for {chem}")
        if not np.allclose(delays[1:], np.diff(onsets), atol=0.15):
            raise ValueError(f"Cell delays inconsistent with onset differences for {chem}")

    energy_cols = ["Gravimetric energy density [Wh/kg]", "Volumetric energy density [Wh/L]"]
    energy = criteria.loc[CHEMS, energy_cols].T
    scenario_cols = ["Safety-oriented", "Energy-density-oriented"]
    for col in scenario_cols:
        if col not in scores.columns or col not in weights.columns:
            raise ValueError(f"Missing MCDA scenario: {col}")
    if set(weights.index) != set(criteria.columns):
        raise ValueError("MCDA weight criteria do not match raw criteria")
    for scenario in scenario_cols:
        w = pd.to_numeric(weights[scenario], errors="raise")
        if not np.isclose(w.sum(), 1.0, atol=1e-8) or (w < 0).any():
            raise ValueError(f"Invalid weights for {scenario}")
        for chem in CHEMS:
            reconstructed = sum(number(normalized, chem, criterion) * float(w[criterion]) for criterion in weights.index)
            if not np.isclose(reconstructed, number(scores, chem, scenario), atol=2e-4):
                raise ValueError(f"MCDA score inconsistent: {chem}/{scenario}")
    switching = derive_switch(data["sensitivity"])
    switching_text = (f"**α = {switching:.4f}**" if switching is not None else
                      "No crossover found along the selected interpolation path")
    winners = {scenario: scores.loc[CHEMS, scenario].astype(float).idxmax() for scenario in scenario_cols}
    dependencies, dep_warnings = inspect_dependencies()
    if not SCRIPTS.is_dir():
        raise FileNotFoundError(f"Scripts directory missing: {SCRIPTS}")

    cell_view = cells.sort_values(["Chemistry", "Cell"]).set_index(["Chemistry", "Cell"])
    cell_cols = ["TR onset candidate [s]", "Propagation delay [s]", "Temperature sensor",
                 "Cell-associated Tmax [C]", "Venting sensor", "Venting Tmax [C]", "TR onset status"]
    cell_cols = [c for c in cell_cols if c in cell_view.columns]
    cell_table = cell_view[cell_cols].reset_index().set_index("Chemistry")
    raw_table = criteria.loc[CHEMS].T
    norm_table = normalized.loc[CHEMS].T
    score_table = scores.loc[CHEMS, scenario_cols]

    readme = f"""# Battery Thermal Safety Assessment | LFP vs NMC811

**Python-based secondary analysis of published automotive lithium-ion battery thermal-runaway experiments.**

## Project at a glance

This independent engineering project processes published experimental signals, estimates thermal-runaway (TR) onset events, compares propagation timing and thermal exposure, and evaluates the trade-off between thermal-safety indicators and cell-level energy density. **No new experiments were conducted.**

**Original study:** {SOURCE} DOI: `{ARTICLE_DOI}`. Dataset DOI: `{DATA_DOI}`.

### Main findings

{md_table(values)}

- LFP/NMC811 propagation-duration ratio: **{ratio:.2f}**. This is an elapsed-time ratio, **not** a thermal-front speed or guaranteed intervention time.
- Safety-prioritized illustrative MCDA: **{winners['Safety-oriented']}** ranks first.
- Energy-density-prioritized illustrative MCDA: **{winners['Energy-density-oriented']}** ranks first.
- TR onset times are **preliminary algorithmic candidates**, not fully manually validated events.

## Engineering workflow

1. Inspect and clean thermocouple and voltage signals; identify suspect sensor readings.
2. Identify preliminary TR onset candidates using temperature-rise and voltage-slope evidence.
3. Calculate delays between consecutive cell events and compare thermal severity.
4. Analyze venting thermocouples and spatially resolved infrared camera channels separately.
5. Produce IR heatmaps and 100/150/200 °C threshold-crossing diagnostics.
6. Integrate results and evaluate engineering implications.
7. Apply transparent weighted-sum MCDA to safety-oriented and energy-density-oriented scenarios.
8. Generate this README and the technical report automatically from the CSV outputs.

## Selected visualizations

{optional_image('integrated_experimental_comparison.png', 'Integrated experimental comparison')}
{optional_image('engineering_safety_comparison.png', 'Engineering safety indicators')}
{optional_image('technology_decision_matrix.png', 'Technology decision matrix')}
{optional_image('technology_sensitivity_analysis.png', 'MCDA sensitivity analysis')}

## Repository layout

```text
project-root/
├── data/                   # Original measurement files, if included
├── data_cleaned/           # Cleaned thermocouple signals
├── scripts/                # Numbered Python analysis scripts
├── results/                # Derived CSV tables and figures
├── README.md
├── FINAL_ENGINEERING_REPORT.md
└── requirements.txt
```

## Reproducing the analysis

Requires Python 3.10+ and access to the original dataset. On Windows PowerShell:

```powershell
python -m venv .venv
.\\.venv\\Scripts\\Activate.ps1
python -m pip install -r requirements.txt
```

Run the numbered scripts **in dependency order**, starting from raw-data inspection. To regenerate documentation after the preceding analysis has completed:

```powershell
python scripts/25_generate_project_documentation.py
```

**Note:** The requirements file is inferred from static imports in the local scripts and installed package versions. Inspect it before claiming a clean-environment reproduction; a full fresh-install test has not yet been performed.

## Scope and limitations

- Five prismatic cells per chemistry, **not electrically interconnected**; cell geometry and heater power differ between configurations.
- TR onset candidates have not been independently verified against every event.
- NMC811 T5 artifacts and LFP T1 damage timing limit thermocouple comparisons.
- IR surface measurements above approximately 250 °C carry increased uncertainty.
- Vent temperatures are not calorimetric heat-release or gas-volume measurements.
- Cell-level energy density does not equal pack-level energy density.
- MCDA weights are illustrative and the thermal criteria may be correlated.

For methodology, validation caveats and detailed results, see **[FINAL_ENGINEERING_REPORT.md](FINAL_ENGINEERING_REPORT.md)**.

## Attribution

All experimental measurements originate from the cited publication and dataset. Processing, visualization, derived metrics and MCDA are independent secondary work. Check the original dataset license and attribution conditions before redistributing raw measurements.
"""

    report = f"""# Final Engineering Report

## Battery Thermal Safety Assessment: LFP vs NMC811

**Type:** Independent secondary experimental-data analysis  
**Generated:** {STAMP}  
**Reference:** {SOURCE} DOI: `{ARTICLE_DOI}`. Dataset DOI: `{DATA_DOI}`.

## 1. Executive summary

The published five-cell arrangements exhibit distinct thermal-runaway behavior. The LFP configuration has an estimated total propagation duration of **{prop_lfp:.2f} s**, compared with **{prop_nmc:.2f} s** for NMC811, a duration ratio of **{ratio:.2f}**. NMC811 has higher measured peak thermocouple and venting-sensor temperatures in the analyzed records, while offering higher cell-level energy density. The MCDA ranking changes with the chosen application priorities. These conclusions are specific to the tested cells and conditions.

## 2. Research question and experimental basis

**Question:** How do observed thermal-runaway propagation, temperature severity and cell energy density influence an automotive battery-chemistry selection?

The original experiments used **five prismatic cells of each chemistry**, 100% state of charge, with external heating to initiate TR. Cells were **not electrically connected**. The heater power differed: **600 W for NMC811** and **1200 W for LFP**, reflecting different heater/contact geometries; the study reports a nominal heat flux of approximately **6.2 W/cm²** in both setups. Similar nominal heat flux does not make all boundary conditions equivalent. Cell sizes, capacities and geometry differ, so chemistry alone cannot explain every observed difference.

**Measurement modalities must remain separate:** cell thermocouples (T channels), cell voltages (U channels), vent thermocouples (V channels), and spatially resolved infrared camera profiles.

## 3. Data processing and analytical methodology

### 3.1 Inspection and cleaning

The project inspects time-series headers, sampling intervals and signal ranges, removes physically invalid negative thermocouple readings where applicable, and diagnoses implausible spikes and abrupt changes. These steps are documented in the numbered Python scripts; **not every anomaly is automatically corrected**. The NMC811 T5 thermocouple is excluded from the reported cell-associated Tmax average because of artifacts. LFP T1 is affected by reported sensor damage after approximately 185.9 s, and its maximum occurs close to that limit.

### 3.2 Preliminary TR onset identification

The reference study combines **dU/dt < 0** and a **simultaneous cell-surface temperature rate dT/dt > 10 K/s**. The Python workflow estimates candidates using voltage and temperature evidence, but smoothing, numerical differentiation, time alignment and sensor choice can affect exact onset times. **These event times are not independently validated against all original video/event annotations.** A late voltage collapse to zero is not automatically the onset of TR.

### 3.3 Propagation timing

For ordered candidate onset times t₁ ... t₅:

- Inter-cell delay: **Δtᵢ = tᵢ₊₁ − tᵢ** for i = 1...4.
- Total propagation duration: **Δt_total = t₅ − t₁**.
- Mean delay: **(Σ Δtᵢ) / 4**.

These quantities represent elapsed event times, **not** physical propagation velocities. They do not provide guaranteed detection or evacuation margins.

### 3.4 Thermal and venting severity

For each usable thermocouple, the workflow evaluates the maximum recorded temperature and related exposure indicators. The chemistry-level mean cell-associated Tmax is an **average of individual sensor maxima**, not a maximum of the average temperature trace and not a cell-core temperature. Usable sensor counts differ (LFP n=5; NMC811 n=4). Venting Tmax similarly summarizes local vent thermocouple maxima, not gas enthalpy, released heat or total gas volume.

### 3.5 Infrared surface analysis

Spatially resolved IR channels are analyzed as temperature-versus-time signals. Heatmaps and first crossings of **100, 150 and 200 °C** provide spatial heating diagnostics. The emissivity-control tape was rated to approximately **250 °C**; above this level, measurement uncertainty increases. Channel indices must not be converted into millimetres or m/s without independently verified spatial calibration. IR threshold-crossing spans are not equivalent to cell TR propagation durations.

### 3.6 Multi-criteria decision analysis (MCDA)

Five criteria are used: TR propagation duration (larger preferred), cell-associated Tmax (smaller), venting Tmax (smaller), gravimetric energy density (larger), and volumetric energy density (larger).

- Benefit criterion: **sᵢⱼ = xᵢⱼ / maxᵢ(xᵢⱼ)**.
- Cost criterion: **sᵢⱼ = minᵢ(xᵢⱼ) / xᵢⱼ**.
- Weighted score: **Sᵢ = Σⱼ wⱼ sᵢⱼ**, with **Σⱼ wⱼ = 1**.

Scores are relative to the compared alternatives. Weights are **illustrative engineering preferences**, not experimentally established risk weights. The three thermal criteria may partially overlap, so the MCDA must not be read as a validated safety index.

## 4. Integrated experimental results

{md_table(values)}

### 4.1 Cell-level onset candidates and sensor observations

{md_table(cell_table, 'Chemistry')}

**N/A** denotes a missing or excluded measurement, not zero. Cell 1 has no preceding inter-cell propagation delay.

### 4.2 Interpretation

LFP propagation is **{ratio:.2f} times longer** than NMC811 in the analyzed configuration. The NMC811 cell-associated mean Tmax exceeds the LFP mean by **{number(integrated, 'NMC811', 'Mean cell-associated Tmax [C]') - number(integrated, 'LFP', 'Mean cell-associated Tmax [C]'):.2f} °C**; the corresponding venting-sensor difference is **{number(integrated, 'NMC811', 'Mean venting Tmax [C]') - number(integrated, 'LFP', 'Mean venting Tmax [C]'):.2f} °C**. These differences describe measured local sensor indicators, not total heat release or a universal chemistry effect.

## 5. Comparison with published reference results

The reference article reports a **mean LFP inter-cell TR interval of approximately 460.8 s**. Our derived mean is **{number(integrated, 'LFP', 'Mean inter-cell TR delay [s]'):.2f} s**, a difference of **{number(integrated, 'LFP', 'Mean inter-cell TR delay [s]') - 460.8:+.2f} s**. This is a useful consistency check, **not independent proof** that each onset time is correct.

The article also reports mean maximum thermocouple temperatures of approximately **524.4 °C (LFP)** and **899.3 °C (NMC811)**. Our means are **{number(integrated, 'LFP', 'Mean cell-associated Tmax [C]'):.2f} °C** and **{number(integrated, 'NMC811', 'Mean cell-associated Tmax [C]'):.2f} °C**, respectively. **Do not interpret the numerical differences as measurement error without matching the article's sensor selection, exclusion rules and aggregation definition.** LFP T1 near damage and NMC811 T5 artifacts are particularly relevant.

The study's qualitative findings about reaction duration, ignition, mass loss and venting behavior must not be confused with direct calorimetry or measurements reproduced by this project. A complete independent event-by-event validation remains future work.

## 6. Technology decision matrix

### 6.1 Raw decision criteria

{md_table(raw_table, 'Criterion')}

The energy-density values refer to the **specific tested cells**, not complete automotive battery packs.

### 6.2 Normalized scores

{md_table(norm_table, 'Criterion', 4)}

### 6.3 Scenario weights

{md_table(weights, 'Criterion', 2)}

### 6.4 Scenario outcomes

{md_table(score_table, 'Chemistry', 4)}

The **safety-oriented** scenario favors **{winners['Safety-oriented']}**, while the **energy-density-oriented** scenario favors **{winners['Energy-density-oriented']}** under the selected weights.

### 6.5 Sensitivity analysis

The interpolation between the two specified weight vectors produces a decision crossover at **{switching_text}**. Alpha measures the position along this interpolation, **not** the percentage weight of energy density. The 101-point series is retained in `results/technology_sensitivity_analysis.csv` rather than repeated as a long report table. This one-dimensional sensitivity path does not exhaust all possible weighting combinations or quantify uncertainty in the experimental inputs.

{optional_image('technology_sensitivity_analysis.png', 'MCDA sensitivity analysis')}

## 7. Engineering implications

The tested LFP configuration shows slower sequential TR propagation and lower local peak-temperature indicators. These observations may be useful when considering thermal barriers, vent routing, monitoring and protection strategies. NMC811 offers higher reported cell-level gravimetric and volumetric energy density. Neither chemistry is universally superior: cell format, pack structure, cooling, electrical connections, cost, aging and charging requirements matter.

## 8. Limitations and uncertainty

1. **Experimental comparability:** Different cell geometries, capacities and heater powers; the tests are not chemistry-only controlled comparisons.
2. **Onset uncertainty:** Algorithmic TR candidates remain preliminary and are sensitive to derivative calculation, noise and time alignment.
3. **Thermocouple quality:** NMC811 T5 is excluded; LFP T1 is near its damage time; averages have unequal sensor counts.
4. **Venting measurements:** Sensor movement can affect maxima; vent temperatures do not quantify gas amount or released energy.
5. **Infrared uncertainty:** Tape is rated to ~250 °C; readings above that are less reliable.
6. **Spatial calibration:** No physically calibrated front speed in m/s is claimed.
7. **Statistical scope:** No repeated-experiment confidence intervals or full uncertainty propagation have been established.
8. **Pack-level extrapolation:** Five-cell, electrically disconnected arrangements do not establish pack-level hazard severity or performance.
9. **MCDA scope:** Illustrative weights, potentially correlated criteria, and omission of cost, life, fast charging and other system requirements.

## 9. Conclusions

1. Estimated LFP total TR propagation duration is **{prop_lfp:.2f} s**, compared with **{prop_nmc:.2f} s** for NMC811.
2. NMC811 shows higher observed local cell-associated and venting thermocouple peak-temperature means in this dataset.
3. NMC811 has higher reported **cell-level** gravimetric and volumetric energy density.
4. The illustrative decision ranking changes with engineering priorities; the MCDA is **decision support**, not a certified safety assessment.
5. The project remains a secondary analysis requiring event-level onset verification and further uncertainty work before stronger generalization.

## 10. Reproducibility and future work

The numbered scripts in `scripts/` form the analysis pipeline. Raw measurements, intermediate outputs, derived CSVs and figures should be retained for full reproduction. The dependency inventory is a static best-effort scan, not a substitute for testing in a fresh environment.

Priority next checks: manually validate all onset candidates against original event records; verify IR channel-to-distance mapping; test uncertainty sensitivity; rerun the complete pipeline from raw data in a clean virtual environment; assess broader MCDA weight combinations and correlated criteria.

## 11. Selected figures

{optional_image('integrated_experimental_comparison.png', 'Integrated experimental comparison')}
{optional_image('engineering_safety_comparison.png', 'Engineering safety assessment')}
{optional_image('technology_decision_matrix.png', 'Technology decision matrix')}

## 12. Attribution

**Article:** {SOURCE} DOI: `{ARTICLE_DOI}`.  
**Original dataset:** DOI `{DATA_DOI}`.  
The original experiments belong to the cited research team; this repository documents independent secondary analysis.
"""

    # Validation report is created only after all required data checks pass.
    validation_rows = [
        ("Required CSV inputs", "PASS", f"{len(INPUTS)} inputs loaded"),
        ("Chemistry coverage", "PASS", "LFP and NMC811 present"),
        ("Cell onset chronology", "PASS", "Five increasing onset candidates per chemistry"),
        ("Propagation arithmetic", "PASS", "Onset differences agree with integrated summary"),
        ("MCDA weight sums", "PASS", "Each scenario sums to one"),
        ("MCDA score reconstruction", "PASS", "Scores match normalized values and weights"),
        ("Independent TR onset validation", "PENDING", "Manual event-by-event review needed"),
        ("IR spatial calibration", "PENDING", "Physical front velocity not claimed"),
        ("Full clean-environment reproduction", "PENDING", "Not performed by documentation generator"),
    ]
    for warning in dep_warnings:
        validation_rows.append(("Dependency inventory", "REVIEW", warning))
    validation = pd.DataFrame(validation_rows, columns=["Check", "Status", "Details"])

    # Write outputs after all validation and rendering succeeds.
    outputs = {
        ROOT / "README.md": readme,
        ROOT / "FINAL_ENGINEERING_REPORT.md": report,
        ROOT / "requirements.txt": dependencies,
        RESULTS / "documentation_validation.csv": validation.to_csv(index=False),
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    for path, content in outputs.items():
        backup_and_write(path, content)

    print("\nGenerated files:")
    for path in outputs:
        print(f"  {path.relative_to(ROOT)} ({path.stat().st_size:,} bytes)")
    print(f"\nPropagation ratio: {ratio:.4f}")
    print(f"Safety-oriented winner: {winners['Safety-oriented']}")
    print(f"Energy-density-oriented winner: {winners['Energy-density-oriented']}")
    print(f"MCDA switching alpha: {switching:.4f}" if switching is not None else "No switching point found")
    if dep_warnings:
        print("\nDependency inventory warnings:")
        for warning in dep_warnings:
            print("  - " + warning)
    print("\nStep 25 v2 completed successfully.")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError, KeyError, OSError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)
