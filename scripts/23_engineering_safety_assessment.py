from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# STEP 23 - ENGINEERING SAFETY ASSESSMENT
# LFP vs NMC811
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

INPUT = RESULTS / "integrated_chemistry_summary.csv"

ASSESSMENT_OUTPUT = RESULTS / "engineering_safety_assessment.csv"
REPORT_OUTPUT = RESULTS / "engineering_safety_report.md"
FIGURE_OUTPUT = RESULTS / "engineering_safety_comparison.png"

CHEMISTRIES = ["LFP", "NMC811"]


# ============================================================
# 1. LOAD AND VALIDATE DATA
# ============================================================

print("=" * 75)
print("STEP 23 - ENGINEERING SAFETY ASSESSMENT")
print("=" * 75)

if not INPUT.exists():
    raise FileNotFoundError(
        f"Missing input file:\n{INPUT}\n"
        "Run step 22 before running step 23."
    )

df = pd.read_csv(INPUT)

required_columns = [
    "Chemistry",
    "Total TR propagation [s]",
    "Mean inter-cell TR delay [s]",
    "Mean cell-associated Tmax [C]",
    "Mean venting Tmax [C]",
    "Temperature sensors used",
]

missing = [c for c in required_columns if c not in df.columns]

if missing:
    raise ValueError(
        f"Missing required columns: {missing}"
    )

if df["Chemistry"].duplicated().any():
    raise ValueError("Duplicate chemistry rows detected.")

df = df.set_index("Chemistry")

if not set(CHEMISTRIES).issubset(df.index):
    raise ValueError("Both LFP and NMC811 must be present.")

numeric_columns = required_columns[1:]

for column in numeric_columns:
    df[column] = pd.to_numeric(df[column], errors="raise")

if df.loc[CHEMISTRIES, numeric_columns].isna().any().any():
    raise ValueError("Missing numerical values in input.")

if (df.loc[CHEMISTRIES, [
    "Total TR propagation [s]",
    "Mean inter-cell TR delay [s]",
]] <= 0).any().any():
    raise ValueError("Propagation times must be positive.")

print("\nInput validation passed.")


# ============================================================
# 2. EXTRACT EXPERIMENTAL RESULTS
# ============================================================

lfp = df.loc["LFP"]
nmc = df.loc["NMC811"]

lfp_prop = float(lfp["Total TR propagation [s]"])
nmc_prop = float(nmc["Total TR propagation [s]"])

lfp_delay = float(lfp["Mean inter-cell TR delay [s]"])
nmc_delay = float(nmc["Mean inter-cell TR delay [s]"])

lfp_temp = float(lfp["Mean cell-associated Tmax [C]"])
nmc_temp = float(nmc["Mean cell-associated Tmax [C]"])

lfp_vent = float(lfp["Mean venting Tmax [C]"])
nmc_vent = float(nmc["Mean venting Tmax [C]"])

lfp_n = int(lfp["Temperature sensors used"])
nmc_n = int(nmc["Temperature sensors used"])

prop_ratio = lfp_prop / nmc_prop
delay_ratio = lfp_delay / nmc_delay

temperature_difference = nmc_temp - lfp_temp
venting_difference = nmc_vent - lfp_vent

print("\nEXPERIMENTAL COMPARISON")
print("-" * 75)

print(f"LFP propagation: {lfp_prop:.2f} s")
print(f"NMC811 propagation: {nmc_prop:.2f} s")
print(f"Propagation time ratio: {prop_ratio:.2f}")

print(f"\nLFP mean inter-cell delay: {lfp_delay:.2f} s")
print(f"NMC811 mean inter-cell delay: {nmc_delay:.2f} s")
print(f"Inter-cell delay ratio: {delay_ratio:.2f}")

print(f"\nLFP mean sensor Tmax: {lfp_temp:.2f} °C")
print(f"NMC811 mean sensor Tmax: {nmc_temp:.2f} °C")

print(f"\nLFP mean venting Tmax: {lfp_vent:.2f} °C")
print(f"NMC811 mean venting Tmax: {nmc_vent:.2f} °C")


# ============================================================
# 3. ENGINEERING ASSESSMENT
# ============================================================

assessment = [
    {
        "Criterion": "TR propagation",
        "LFP observation": f"{lfp_prop:.2f} s",
        "NMC811 observation": f"{nmc_prop:.2f} s",
        "Experimental interpretation":
            f"LFP propagation duration is {prop_ratio:.2f} times longer.",
        "Engineering implication":
            "Slower propagation may provide more time for protective "
            "actions, but usable response time is not established.",
        "Limitation":
            "TR onset estimates remain preliminary. "
            "Cell geometries and heating powers differ."
    },
    {
        "Criterion": "Inter-cell propagation delay",
        "LFP observation": f"{lfp_delay:.2f} s",
        "NMC811 observation": f"{nmc_delay:.2f} s",
        "Experimental interpretation":
            f"LFP mean inter-cell delay is {delay_ratio:.2f} times longer.",
        "Engineering implication":
            "Thermal barriers and cooling strategies should be "
            "evaluated against chemistry-specific propagation times.",
        "Limitation":
            "Mean delay is based on four transitions per chemistry. "
            "It is not a physical front velocity."
    },
    {
        "Criterion": "Thermal severity",
        "LFP observation": f"{lfp_temp:.2f} °C (n={lfp_n})",
        "NMC811 observation": f"{nmc_temp:.2f} °C (n={nmc_n})",
        "Experimental interpretation":
            f"NMC811 mean sensor Tmax is {temperature_difference:.2f} °C higher.",
        "Engineering implication":
            "Higher measured temperatures indicate a greater "
            "local thermal load on nearby components under the "
            "tested conditions.",
        "Limitation":
            "Different numbers of thermocouples are used. "
            "NMC811 T5 is excluded; LFP T1 is close to its "
            "reported damage time. Tmax is not total heat release."
    },
    {
        "Criterion": "Venting severity",
        "LFP observation": f"{lfp_vent:.2f} °C",
        "NMC811 observation": f"{nmc_vent:.2f} °C",
        "Experimental interpretation":
            f"NMC811 mean venting-sensor Tmax is "
            f"{venting_difference:.2f} °C higher.",
        "Engineering implication":
            "Gas routing, vent placement, thermal shielding and "
            "ignition prevention are important design considerations "
            "for both chemistries.",
        "Limitation":
            "Sensor position can change during venting. "
            "Gas volume, composition and total released heat "
            "are not determined by this temperature metric."
    },
    {
        "Criterion": "Detection and response",
        "LFP observation": "Not directly quantified",
        "NMC811 observation": "Not directly quantified",
        "Experimental interpretation":
            "Propagation delays alone do not establish detection "
            "performance or intervention effectiveness.",
        "Engineering implication":
            "Evaluate detection thresholds, sensor placement, "
            "alarm latency, shutdown logic and venting response "
            "in dedicated system-level tests.",
        "Limitation":
            "No validated detection latency, suppression time "
            "or evacuation margin is available in step 22."
    },
    {
        "Criterion": "Battery pack integration",
        "LFP observation": "Module-level experimental evidence",
        "NMC811 observation": "Module-level experimental evidence",
        "Experimental interpretation":
            "Different thermal behaviors create different "
            "requirements for pack-level protection.",
        "Engineering implication":
            "Assess thermal barriers, cooling, gas management, "
            "structural protection and packaging together with "
            "energy-density requirements.",
        "Limitation":
            "The experimental cells were not electrically connected. "
            "Results cannot be transferred directly to a complete "
            "automotive battery pack."
    }
]

assessment_df = pd.DataFrame(assessment)
assessment_df.to_csv(ASSESSMENT_OUTPUT, index=False)

print("\nEngineering assessment table created.")


# ============================================================
# 4. COMPARATIVE FIGURE
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

labels = CHEMISTRIES

metrics = [
    (
        "Total TR propagation",
        [lfp_prop, nmc_prop],
        "Time [s]"
    ),
    (
        "Mean cell-associated Tmax",
        [lfp_temp, nmc_temp],
        "Temperature [°C]"
    ),
    (
        "Mean venting-sensor Tmax",
        [lfp_vent, nmc_vent],
        "Temperature [°C]"
    )
]

for ax, (title, values, ylabel) in zip(axes, metrics):

    bars = ax.bar(labels, values)

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)

    ax.set_ylim(0, max(values) * 1.18)

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{value:.1f}",
            ha="center",
            va="bottom"
        )

fig.suptitle(
    "Engineering Safety Assessment: LFP vs NMC811\n"
    "Experimental comparison; preliminary TR onset estimates",
    fontsize=13
)

fig.tight_layout(rect=[0, 0.02, 1, 0.88])

fig.savefig(FIGURE_OUTPUT, dpi=300, bbox_inches="tight")
plt.close(fig)

print("Comparative figure created.")


# ============================================================
# 5. ENGINEERING REPORT
# ============================================================

report = f"""# Engineering Safety Assessment: LFP vs NMC811

## 1. Scope

This assessment is based on the experimental indicators
processed in Step 22.

The comparison concerns the tested five-cell arrangements,
not complete automotive battery packs.

TR onset times remain preliminary candidates.

## 2. Thermal Runaway Propagation

LFP total propagation duration: {lfp_prop:.2f} s.

NMC811 total propagation duration: {nmc_prop:.2f} s.

LFP/NMC811 duration ratio: {prop_ratio:.2f}.

Mean inter-cell delay:

- LFP: {lfp_delay:.2f} s.
- NMC811: {nmc_delay:.2f} s.

### Engineering interpretation

The tested LFP arrangement exhibited slower TR propagation.

This may create more time for detection and protective actions.

However, the available results do not establish a validated
warning time, intervention time or evacuation margin.

Propagation duration must not be interpreted as a physical
thermal-front velocity.

## 3. Thermal Severity

Mean cell-associated thermocouple Tmax:

- LFP: {lfp_temp:.2f} °C (n={lfp_n}).
- NMC811: {nmc_temp:.2f} °C (n={nmc_n}).

Difference: {temperature_difference:.2f} °C.

### Engineering interpretation

The NMC811 arrangement exhibited higher measured
cell-associated peak temperatures.

This indicates more severe local thermal exposure under
the tested conditions.

The values are not measurements of total heat release.

NMC811 T5 was excluded because of known artifacts.

LFP T1 reached its reported maximum shortly before
the known sensor damage time.

The temperature means use different numbers of sensors.

## 4. Venting Severity

Mean venting-sensor Tmax:

- LFP: {lfp_vent:.2f} °C.
- NMC811: {nmc_vent:.2f} °C.

Difference: {venting_difference:.2f} °C.

### Engineering interpretation

The higher NMC811 venting-sensor temperatures indicate
greater measured local thermal exposure near the vents.

Both chemistries require appropriate vent-gas management.

Potential design measures include:

- Controlled gas routing.
- Thermal shielding of adjacent components.
- Prevention of gas accumulation.
- Assessment of ignition and pressure hazards.

The temperature measurements do not quantify total
gas volume, gas composition or total released heat.

## 5. Detection and Response

Propagation delays provide useful information about
the temporal development of a failure.

They do not independently determine the performance
of a battery management system.

Further validation is needed for:

- Detection thresholds.
- Sensor placement.
- Alarm latency.
- Electrical isolation.
- Protective response effectiveness.

No quantitative detection or intervention margin
is claimed in this assessment.

## 6. Battery Pack Integration

The tested cells were not electrically connected.

The chemistries also differed in cell geometry and
applied heater power.

Therefore, the experimental results cannot be directly
generalized to complete automotive battery packs.

Pack-level engineering should consider:

- Thermal barriers between cells.
- Cooling-system design.
- Vent-gas routing.
- Mechanical and structural protection.
- Electrical fault isolation.
- Energy-density and packaging trade-offs.

## 7. Measurement Limitations

1. TR onset times remain preliminary candidates.
2. LFP and NMC811 use different cell geometries.
3. The applied heater powers differ.
4. NMC811 T5 is excluded because of sensor artifacts.
5. LFP T1 has a known damage-time limitation.
6. Temperature averages use different sensor counts.
7. Venting-sensor positions may change during events.
8. IR measurements above 250 °C have increased uncertainty.
9. Surface threshold crossing spans are not TR durations.
10. Results represent the tested experimental configuration.

## 8. Overall Engineering Conclusion

Under the tested conditions, LFP exhibited slower
thermal runaway propagation and lower measured
cell-associated and venting-sensor temperatures
than NMC811.

These observations favor LFP with respect to the
specific thermal-safety indicators examined.

However, neither chemistry can be classified as
universally safe or unsafe from these results alone.

A complete technology selection requires additional
evaluation of energy density, system integration,
cost, performance and pack-level safety.
"""

REPORT_OUTPUT.write_text(report, encoding="utf-8")

print("Engineering report created.")


# ============================================================
# 6. FINAL SUMMARY
# ============================================================

print("\n" + "=" * 75)
print("FINAL ENGINEERING CONCLUSIONS")
print("=" * 75)

print(
    f"\n1. LFP TR propagation duration is "
    f"{prop_ratio:.2f} times longer."
)

print(
    f"2. NMC811 mean cell-associated Tmax is "
    f"{temperature_difference:.2f} °C higher."
)

print(
    f"3. NMC811 mean venting-sensor Tmax is "
    f"{venting_difference:.2f} °C higher."
)

print(
    "4. Detection and response margins are not "
    "quantified by the available data."
)

print(
    "5. Pack-level safety cannot be established "
    "from this experimental comparison alone."
)

print("\n" + "=" * 75)
print("SAVED FILES")
print("=" * 75)

print(ASSESSMENT_OUTPUT.name)
print(REPORT_OUTPUT.name)
print(FIGURE_OUTPUT.name)

print("\nStep 23 completed successfully.")