# Engineering Safety Assessment: LFP vs NMC811

## 1. Scope

This assessment is based on the experimental indicators
processed in Step 22.

The comparison concerns the tested five-cell arrangements,
not complete automotive battery packs.

TR onset times remain preliminary candidates.

## 2. Thermal Runaway Propagation

LFP total propagation duration: 1838.20 s.

NMC811 total propagation duration: 395.80 s.

LFP/NMC811 duration ratio: 4.64.

Mean inter-cell delay:

- LFP: 459.55 s.
- NMC811: 98.95 s.

### Engineering interpretation

The tested LFP arrangement exhibited slower TR propagation.

This may create more time for detection and protective actions.

However, the available results do not establish a validated
warning time, intervention time or evacuation margin.

Propagation duration must not be interpreted as a physical
thermal-front velocity.

## 3. Thermal Severity

Mean cell-associated thermocouple Tmax:

- LFP: 534.63 °C (n=5).
- NMC811: 891.54 °C (n=4).

Difference: 356.91 °C.

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

- LFP: 297.86 °C.
- NMC811: 900.73 °C.

Difference: 602.87 °C.

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
