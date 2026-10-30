# Technology Decision Matrix: LFP vs NMC811

## Objective

Compare LFP and NMC811 battery cells using a
multi-criteria engineering decision framework.

## Data sources

Thermal-runaway indicators:
Derived from the experimental analysis in Steps 21-23.

Energy-density indicators:
Reported characteristics of the tested cells in
Schoberl et al., eTransportation 19 (2024), 100305.

## Method

Weighted-sum multi-criteria decision analysis.

Benefit criteria:
score = value / maximum value

Cost criteria:
score = minimum value / value

The normalized scores are relative to the two
technologies compared here.

## Scenario results

Safety-oriented:
LFP = 0.8968
NMC811 = 0.5449
Preferred technology = LFP

Energy-density-oriented:
LFP = 0.7592
NMC811 = 0.8146
Preferred technology = NMC811

## Sensitivity analysis

A ranking reversal occurs at alpha = 0.8639.

## Important limitations

1. TR onset estimates remain preliminary.

2. The experiments used different cell geometries
   and different heater powers.

3. Thermal-runaway propagation time is not equivalent
   to guaranteed detection or response time.

4. Temperature measurements have sensor-specific
   limitations.

5. The five criteria are not statistically independent.
   In particular, the three thermal-safety indicators
   may partially describe the same underlying behavior.

6. Energy densities refer to the tested cells,
   not complete automotive battery packs.

7. Criterion weights are illustrative engineering
   assumptions, not experimentally measured quantities.

8. Weighted scores are relative decision-support
   indicators, not absolute safety probabilities.

9. No uncertainty intervals or repeated-experiment
   variability are included in the weighted scores.

10. A real vehicle technology decision would require
    additional information about cost, cycle life,
    charging performance, operating temperatures,
    pack architecture and system-level safety.

## Engineering interpretation

The preferred battery chemistry depends on the
application requirements and chosen priorities.

This decision matrix demonstrates a transparent
trade-off between measured thermal-runaway behavior
and cell-level energy density.

It does not establish universal superiority
of either chemistry.
