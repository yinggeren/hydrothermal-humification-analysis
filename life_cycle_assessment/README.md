# Life-cycle assessment calculations

This folder contains the reproducible life-cycle assessment (LCA) calculations
for the sludge-straw hydrothermal pathways reported in the manuscript and
Supplementary Information. All calculations are expressed per functional unit
of 1 t dry mixed feedstock.

## What this folder reproduces

- Scenario net GWP results for incineration, HC and PHC under OIS, EFS and CSS
  (Table S12 and Fig. 6b).
- Optimized-integration component contributions (Table S11).
- One-at-a-time sensitivity ranges for PHC (Table S13 and Fig. S10).
- PHC heat-recovery break-even thresholds (Table S14 and Fig. 6c).
- Monte Carlo uncertainty summary for PHC (Table S15 and Fig. S11).

## How to run

The calculation uses base R only.

```bash
Rscript run_lca.R
```

The script writes reproducible outputs to `outputs/` and stops with an error if
the reported values are not reproduced at manuscript rounding precision.

## Folder structure

```text
life_cycle_assessment/
  README.md
  run_lca.R
  R/
    lca_model.R
  data/
    scenario_components.csv
    scenario_assumptions.csv
    sensitivity_inputs.csv
    break_even_inputs.csv
    monte_carlo_summary_expected.csv
    expected_results.csv
  outputs/
    scenario_results.csv
    table_s11_contributions.csv
    sensitivity_results.csv
    break_even_results.csv
    monte_carlo_summary.csv
    validation_report.txt
  docs/
    lca_equations.md
```

## Important interpretation note

The Monte Carlo 95% interval is -264.0 to +24.0 kg CO2-eq t-1. Because this
interval crosses zero, the PHC pathway is likely, but not guaranteed, to remain
net negative under the stated uncertainty ranges.

## Data boundary

The folder does not distribute Ecoinvent datasets. The scenario-specific
electricity and natural-gas characterization factors used in the manuscript are
reported directly in `data/scenario_assumptions.csv`.

