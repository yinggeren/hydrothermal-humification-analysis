# Hydrothermal humification analysis

This repository contains the processed dataset and Python code used for the
Random Forest analysis linking hydrochar properties to the soil quality index
(SQI). The analysis corresponds to Fig. 5c in the manuscript.

## Repository contents

```text
data/
  tableS7_hydrochar_properties.csv
  tableS8_soil_quality_indicators.csv
  sqi_dataset_for_random_forest.csv
scripts/
  random_forest_sqi_no_ar.py
outputs/
  fig5c_feature_importance.csv
  fig5c_feature_importance.svg
requirements.txt
CITATION.cff
LICENSE
```

## Analysis scope

The machine-learning analysis used 150 paired hydrochar-soil records compiled
from published studies. The predictor variables were intrinsic hydrochar
properties:

- pH
- C (%)
- N (%)
- P (g/kg)
- K (g/kg)
- O/C
- H/C

Application rate (AR, %) was not included in the final model. The target
variable was SQI, calculated from soil physicochemical and biological
indicators as described in the Supplementary Information.

## How to run

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the Random Forest analysis:

```bash
python scripts/random_forest_sqi_no_ar.py
```

The script reads `data/tableS7_hydrochar_properties.csv` and
`data/tableS8_soil_quality_indicators.csv`, recalculates SQI, and writes:

```text
outputs/sqi_dataset_for_random_forest.csv
outputs/fig5c_feature_importance.csv
outputs/fig5c_feature_importance.svg
outputs/run_metadata.csv
```

## Expected feature importance

The default settings reproduce the Random Forest mean-decrease-in-impurity
importance used in the manuscript:

| Rank | Feature | Importance |
|---:|---|---:|
| 1 | C (%) | 0.4701 |
| 2 | pH | 0.1622 |
| 3 | O/C | 0.1037 |
| 4 | P (g/kg) | 0.0875 |
| 5 | N (%) | 0.0866 |
| 6 | H/C | 0.0575 |
| 7 | K (g/kg) | 0.0326 |

## Notes for reviewers

The repository intentionally does not include the unpublished manuscript PDF or
the full Supplementary Information Word file. Instead, the relevant tabular
data are provided as CSV files for transparency and reuse.
