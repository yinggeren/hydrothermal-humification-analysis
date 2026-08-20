# Data description

`tableS7_hydrochar_properties.csv` contains the hydrochar-property records
compiled from the literature-derived dataset in Table S7.

`tableS8_soil_quality_indicators.csv` contains the corresponding soil-response
records in Table S8.

`sqi_dataset_for_random_forest.csv` is the paired modelling dataset used by
`scripts/random_forest_sqi_no_ar.py`. It contains the final predictor variables
and the calculated SQI target variable for 150 paired records.

Missing predictor values are mean-imputed during model fitting, following the
statistical procedure described in the manuscript.

