# LC–MS Missingness Analysis

This repository contains downstream analysis scripts for studying feature missingness in untargeted LC–MS metabolomics data. The project uses repeated technical injections to estimate how consistently each feature is detected after data processing.

Raw vendor files were converted to `.mzML` using [**ProteoWizard/MSConvert**](https://proteowizard.sourceforge.io/download.html), then processed using [**MassCube**](https://github.com/huaxuyu/masscube). 

## Overview

The main metric used in this project is the probability of detection, $$P_{\mathrm{detection}}$$. For each feature, $$P_{\mathrm{detection}}$$ is calculated as:

$$P_{\mathrm{detection}, i} = \frac{N_{\mathrm{detected}, i}}{N_{\mathrm{replicates}}}$$

where:

$$N_{\mathrm{detected}, i}$$ is the number of replicates in which the feature was detected and $$N_{\mathrm{replicates}}$$ the total number of technical replicates

In this project, a feature is considered detected in a processed feature table when its intensity is greater than zero.

## Processing workflow

The overall workflow is:

```text id="dk159l"
Raw vendor files
    ↓
mzML conversion using ProteoWizard/MSConvert
    ↓
Feature detection and alignment using MassCube
    ↓
MassCube feature tables
    ↓
Downstream missingness analysis in this repository
```

## Repository status

**XXXXXX CLEAN UP CODE AND UPLOAD**

The analysis scripts are currently being cleaned and reorganized. Script names and usage examples will be updated as the repository is finalized.

Planned structure:

```text id="51f7xi"
missingness/
│
├── README.md
├── data/
│   └── README.md
│
├── scripts/
│   ├── calculate_pdetection.py
│   ├── compare_concentration_tables.py
│   ├── extract_raw_signal.py
│   ├── analyze_local_environment.py
│   ├── regression_analysis.py
│   └── generate_figures.py
│
├── results/
│   └── README.md
│
└── figures/
    └── README.md
```

## Scripts

Placeholder descriptions are listed below. Detailed usage will be added after the scripts are standardized.

| Script                            | Purpose                                                             |
| --------------------------------- | ------------------------------------------------------------------- |
| `calculate_pdetection.py`         | Calculate (P_{\mathrm{detection}}) from MassCube feature tables.    |
| `compare_concentration_tables.py` | Compare aligned features across dilution/concentration tables.      |
| `extract_raw_signal.py`           | Check raw `.mzML` files for signal near expected feature locations. |
| `analyze_local_environment.py`    | Estimate local co-eluting MS1 signal density around each feature.   |
| `regression_analysis.py`          | Relate feature-level metrics to (P_{\mathrm{detection}}).           |
| `generate_figures.py`             | Generate manuscript and supplementary figures.                      |

## Input data

This repository expects processed feature tables from MassCube and, for selected analyses, converted `.mzML` files.

Example inputs:

```text id="c9gz5q"
feature_table.csv
feature_table_missing_annotated.csv
common_concentration_feature_table.csv
*.mzML
```

Exact input formats will be documented with the finalized scripts.

## Outputs

Typical outputs include:

```text id="xfatq7"
P_detection tables
aligned concentration comparison tables
raw-signal extraction summaries
local-environment summaries
regression results
figure panels
```

## Dependencies

Package versions will be finalized after script cleanup.

Main Python packages:

```text id="dg05t4"
pandas
numpy
matplotlib
seaborn
scipy
scikit-learn
statsmodels
pymzml
```

## Citation

Citation information will be added after manuscript submission or publication.

## License

License information will be added before public release.
