# LC–MS Missingness Analysis

This repository contains downstream analysis scripts for studying feature missingness in untargeted LC–MS metabolomics data. We used repeated technical injections to estimate how consistently each feature is detected after data processing and to calculate peak-level attributes that may explain stochastic detection.

Raw vendor files were converted to `.mzML` using [**ProteoWizard/MSConvert**](https://proteowizard.sourceforge.io/download.html), then processed using [**MassCube (v1.2.13)**](https://github.com/huaxuyu/masscube). We also tested this workflow on [MZmine (4.10.6)](https://github.com/mzmine/mzmine/releases/tag/v4.10.6) and [MS-DIAL (v5.5.260323)](https://github.com/systemsomicslab/MsdialWorkbench/releases). The scripts in this repository start from the raw aligned feature tables exported by each processing software and apply the downstream filtering, detection probability calculation, raw signal extraction, and peak attribute calculation.

The sample data included in this repository contains a small subset of files for testing the workflow. Each software folder currently includes 3 analytical sample files plus 1 method blank (MB) file. A small sample `mzML.zip` is also included. The full dataset will be available on [Zenodo](https://doi.org/10.5281/zenodo.21460829).

---

## 📌 Overview

The main metric used in this project is the probability of detection, $$P_{\mathrm{detection}}$$. For each feature $$i$$, $$P_{\mathrm{detection}}$$ is calculated as:

$$P_{\mathrm{detection}, i} = \frac{N_{\mathrm{detected}, i}}{N_{\mathrm{replicates}}}$$

where $$N_{\mathrm{detected}, i}$$ is the number of replicates in which the feature was detected and $$N_{\mathrm{replicates}}$$ the total number of technical replicates.

In this project, a feature is considered detected in a processed feature table when its intensity is greater than the configured detection threshold:

```python
INT_THRESHOLD = 1000
```

With 25 technical replicates, $$P_{\mathrm{detection}}$$ is reported on a 0 to 100 scale in 4% increments.

The extraction scripts also calculate peak attributes for each replicate, including intensity, area, peak width, scan count, signal to noise ratio (S/N), smoothness, sharpness, symmetry, local peak density, local rank, and RT shift. These attributes are used to examine how chromatographic behavior, local chemical environment, and signal variability relate to stochastic feature detection.

---

## 🧭 Processing Workflow

The overall workflow is:

```text
Raw vendor files
    ↓
mzML conversion using ProteoWizard/MSConvert
    ↓
Feature detection and alignment using MassCube, MZmine, or MS-DIAL
    ↓
Software-specific aligned feature tables
    ↓
Filtering, detection probability, and peak attribute extraction
    ↓
Feature detection and peak attribute summaries
```

---

## ⚙️ Included Filtering Steps

Users do **NOT** need to pre-calculate $$P_{\mathrm{detection}}$$ bins or manually filter the feature table before running these scripts. The scripts are designed to start from the raw software outputs. Filtering is performed inside the workflow before $$P_{\mathrm{detection}}$$ and peak attributes are calculated.

The included filters are:

```text
m/z >= 65
RT <= 23 min
max sample intensity >= 3 × max blank intensity
```

The m/z filter is included because the MS method was set to acquire ions above 65 m/z, but lower m/z features can still appear in software outputs.

You can edit these values near the top of each analysis script:

```python
MIN_MZ = 65
MAX_RT = 23
BLANK_MULTIPLIER = 3
INT_THRESHOLD = 1000
```

For example, to use a stricter blank filter:

```python
BLANK_MULTIPLIER = 5
```

---

## 📁 Repository Structure

```text
missingness/
│
├── README.md
│
├── scripts/
│   ├── masscube_analysis.py
│   ├── mzmine_analysis.py
│   ├── msdial_analysis.py
│   ├── peak_attributes.py
│   ├── run_masscube.sh
│   ├── run_msdial.sh
│   └── run_mzmine.sh
│
└── sample_data/
    ├── masscube_sample/
    │   ├── aligned_feature_table.csv
    │   └── single_files/
    │       ├── 100-1_P1-C-C-1_01_13558.txt
    │       ├── 100-2_P1-C-C-1_01_13559.txt
    │       ├── 100-3_P1-C-C-1_01_13560.txt
    │       └── MB_P1-A-4_01_13240.txt
    │
    ├── mzmine_sample/
    │   └── aligned_feature_table.csv
    │
    ├── msdial_sample/
    │   ├── feature_table_mb.csv
    │   ├── feature_table_no_mb.csv
    │   └── single_files/
    │       ├── 100-1_P1-C-C-1_01_13558.txt
    │       ├── 100-2_P1-C-C-1_01_13559.txt
    │       ├── 100-3_P1-C-C-1_01_13560.txt
    │       └── MB_P1-A-4_01_13240.txt
    │
    └── mzML/
        └── mzML.zip
```

The sample data are intentionally small so users can test whether the scripts run correctly. For full-scale analysis, replace the sample files with complete aligned feature tables, single-file outputs, and `.mzML` files.

The scripts use relative paths by default. If you run them from the repository root, they will look for input files inside `sample_data/`.

---

## 📦 Dependencies

Install the main dependencies with:

```bash
pip install pandas numpy pyteomics tqdm scipy matplotlib
```

---

## 🚀 Running the Workflows

Each software workflow has a small shell script in `scripts/`. The easiest way to run the analysis is to edit the paths at the top of the shell script, then run it from the [VS Code](https://code.visualstudio.com/thank-you?dv=win64user) terminal.

```text
scripts/
├── masscube_analysis.py
├── msdial_analysis.py
├── mzmine_analysis.py
├── peak_attributes.py
├── run_masscube.sh
├── run_msdial.sh
└── run_mzmine.sh
```

On Windows, use Git Bash or WSL in the VS Code terminal. Paths should use forward slashes (`/`).

For Git Bash on Windows:

```bash
PROJECT_DIR="E:/Code/Python/MassCube"
```

For WSL:

```bash
PROJECT_DIR="/mnt/e/Code/Python/MassCube"
```

Avoid Windows backslash paths inside `.sh` files:

```bash
PROJECT_DIR="E:\Nhi\Python\MassCube"
```

Each `.sh` file automatically finds the repository root, so it can be launched from VS Code without manually changing directories.

Run the workflows with:

```bash
bash scripts/run_masscube.sh
bash scripts/run_msdial.sh
bash scripts/run_mzmine.sh
```

Or, on macOS/Linux, make the scripts executable once:

```bash
chmod +x scripts/run_masscube.sh scripts/run_msdial.sh scripts/run_mzmine.sh
```

Then run:

```bash
./scripts/run_masscube.sh
./scripts/run_msdial.sh
./scripts/run_mzmine.sh
```

Inside each shell script, leave any variable empty to use the default path defined in the Python script. For example:

```bash
BLANK_MZML=""
```

means that `--blank-mzml` will not be passed, and the Python script will use its default blank mzML path.


### 🟢 MassCube

Run:

```bash
bash scripts/run_masscube.sh
```

Default input:

```text
sample_data/masscube_sample/aligned_feature_table.csv
sample_data/masscube_sample/single_files/
sample_data/mzML/
```

MassCube `single_files/*.txt` are converted automatically to `single_files_csv/*.csv`. Existing CSV files are reused unless reconversion is requested.

Default output:

```text
sample_data/masscube_sample/masscube_pdet/
```

Users can specify new paths by editing `scripts/run_masscube.sh`:

```bash
PROJECT_DIR="sample_data/masscube_sample" 
INPUT_CSV="" 
SINGLE_TXT_DIR="" 
SINGLE_CSV_DIR="" 
MZML_ROOT="sample_data/mzML" 
BLANK_MZML="" 
OUTPUT_ROOT=""
```


### 🔵 MS-DIAL

Run:
```bash
bash scripts/run_msdial.sh
```

MS-DIAL analysis in this project uses two aligned tables because it gap fills automatically and the `Fill %` from a blank-included table also includes the blank in the denominator. Since we calculate $$P_{\mathrm{detection}}$$ based on `Fill %`, this would be inaccurate.  Users would need to manually export individual traces into `single_files/`.

Inputs:

```text
feature_table_mb.csv        used for m/z, RT, and blank filtering
feature_table_no_mb.csv     used for Fill % and P_detection
single_files/               individual MS-DIAL peak lists
sample_data/mzML/           mzML files
```

Default output:

```text
sample_data/msdial_sample/msdial_pdet/
```

The feature table with MB analyzed is used for blank filtering, after which it's matched to the no-blank table to generate a list of filtered features with correct $$P_{\mathrm{detection}}$$. This is saved as:

```text
sample_data/msdial_sample/aligned_feature_table.csv
```

Users can edit these variables in `scripts/run_msdial.sh`:

```bash
PROJECT_DIR="sample_data/msdial_sample"
BLANK_INCLUDED_CSV=""
NOBLANK_CSV=""
INPUT_CSV=""
SINGLE_FILES_DIR=""
SINGLE_FILES_CSV_DIR=""
MZML_DIR="sample_data/mzML/100"
BLANK_MZML=""
OUTPUT_ROOT=""
```

### 🟡 MZmine

Run:

```bash
bash scripts/run_mzmine.sh
```

Default input:

```text
sample_data/mzmine_sample/aligned_feature_table.csv
sample_data/mzML/
```

MZmine stores sample-specific RT ranges directly in the aligned feature table, so no individual trace files are needed.

Expected columns include:

```text
rt
mz
id
datafile:<sample>.d:height
datafile:<sample>.d:area
datafile:<sample>.d:rt_range:min
datafile:<sample>.d:rt_range:max
```

Default output:

```text
sample_data/mzmine_sample/mzmine_pdet/
```

Users can edit these variables in `scripts/run_mzmine.sh`:

```bash
PROJECT_DIR="sample_data/mzmine_sample"
INPUT_CSV=""
MZML_DIR="sample_data/mzML/100"
BLANK_MZML=""
OUTPUT_ROOT=""
```

### 🤖 Advanced: Running Python Scripts

The shell scripts are wrappers around the Python scripts. Advanced users can still pass command line flags directly if they prefer. A list of MassCube flags is below (see software-specific scripts for details):

```text
--project-dir              project folder
--input-csv                aligned feature table
--single-txt-dir           MassCube single-file TXT folder
--single-csv-dir           converted single-file CSV folder
--mzml-root                mzML root folder
--blank-mzml               blank mzML file for S/N
--output-root              output folder
--skip-txt-conversion      skip TXT to CSV conversion
--overwrite-converted-csv  force TXT to CSV reconversion
```

For example:

```bash
python scripts/masscube_analysis.py --project-dir path/to/project --mzml-root path/to/mzML
```

---

## 📊 Output Files

Each analysis script creates a software-specific output folder inside `sample_data/*_sample` with similar structure. The example below expands on the MassCube output only:

```text
masscube_pdet/
├── masscube_pdet_full.csv
├── masscube_extraction_analysis.csv
├── masscube_pdet_summary.csv
├── feature_filter_report.csv
├── original_table.csv
├── attributes/
│   ├── area.csv
│   ├── density.csv
│   ├── int.csv
│   ├── ...
│   ├── symmetry.csv
│   ├── width.csv
├── 0/
├── 4/
├── ...
└── 100/

msdial_pdet/
└── same structure with msdial file prefixes

mzmine_pdet/
└── same structure with mzmine file prefixes
```

The main table, `*_pdet_full.csv`, contains only the analysis-ready attributes: `area`, `density`, `int`, `rank`, `rtshift`, `scancount`, `sharpness`, `smoothness`, `snr`, `symmetry`, and `width`. Each attribute is reported across replicates using columns such as `int_1`, `int_2`, and so on.

The `*_extraction_analysis.csv` file contains diagnostic information such as RT windows, scan IDs, scanlists, raw fallback values, intensity sources, excluded scans, and blank noise estimates.

The `attributes/` folder contains one file per attribute.

The numbered folders contain features grouped by \(P_{\mathrm{detection}}\), from `0/` to `100/`.

---

## 📈 Peak Attribute Summaries

After running one of the extraction workflows, summarize replicate-level attributes using:

```bash
python scripts/peak_attributes.py
```

Choose a software output:

```bash
python scripts/peak_attributes.py --model masscube
python scripts/peak_attributes.py --model msdial
python scripts/peak_attributes.py --model mzmine
```

Or provide a table directly:

```bash
python scripts/peak_attributes.py --input-csv path/to/masscube_pdet_full.csv
```

The output is written to `peak_analysis/` inside the software folder `*_sample/`:

```text
peak_analysis/
├── int_summary.csv
├── area_summary.csv
├── width_summary.csv
├── scancount_summary.csv
├── snr_summary.csv
├── smoothness_summary.csv
├── sharpness_summary.csv
├── symmetry_summary.csv
├── density_summary.csv
├── rank_summary.csv
└── rtshift_summary.csv
```

Each summary file contains feature metadata, the value for that attribute across replicates, and summary statistics: `max`, `min`, `mean`, `median`, `75_perc`, `25_perc`, `sd`, and `rsd`.

---

## 📝 Citation

Citation information will be added after manuscript submission or publication.

---

## 📄 License

This project is licensed under the Apache License 2.0. See the `LICENSE` file for details.
