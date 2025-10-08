# XMM_extended_source

**Automated pipeline for detecting and analyzing extended X-ray sources in XMM-Newton observations**

---

## Overview

`XMM_extended_source` is a Python and shell-based toolkit for **identifying, filtering, and analyzing extended X-ray sources** (such as galaxy clusters and diffuse emission) in **XMM-Newton** observations. It streamlines repetitive steps in the XMM Science Analysis System (SAS) and XMM-ESAS workflow, enabling reproducible large-scale studies of extended sources.

This repository provides:

* Modular scripts to automate event filtering, source region creation, and spectral extraction.
* Tools for multi-observation spectral analysis and background modeling.
* Integration with **XSPEC** for model fitting and result interpretation.

Originally developed for extended source studies in archival XMM data, the pipeline can be adapted for custom surveys or instrument calibration studies.

---

## Key Features

* ✅ End-to-end pipeline for extended source analysis
* 🧩 Modular design (each step can be run independently)
* 📁 Multi-observation stacking and joint fitting
* 🔬 Automated background estimation and region filtering
* ⚙️ Integration with XSPEC and XMM-ESAS
* 📊 Output of spectra, response files, and diagnostic plots

---

## Table of Contents

1. [Installation](#installation)
2. [Quickstart](#quickstart)
3. [Workflow](#workflow)
4. [Input and Output](#input-and-output)
5. [Code Structure](#code-structure)
6. [Configuration](#configuration)
7. [Example](#example)
8. [Performance and Limitations](#performance-and-limitations)
9. [References](#references)
10. [Contributing](#contributing)
11. [License](#license)
12. [Acknowledgements](#acknowledgements)

---

## Installation

### Requirements

* **Operating System**: Linux (recommended), macOS (with SAS configured)
* **Python**: ≥ 3.8
* **External tools**:

  * [XMM-Newton SAS](https://www.cosmos.esa.int/web/xmm-newton/sas)
  * [XMM-ESAS](https://www.cosmos.esa.int/web/xmm-newton/xmm-esas)
  * [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/)

### Python dependencies

Install required packages:

```bash
pip install -r requirements.txt
```

or manually:

```bash
pip install numpy astropy matplotlib scipy
```

### Environment setup

Before running the scripts, initialize your SAS environment:

```bash
source set_env.sh
```

This script should configure the SAS paths, calibration files (CCF), and environment variables necessary for XMM data processing.

---

## Quickstart

Assuming SAS and XMM data are properly configured:

```bash
bash run_extended_source_analysis.sh <obsid> <source_name>
```

This will:

1. Filter raw event files
2. Generate extended source regions
3. Extract spectra and responses
4. Perform spectral fitting in XSPEC
5. Save results to the output directory

You can modify `run_extended_source_analysis.sh` to adapt parameters to your target source.

---

## Workflow

A typical workflow involves the following scripts:

| Step | Script                          | Description                                         |
| ---- | ------------------------------- | --------------------------------------------------- |
| 1    | `run_filtering.sh`              | Filters event files, removes background flares      |
| 2    | `generate_regions.sh`           | Generates source and background extraction regions  |
| 3    | `spectral_analysis.py`          | Extracts and groups spectra, creates response files |
| 4    | `xspec_analysis.py`             | Runs XSPEC fits, stores best-fit parameters         |
| 5    | `multi_observation_analysis.py` | Combines results from multiple observations         |
| 6    | `common_functions.py`           | Utility library used across analysis scripts        |

Each script can be run independently or via `run_extended_source_analysis.sh` for full automation.

---

## Input and Output

### Input

* **XMM-Newton Observation Data Files (ODF)**
* **Event lists** (`*.EVLI.FIT`)
* **Calibration files** (via `set_env.sh`)
* **Region definitions** (optional custom DS9 regions)

Expected directory structure:

```
data/
 └── <OBS_ID>/
      ├── ODF/
      ├── event/
      ├── spectra/
      └── results/
```

### Output

* Cleaned event files (`*_clean_evt.fits`)
* Region files (`*.reg`)
* Spectra and response files (`.pha`, `.arf`, `.rmf`)
* XSPEC log and fit results (`.xcm`, `.txt`)
* Diagnostic plots (`.png`, `.pdf`)

All outputs are stored in `results/` subdirectories within each observation.

---

## Code Structure

```
XMM_extended_source/
│
├── common_functions.py         # Shared utilities (file handling, logging, plotting)
├── spectral_analysis.py        # Spectrum extraction and response creation
├── xspec_analysis.py           # XSPEC-based spectral fitting
├── multi_observation_analysis.py # Combines multiple obs results
│
├── run_filtering.sh            # Event filtering and cleaning
├── run_extended_source_analysis.sh  # Full pipeline runner
├── set_env.sh                  # Environment setup for SAS and CCF
│
└── README.md                   # This file
```

---

## Configuration

### Environment Variables

Ensure the following SAS-related environment variables are defined (usually handled in `set_env.sh`):

```bash
export SAS_DIR=/path/to/sas
export SAS_CCFPATH=/path/to/ccf
export SAS_ODF=/path/to/odf
export SAS_VERBOSITY=1
```

### Custom parameters

Scripts such as `spectral_analysis.py` and `xspec_analysis.py` may accept configuration options for:

* Energy bands (e.g., 0.3–10 keV)
* Region sizes
* Background model choice
* Fit statistic and energy range
* XSPEC model template (e.g., `apec + powerlaw`)

---

## Example

Example workflow on a single XMM observation:

```bash
# 1. Set environment
source set_env.sh

# 2. Run filtering
bash run_filtering.sh 0830190101

# 3. Run analysis
python spectral_analysis.py --obsid 0830190101 --src-reg src.reg --bkg-reg bkg.reg

# 4. Run XSPEC fit
python xspec_analysis.py --input results/0830190101/spectra

# 5. Combine multiple observations
python multi_observation_analysis.py --list obs_list.txt
```

---

## Performance and Limitations

* Designed for **moderate data volumes** (tens to hundreds of observations).
* Performance scales with the number of event files and size of extraction regions.
* Requires functional SAS installation and valid calibration files.
* Does **not** perform source detection — assumes regions are predefined or externally generated.
* Results depend on background modeling and region selection quality.

---

## References

* Snowden, S. L., et al. (2008), *The XMM-Newton Extended Source Analysis Software (ESAS)*, A&A, 478, 615
* Nevalainen, J., et al. (2005), *Background modeling for extended sources in XMM-Newton*, ApJ, 620, 128
* XMM-Newton SAS Documentation: [https://www.cosmos.esa.int/web/xmm-newton/sas](https://www.cosmos.esa.int/web/xmm-newton/sas)
* XSPEC User Guide: [https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/)

---

## Contributing

Contributions and improvements are welcome! To contribute:

1. Fork the repository
2. Create a new branch (`git checkout -b feature/your-feature`)
3. Commit changes with clear messages
4. Push and open a pull request

Please ensure any new features are documented and tested.

---

## License

This project is licensed under the **MIT License**.
See the [LICENSE](LICENSE) file for details.

---

## Acknowledgements

Developed by **Ruo Shang**.
This work was carried out at **Barnard College, Columbia University**, in collaboration with **Prof. Snow White**.
Data analysis builds upon the **XMM-Newton SAS** and **ESAS** toolkits provided by ESA.
