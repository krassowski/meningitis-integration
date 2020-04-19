# Meningitis omics integration with PLS, O2-PLS and LASSO

[![Build Status](https://travis-ci.org/krassowski/meningitis-integration.svg?branch=master)](https://travis-ci.org/krassowski/meningitis-integration)


### Scope

Data: Transcriptomics (mRNA, gene-level) and proteomics (SOMAScan) integration + clinical data including survival, meningitis subtype and HIV status.

For the list of analyses, please see the [Research_plan.md](Research_plan.md)

### Code organisation overview

There will be three types of files in the repository:
 - Jupyter Notebooks (`.ipynb`) which contain high-level analysis code, table summaries and visualisations.
  The idea is similar to R markdown (`.Rmd`) files.
  The notebooks can be viewed in the browser by clicking on them in the repository.
 - Implementation details are stored in:
   - Python modules (`.py`)
   - R scripts (`.R`)

Most of the notebooks use the Python kernel, but some also contain R code. The R code in Python notebooks is marked by the `%%R` at the beginning of the cell. An exclamation mark represents a bash command, e.g. `!ls`.

The generic helpers are stored in [helpers](helpers) directory, regardless of the implementation language.


#### Automation for reproducibility

The order of notebooks execution is described by the workflow rules in the [pipeline.py](pipeline.py) file.
Following the installation with:

```bash
# Download the code
git clone https://github.com/krassowski/meningitis-integration.git
# Enter the cloned directory; if you use conda/virtualenv, activate it now
cd meningitis-integration
# Install Python and R dependencies
./setup.sh
# Download the data that we depend on
./download.sh
```

execute nbpipline to reproduce all results of our study with:

```bash
# -i will generate an interactive graph featuring reproducibility reports
PYTHONPATH=$(pwd):$PYTHONPATH nbpipeline -i
```

#### The analyses plot

For interactive plot use:
```bash
nbpipline pipeline.py -i
```

For static SVG plot use:
```bash
nbpipline pipeline.py -g
```

Append `-n` switch to skip execution of the pipeline and generate minimal plots (without reproducibility reports and code analysis):

For static SVG plot use:
```bash
nbpipline pipeline.py -g -n
```

### Version compatibility

Developed and tested with Python 3.7.6, R 3.6.3 and Ubuntu 19.10.

All dependencies are declared in:

  - `requirements.txt` (Python)
  - `install.R` (R)


### About

The code in this repository was written as a part of MRes research project at Imperial College London.
