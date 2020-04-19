# Meningitis omics integration with PLS, O2-PLS and LASSO

[![Build Status](https://travis-ci.org/krassowski/meningitis-integration.svg?branch=master)](https://travis-ci.org/krassowski/meningitis-integration)

### Scope

Data: Transcriptomics (mRNA, gene-level) and proteomics (SOMAScan) integration + clinical data including survival, meningitis subtype and HIV status.

For the list of analyses, please see the [Research_plan.md](Research_plan.md)

### Code organisation overview

Code is organised into:
 - Jupyter Notebooks (`.ipynb`) containing high-level analysis code, table summaries and visualisations.
  The notebooks can be viewed in the browser on the GitHub page of this repository.
 - Implementation details are stored in:
   - Python modules (`.py`)
   - R scripts (`.R`)

Most of the notebooks use the Python kernel, but some also contain R code. The R code in Python notebooks is marked by the `%%R` at the beginning of the cell (thanks to [rpy2](https://github.com/rpy2/rpy2)). An exclamation mark indicates a bash command, e.g. `!ls`.

The generic utilities and helper functions are stored in the [helpers](helpers) directory.


### Automation for reproducibility

The order of notebooks execution is described by the workflow rules in the [pipeline.py](pipeline.py) file.
Following the installation with:

```bash
# download the code
git clone git@github.com:krassowski/meningitis-integration.git
# enter the directory (if you use conda/virtualenv, activate it now)
cd meningitis-integration
# install Python and R dependencies
./setup.sh
# download the data
./download.sh
```

execute [nbpipeline](https://github.com/krassowski/nbpipeline) to reproduce all results of our study with:

```bash
# -i will generate an interactive graph featuring reproducibility reports
PYTHONPATH=$(pwd):$PYTHONPATH nbpipeline -i
```

#### The analyses plot

For interactive plot use `-i` (as shown above), for static use `-g`:

```bash
nbpipline pipeline.py -g
```

Append `-n` switch to skip execution of the pipeline and generate minimal plots (without reproducibility reports and code analysis):

```bash
nbpipline pipeline.py -i -n
```


### Version compatibility

Developed and tested with Python 3.7.6, R 3.6.3 and Ubuntu 19.10.

All dependencies are declared in:

  - `requirements.txt` for Python,
  - `install.R` for R


### About

The code in this repository was written as a part of MRes research project at Imperial College London.
