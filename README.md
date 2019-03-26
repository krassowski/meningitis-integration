# Tuberculosis meningitis integration project

### Code organisation overview

There will be three types of files in the repository:
 - Jupyter Notebooks (`.ipynb` extension) which contain high-level analysis code, table summaries and visualisations.
  The idea is similar to Rmd files used in the R environment.
  The notebooks can be viewed as HTML files (in browser) by simply clicking on them in the repository.
 - Python modules (`.py` extension) with implementation details for the Python functions.
 - R scripts (`.R` extension) with implementation details for the R functions.

Jupyter Notebooks be default use the Python kernel, though some cells may contain R code.
This will be marked by the `%%R` at the beginning of the cell.
Similarly, an exclamation mark represents a bash command, e.g. `!ls`.

The generic helpers are stored in [helpers](helpers) directory, regardless of the implementation language.


## Setup

Clone with:

```bash
git clone https://github.com/krassowski/tuberculosis-meningitis-integration.git
```

Use pip to install Python requirements:

```bash
pip install -r requirements.txt
```

### Version compatibility

Developed and tested with Python 3.7.2 and R 3.5.1.
