# Meningitis integration with PLS and network analysis

[![Build Status](https://travis-ci.org/krassowski/meningitis-integration.svg?branch=master)](https://travis-ci.org/krassowski/meningitis-integration)

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


#### Automation for reproducibility

The order of execution of each of the notebooks and scripts is described by the Snakemake-like workflow in the [pipeline.py](pipeline.py) file.
Following a three-commands installation:

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

Execute nbpipline to reproduce all results of our study with:

```bash
# -i will generate an interactive graph featuring reproducibility reports
nbpipline pipeline.py -i
```

There is an excellent blog post on merits of a similar approach approach [here](http://ivory.idyll.org/blog/2018-repeatability-in-practice.html). Importantly we do not use snakemake as it is not ready yet to handle Jupyter notebooks the way we want to. 

Note for future: we could use snakemake [remote files](https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html) feature to move the data files out of this repository and and grant access to the data by providing authentication credential for the data access. 

### Analyses

For the list of conducted and planned analyses, please see the [Research_plan.md](Research_plan.md)

#### The analyses plot

For interactive plot use:
```bash
nbpipline pipeline.py -i
```

For static SVG plot use:
```bash
nbpipline pipeline.py -g
```

You can also append `-n` switch witch will skip execution of the pipeline and generate minimal plots (without reproducibility reports and code analysis):

For static SVG plot use:
```bash
nbpipline pipeline.py -g -n
```

### Version compatibility

Developed and tested with Python 3.7.2, R 3.5.1 and Ubuntu 19.04.

We declared all our Python dependencies in `requirements.txt` file and R dependencies in `install.R` file.
