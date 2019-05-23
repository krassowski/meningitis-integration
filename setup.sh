#!/usr/bin/env bash

sudo apt-get install libproj-dev libgdal-dev  # Needed for some R dependencies

pip3 install -r requirements.txt
Rscript install.R
