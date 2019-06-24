#!/usr/bin/env bash

sudo apt-get install libproj-dev libgdal-dev  # Needed for some R dependencies

pip3 install -r requirements.txt
cd helpers
wget https://raw.githubusercontent.com/krassowski/pyvenn/master/venn.py -o _venn.py
cd ..
Rscript install.R
