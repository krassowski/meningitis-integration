#!/usr/bin/env bash

sudo apt-get install libproj-dev libgdal-dev  # Needed for some R dependencies
sudo apt-get install graphviz libgraphviz-dev graphviz-dev # (not used directly, required for static plots of nbpipeline)

pip3 install -r requirements.txt
cd helpers
wget https://raw.githubusercontent.com/krassowski/pyvenn/master/venn.py -O _venn.py
cd ..
Rscript install.R

cd thirdparty
wget https://data.broadinstitute.org/gsea-msigdb/gsea/software/desktop/3.0/gsea-3.0.jar
cd ..
