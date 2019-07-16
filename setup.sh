#!/usr/bin/env bash

# see https://askubuntu.com/questions/1056630/r-3-5-0-not-working-on-ubuntu-18-04
sudo apt-get update

cat <<EOF | sudo tee /etc/apt/preferences.d/pin-r35
Package: r-*
Pin: release a=bionic-cran35
Pin: version 3.5*
Pin-Priority: 800

Package: r-cran-nlme
Pin: release a=bionic-cran35
Pin: version 3.1.139-1bionic0
Pin-Priority: 800
EOF

sudo apt-get update
sudo apt-get install r-base-core
sudo apt-get install r-cran-cluster
sudo apt-get install r-recommended

sudo apt-get install libproj-dev libgdal-dev  # Needed for some R dependencies
#sudo apt-get install r-recommended=3.5.3-1bionic_all
sudo apt-get install r-base r-base-dev
sudo apt-get install graphviz libgraphviz-dev graphviz-dev # (not used directly, required for static plots of nbpipeline)

pip3 install -r requirements.txt
cd helpers
wget https://raw.githubusercontent.com/krassowski/pyvenn/master/venn.py -o _venn.py
cd ..
sudo Rscript install.R
