sudo sh -c 'echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list'
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update

# see https://askubuntu.com/questions/1056630/r-3-5-0-not-working-on-ubuntu-18-04
cat <<EOF | sudo tee /etc/apt/preferences.d/pin-r35
Package: r-*
Pin: release a=bionic-cran35
Pin: version 3.5*
Pin-Priority: 800

Package: r-cran-nlme
Pin: release a=bionic-cran35
Pin: version 3.1.139-1bionic0
Pin-Priority: 800

Package: r-cran-cluster
Pin: release a=bionic-cran35
Pin: version 2.0.8-1bionic0
Pin-Priority: 800
EOF


sudo apt-get install r-base r-base-dev

