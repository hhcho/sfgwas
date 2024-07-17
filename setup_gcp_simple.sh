#!/bin/bash

cd ..

# Update package list
sudo apt-get update -y

# Install python3-pip, wget, and git
sudo apt-get install python3-pip wget git -y

# Download Go
wget https://golang.org/dl/go1.21.5.linux-amd64.tar.gz

# Remove any existing Go installation and extract the new one
sudo rm -rf /usr/local/go
sudo tar -C /usr/local -xzf go1.21.5.linux-amd64.tar.gz

# Update PATH and set alias for python
echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
echo 'alias python=python3' >> ~/.bashrc
source ~/.bashrc

# Install numpy using pip3
pip3 install numpy

# Install screen
sudo apt install -y screen

# Clone the repositories
git clone https://github.com/hhcho/mpc-core
git clone https://github.com/hcholab/lattigo.git

# Checkout the specific branch in lattigo repository
cd lattigo/
git checkout lattigo_pca
cd ..
