#!/usr/bin/env bash

# you can add this script as init for AWS server

sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libseccomp-dev \
    libglib2.0-dev \
    uuid-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup \
    runc

export VERSION=1.19.2 OS=linux ARCH=amd64 && \
  wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
  sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
  rm go$VERSION.$OS-$ARCH.tar.gz

echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

git clone --recurse-submodules https://github.com/sylabs/singularity.git
cd singularity

git checkout --recurse-submodules v3.10.3

./mconfig && \
make -C ./builddir && \
sudo make -C ./builddir install
