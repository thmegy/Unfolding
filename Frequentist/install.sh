#!/bin/bash

if ! [[ $# -eq 1 ]] ; then
    echo "Usage: source ${0} <architecture>"
    echo "Where <architecture> is the architecture specification used by lsetup from cvmfs"
    echo "Example architecture: x86_64-centos7-gcc8-opt"
    return 1
fi

# store setup gcc + linux version for future setup and for setup on grid
echo ${1} > architecture.txt
export LCG_ARCH=${1}
# store absolute path to fbuenv (with traversed symlinks) for migating
# to migrate fbuenv we need to know the old absolute paths to overwrite in new location
# this is needed for running on grid, where we send the virtualenv in a tarball
echo `readlink -f $PWD` >> architecture.txt

shopt -s expand_aliases
if [[ -z ${ATLAS_LOCAL_ROOT_BASE+x} ]] ; then
  export ATLAS_LOCAL_ROOT_BASE='/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase'
fi
alias setupATLAS="source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh"

echo "> setupATLAS"
setupATLAS --quiet
echo "> Setup Python3+ROOT"
lsetup "lcgenv -p LCG_96python3 ${LCG_ARCH} ROOT" --quiet
export UNFOLDINGDIR=$PWD

echo "> Setting up virtualenv 'plenv'"
mkdir plenv
cd plenv
python3 -m venv .
cd ..
source plenv/bin/activate
echo "> Installing extra dependencies via pip"
pip install -r install_packages_pip.txt

echo "> Setup TRExFitter..."
git submodule init
git submodule update
cd TRExFitter
git submodule init
git submodule update
source setup.sh
mkdir -p build && cd build && cmake ../ && cmake --build ./ && cd ..
cd ..
