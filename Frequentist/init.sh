#!/bin/bash

shopt -s expand_aliases
if [[ -z ${ATLAS_LOCAL_ROOT_BASE+x} ]] ; then
  export ATLAS_LOCAL_ROOT_BASE='/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase'
fi
alias setupATLAS="source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh"

export LCG_ARCH=`head -n1 architecture.txt`

echo "> setupATLAS"
setupATLAS --quiet
echo "> Setup Python3+ROOT"
lsetup "lcgenv -p LCG_96python3 ${LCG_ARCH} ROOT" --quiet
export UNFOLDINGDIR=$PWD

echo "> Loading virtualenv 'plenv'"
source plenv/bin/activate

echo "> Setup TRExFitter"
cd TRExFitter
source setup.sh
cd ..
