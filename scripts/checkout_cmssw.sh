#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

scramv1 project CMSSW CMSSW_9_0_0_pre1;
cd CMSSW_9_0_0_pre1/src
eval `scramv1 runtime -sh`

git cms-addpkg Validation/RecoTau
cd Validation/RecoTau
git remote add -f cms-tau-pog https://github.com/cms-tau-pog/cmssw.git
git fetch cms-tau-pog/CMSSW_9_0_X_tau-pog_validation-tools
git checkout cms-tau-pog/CMSSW_9_0_X_tau-pog_validation-tools
cd $CMSSW_BASE/src

git clone https://gitlab.cern.ch/knutzen/tau_validation_tools.git Validation/tau_validation_tools

scramv1 b -j 4
