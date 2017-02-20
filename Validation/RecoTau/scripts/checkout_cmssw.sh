#!/bin/bash
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuidePFTauIDDevelopers#How_to_set_up_CMSSW_for_my_devel

export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

scramv1 project CMSSW CMSSW_9_0_0_pre1;
cd CMSSW_9_0_0_pre1/src
eval `scramv1 runtime -sh`
git cms-init

git cms-addpkg Validation/RecoTau
git remote add tau-pog https://github.com/cms-tau-pog/cmssw.git
git remote add my-tau-pog https://github.com/thomas-mueller/cmssw.git
git fetch my-tau-pog CMSSW_9_0_X_tau-pog_validation-tools:CMSSW_9_0_X_tau-pog_validation-tools
git checkout CMSSW_9_0_X_tau-pog_validation-tools

git clone https://gitlab.cern.ch/knutzen/tau_validation_tools.git Validation/tau_validation_tools

scramv1 b -j 8
