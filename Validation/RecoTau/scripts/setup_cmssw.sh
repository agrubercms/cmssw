#!/bin/bash

if [ $# -ne 1 ]
then
	echo
	echo "Usage ./setup_cmssw.sh <cmssw_version>"
	echo " -> will setup the correct cmssw environment"
	echo
	exit
else
	cd $TAU_VAL_BASE/CMSSW/
	cmsrel $1
	cd $1/src
	cmsenv
	git cms-addpkg Validation/RecoTau
	scram b -j6
	cd Validation/RecoTau/test
	mkdir TauID
fi
