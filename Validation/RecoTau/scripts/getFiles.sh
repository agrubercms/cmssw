#!/usr/bin/bash

for campaign in 11_1_0_pre8
do
    echo python searchAndCopyInputFileFromWeb.py ${campaign} -n 10
done

echo python skimSamplesInFolder.py $CMSSW_BASE/src/Validation/RecoTau/data/relval_inputs/CMSSW_${campaign}_

python makeValPlots.py /afs/cern.ch/work/j/jaandrej/TauValidation2020/CMSSW_9_0_0_pre1/src/Validation/RecoTau/data/relval_skims/CMSSW_11_1_0_pre8_/DQM_V0001_R000000001__RelValZTT_14__CMSSW_11_1_0_pre8-111X_mcRun3_2021_realistic_v4_PU25nsPMX-v1__DQMIO.root /afs/cern.ch/work/j/jaandrej/TauValidation2020/CMSSW_9_0_0_pre1/src/Validation/RecoTau/data/relval_skims/CMSSW_11_1_0_pre8_/DQM_V0001_R000000001__RelValZTT_14__CMSSW_11_1_0_pre8-PU_111X_mcRun3_2021_realistic_v4-v1__DQMIO.root -t 1234_v2 -n 10