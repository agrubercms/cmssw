import FWCore.ParameterSet.Config as cms

from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..modules.hltL1SeedForDoublePuppiTau_cfi import *
from ..modules.hltL1P2GTCandTau_cfi import *
from ..modules.hltL1SeedForSinglePuppiTau_cfi import *

L1P2GT_SingleNNTau52 = cms.Path(HLTBeginSequence+hltL1SeedForSinglePuppiTau+hltL1P2GTCandTau+HLTEndSequence)
