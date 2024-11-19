import FWCore.ParameterSet.Config as cms

from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..modules.hltL1SeedForDoublePuppiTau_cfi import *
from ..modules.hltL1P2GTTau_cfi import *

L1P2GT_DoubleNNTau52 = cms.Path(HLTBeginSequence+hltL1P2GTTau+HLTEndSequence)
#L1P2GT_DoubleNNTau52 = cms.Path(HLTBeginSequence+hltL1SeedForDoublePuppiTau+hltL1P2GTCandTau+HLTEndSequence)
