import FWCore.ParameterSet.Config as cms

from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from ..sequences.HLTJetFlavourTagParticleNetSequencePF_cfi import *
from ..modules.hltL1SeedForDoublePuppiTau_cfi import *
from ..modules.hltL1sDoubleTauWIP_cfi import *

HLT_DoublePNetTauh = cms.Path( 
    HLTBeginSequence + 
    hltL1SeedForDoublePuppiTau +
    hltAK4PFPuppiJets +
    HLTAK4PFPuppiJetsReconstruction +
    HLTJetFlavourTagParticleNetSequencePF +
    hltDoublePFJets30PNetTauhTagMediumWPL2DoubleTau +
    HLTEndSequence
    )