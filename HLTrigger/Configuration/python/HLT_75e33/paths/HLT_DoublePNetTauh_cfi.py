import FWCore.ParameterSet.Config as cms

from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from ..sequences.HLTJetFlavourTagParticleTransformerSequencePF_cfi import *
from ..modules.hltL1SeedForDoublePuppiTau_cfi import *
from ..modules.hltL1sDoubleTauWIP_cfi import *
from ..sequences.HLTHgcalLocalRecoSequence_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTTrackingV61Sequence_cfi import *
from ..sequences.HLTLocalrecoSequence_cfi import *
from ..sequences.HLTRawToDigiSequence_cfi import *
from ..modules.hltDoublePFJets30ParTTauhTagMediumWPL2DoubleTau_cfi import *

HLT_DoublePNetTauh = cms.Path( 
    HLTBeginSequence + 
    hltL1SeedForDoublePuppiTau +
    HLTRawToDigiSequence +
    HLTHgcalLocalRecoSequence +
    HLTLocalrecoSequence +
    HLTTrackingV61Sequence +
    HLTMuonsSequence +
    HLTParticleFlowSequence +
    hltAK4PFPuppiJets +
    HLTAK4PFPuppiJetsReconstruction +
    HLTJetFlavourTagParticleTransformerSequencePF +
    hltDoublePFJets30ParTTauhTagMediumWPL2DoubleTau +
    HLTEndSequence
    )