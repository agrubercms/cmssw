import FWCore.ParameterSet.Config as cms

from ..modules.hltPreDoublePFTauHPS_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from ..sequences.HLTPFTauHPS_cfi import *
from ..modules.hltAK4PFJetsForTaus_cfi import *
from ..sequences.HLTHPSMediumChargedIsoPFTauSequence_cfi import *
from ..modules.hltHpsSelectedPFTausTrackPt1MediumChargedIsolation_cfi import *
from ..modules.hltHpsDoublePFTau40TrackPt1MediumChargedIsolation_cfi import *
from ..modules.hltTauPNet_debug_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from RecoTauTag.HLTProducers.applyL2TauTag import *
from ..sequences.HLTAK4PFJetsReconstruction_cfi import *
from ..sequences.HLTPFTauHPS_cfi import *
from ..modules.hltAK4PFJetsForTaus_cfi import *
from HLTrigger.Configuration.dev.P1_HLTJetFlavourTagParticleNetSequencePF_cfi import *
from HLTrigger.Configuration.dev.hltL1sDoubleTauBigOR_cfi import *

HLT_DoublePNetTauh = cms.Path( 
    HLTBeginSequence + 
    HLTParticleFlowSequence +
    hltAK4PFPuppiJets +
    HLTAK4PFPuppiJetsReconstruction +
    HLTAK4PFJetsReconstruction +
    hltAK4PFJetsForTaus +
    HLTJetFlavourTagParticleNetSequencePF +
    hltDoublePFJets30PNetTauhTagMediumWPL2DoubleTau + 
    HLTEndSequence
    )

"""
fragment.HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3_v1 = cms.Path( 
    fragment.HLTBeginSequence + 
    fragment.hltL1sDoubleTauBigOR + 
    fragment.hltPreDoublePNetTauhPFJet30MediumL2NNeta2p3 +
    fragment.HLTL2TauTagNNSequence +
    fragment.hltL2DoubleTauTagNNFilter +
    fragment.HLTAK4PFJetsSequence +
    fragment.HLTJetFlavourTagParticleNetSequencePF +
    fragment.hltDoublePFJets30PNetTauhTagMediumWPL2DoubleTau +
    fragment.HLTEndSequence )


HLT_DoublePNetTauh = cms.Path(
    HLTBeginSequence + 
    #hltPreDoublePFTauHPS +
    HLTParticleFlowSequence +
    HLTAK4PFJetsReconstruction +
    hltAK4PFJetsForTaus +
    HLTPFTauHPS +
    HLTHPSMediumChargedIsoPFTauSequence +
    hltHpsSelectedPFTausTrackPt1MediumChargedIsolation +
    hltHpsDoublePFTau40TrackPt1MediumChargedIsolation +
    HLTEndSequence
)"""
