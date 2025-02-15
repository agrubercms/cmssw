import FWCore.ParameterSet.Config as cms
#from HLTrigger.Configuration.dev.hltVerticesPF_cfi import *
from ..modules.hltDeepInclusiveMergedVerticesPF_cfi import *
from ..modules.hltDeepInclusiveSecondaryVerticesPF_cfi import *
from ..modules.hltDeepInclusiveVertexFinderPF_cfi import *
from ..modules.hltPrimaryVertexAssociation_cfi import *
from ..modules.hltDeepTrackVertexArbitratorPF_cfi import *
from ..modules.hltPFJetForBtagSelector_cfi import *
from ..modules.hltPFJetForBtag_cfi import *
from ..modules.hltDeepBLifetimeTagInfosPF_cfi import *
from ..modules.hltParticleTransformerAK4TagInfos_cfi import *
from ..modules.hltParticleTransformerONNXJetTags_cfi import *
from ..modules.hltParticleTransformerDiscriminatorsJetTags_cfi import *
from ..modules.hltDoublePFJets30PNetTauhTagMediumWPL2DoubleTau_cfi import *

HLTJetFlavourTagParticleTransformerSequencePF = cms.Sequence( 
    hltPFJetForBtagSelector + 
    hltPFJetForBtag + 
    hltDeepBLifetimeTagInfosPF + 
    hltDeepInclusiveVertexFinderPF + 
    hltDeepInclusiveSecondaryVerticesPF + 
    hltDeepTrackVertexArbitratorPF + 
    hltDeepInclusiveMergedVerticesPF + 
    hltPrimaryVertexAssociation + 
    hltParticleTransformerAK4TagInfos + 
    hltParticleTransformerONNXJetTags + 
    hltParticleTransformerDiscriminatorsJetTags 
    )
