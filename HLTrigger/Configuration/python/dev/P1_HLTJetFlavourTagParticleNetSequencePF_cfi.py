import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFTauPrimaryVertexProducer_cfi import *
from HLTrigger.Configuration.dev.hltVerticesPF_cfi import *
from HLTrigger.Configuration.HLT_75e33.modules.hltDeepInclusiveMergedVerticesPF_cfi import *
from HLTrigger.Configuration.HLT_75e33.modules.hltDeepInclusiveSecondaryVerticesPF_cfi import *
from HLTrigger.Configuration.HLT_75e33.modules.hltDeepInclusiveVertexFinderPF_cfi import *
from HLTrigger.Configuration.HLT_75e33.modules.hltPrimaryVertexAssociation_cfi import *
from HLTrigger.Configuration.HLT_75e33.modules.hltDeepTrackVertexArbitratorPF_cfi import *


HLTJetFlavourTagParticleNetSequencePF = cms.Sequence( hltVerticesPF + 
    hltVerticesPFSelector + 
    hltVerticesPFFilter + 
    hltPFJetForBtagSelector + 
    hltPFJetForBtag + 
    hltDeepBLifetimeTagInfosPF + 
    hltDeepInclusiveVertexFinderPF + 
    hltDeepInclusiveSecondaryVerticesPF + 
    hltDeepTrackVertexArbitratorPF + 
    hltDeepInclusiveMergedVerticesPF + 
    hltPrimaryVertexAssociation + 
    hltParticleNetJetTagInfos + 
    hltParticleNetONNXJetTags + 
    hltParticleNetDiscriminatorsJetTags 
    )
