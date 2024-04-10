import FWCore.ParameterSet.Config as cms
from HLTrigger.Configuration.dev.hltVerticesPF_cfi import *

HLTJetFlavourTagParticleNetSequencePF = cms.Sequence( 
    hltVerticesPF + 
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
    hltParticleNetDiscriminatorsJetTags )