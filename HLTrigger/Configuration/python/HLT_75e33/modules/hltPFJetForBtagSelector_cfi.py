import FWCore.ParameterSet.Config as cms

hltPFJetForBtagSelector = cms.EDFilter( "HLT1PFJet",
    saveTags = cms.bool( True ),
    inputTag = cms.InputTag( "hltAK4PFJetsCorrected" ),
    triggerType = cms.int32( 86 ),
    MinE = cms.double( -1.0 ),
    MinPt = cms.double( 30.0 ),
    MinMass = cms.double( -1.0 ),
    MaxMass = cms.double( -1.0 ),
    MinEta = cms.double( -1.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 1 )
)