import FWCore.ParameterSet.Config as cms

hltL1P2GTCandTau = cms.EDFilter("L1P2GTCandidate",
    inputTag = cms.InputTag("l1tGTProducer", "CL2Taus"),
    MinE = cms.double( -1.0 ),
    MinPt = cms.double( 10.0 ),
    MinMass = cms.double( -1.0 ),
    MaxMass = cms.double( -1.0 ),
    MinEta = cms.double( -1.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
