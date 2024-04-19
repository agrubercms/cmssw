import FWCore.ParameterSet.Config as cms

hltDeepBLifetimeTagInfosPF = cms.EDProducer( "CandIPProducer",
    primaryVertex = cms.InputTag( "hltPhase2PixelVertices" ),
    computeProbabilities = cms.bool( True ),
    computeGhostTrack = cms.bool( True ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 3 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 5.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False ),
    jets = cms.InputTag( "hltPFJetForBtag" ),
    candidates = cms.InputTag("particleFlowTmp" ),
    maxDeltaR = cms.double( 0.4 )
)
