import FWCore.ParameterSet.Config as cms

hltDeepInclusiveMergedVerticesPF = cms.EDProducer( "CandidateVertexMerger",
    secondaryVertices = cms.InputTag( "hltDeepTrackVertexArbitratorPF" ),
    maxFraction = cms.double( 0.2 ),
    minSignificance = cms.double( 10.0 )
)