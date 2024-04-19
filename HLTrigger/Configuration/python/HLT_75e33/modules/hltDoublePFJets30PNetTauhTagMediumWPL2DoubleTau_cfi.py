import FWCore.ParameterSet.Config as cms

hltDoublePFJets30PNetTauhTagMediumWPL2DoubleTau = cms.EDFilter( "TauTagFilter",
    saveTags = cms.bool( True ),
    nExpected = cms.int32( 2 ),
    taus = cms.InputTag( "hltPFJetForBtag" ),
    tauTags = cms.InputTag( 'hltParticleNetDiscriminatorsJetTags','TauhvsAll' ),
    tauPtCorr = cms.InputTag( 'hltParticleNetONNXJetTags','ptcorr' ),
    seeds = cms.InputTag( "hltL1sDoubleTauWIP" ),
    seedTypes = cms.vint32( -100 ),
    selection = cms.string( "double t1 = 0.56, t2 = 0.47, t3 = 0.001, t4 = 0, x1 = 30, x2 = 100, x3 = 500, x4 = 1000; if (pt <= x1) return t1; if ((pt > x1) && (pt <= x2)) return (t2 - t1) / (x2 - x1) * (pt - x1) + t1; if ((pt > x2) && (pt <= x3)) return (t3 - t2) / (x3 - x2) * (pt - x2) + t2; if ((pt > x3) && (pt <= x4)) return (t4 - t3) / (x4 - x3) * (pt - x3) + t3; return t4;" ),
    minPt = cms.double( 30.0 ),
    maxEta = cms.double( 2.3 ),
    usePtCorr = cms.bool( True ),
    matchWithSeeds = cms.bool( True ),
    matchingdR = cms.double( 0.5 )
)