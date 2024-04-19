import FWCore.ParameterSet.Config as cms

hltL1sDoubleTauWIP = cms.EDFilter( "HLTL1TSeed",
    saveTags = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleIsoTau32er2p1 OR L1_DoubleIsoTau34er2p1 OR L1_DoubleIsoTau35er2p1 OR L1_DoubleIsoTau36er2p1 OR L1_DoubleTau70er2p1" ),
    L1ObjectMapInputTag = cms.InputTag( "simGtStage2Digis" ),
    L1GlobalInputTag = cms.InputTag( "simGtStage2Digis" ),
    #L1MuonInputTag = cms.InputTag( 'simGmtStage2Digis'),
    #L1MuonShowerInputTag = cms.InputTag( 'simGtStage2Digis','MuonShower' ),
    #L1EGammaInputTag = cms.InputTag( 'simGtStage2Digis','EGamma' ),
    #L1JetInputTag = cms.InputTag( 'simGtStage2Digis','Jet' ),
    #L1TauInputTag = cms.InputTag( "simCaloStage2Digis" ),
    L1TauInputTag = cms.InputTag( 'l1tCaloJet', 'L1CaloTauCollectionBXV'),
    L1EtSumInputTag = cms.InputTag( 'simCaloStage2Digis' ),
    #L1EtSumZdcInputTag = cms.InputTag( 'simGtStage2Digis','EtSumZDC' )
)