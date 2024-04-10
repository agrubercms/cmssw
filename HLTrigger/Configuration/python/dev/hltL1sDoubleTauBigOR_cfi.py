import FWCore.ParameterSet.Config as cms

from HLTrigger.Configuration.HLT_75e33.sequences.HLTDoLocalPixelSequence_cfi import *
from HLTrigger.Configuration.HLT_75e33.sequences.HLTDoLocalStripSequence_cfi import *
from HLTrigger.Configuration.HLT_75e33.tasks.localrecoTask_cfi import *
from HLTrigger.Configuration.dev.P1_hltL2TauTagNNProducer_cfi import *

hltGtStage2Digis = cms.EDProducer( "L1TRawToDigi",
    FedIds = cms.vint32( 1404 ),
    Setup = cms.string( "stage2::GTSetup" ),
    FWId = cms.uint32( 0 ),
    DmxFWId = cms.uint32( 0 ),
    FWOverride = cms.bool( False ),
    TMTCheck = cms.bool( True ),
    CTP7 = cms.untracked.bool( False ),
    MTF7 = cms.untracked.bool( False ),
    InputLabel = cms.InputTag( "rawDataCollector" ),
    lenSlinkHeader = cms.untracked.int32( 8 ),
    lenSlinkTrailer = cms.untracked.int32( 8 ),
    lenAMCHeader = cms.untracked.int32( 8 ),
    lenAMCTrailer = cms.untracked.int32( 0 ),
    lenAMC13Header = cms.untracked.int32( 8 ),
    lenAMC13Trailer = cms.untracked.int32( 8 ),
    debug = cms.untracked.bool( False ),
    MinFeds = cms.uint32( 0 )
)

hltL1sDoubleTauBigOR = cms.EDFilter( "HLTL1TSeed",
    saveTags = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleIsoTau32er2p1 OR L1_DoubleIsoTau34er2p1 OR L1_DoubleIsoTau35er2p1 OR L1_DoubleIsoTau36er2p1 OR L1_DoubleTau70er2p1" ),
    L1ObjectMapInputTag = cms.InputTag( "hltGtStage2ObjectMap" ),
    L1GlobalInputTag = cms.InputTag( "hltGtStage2Digis" ),
    L1MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    L1MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
    L1EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
    L1JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
    L1TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
    L1EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' ),
    L1EtSumZdcInputTag = cms.InputTag( 'hltGtStage2Digis','EtSumZDC' )
)

hltPreDoublePNetTauhPFJet30MediumL2NNeta2p3 = cms.EDFilter( "HLTPrescaler",
    offset = cms.uint32( 0 ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtStage2Digis" )
)

localrecoSequence = cms.Sequence( localrecoTask )

HLTL2TauTagNNSequence = cms.Sequence( HLTDoLocalPixelSequence + 
    localrecoSequence +
    hltL2TauTagNNProducer
    )


"""HLTRecopixelvertexingSequence + 
    HLTDoCaloSequence + 
    cms.ignore(hltL1sDoubleTauBigOR) + 
    cms.ignore(hltL1sSingleTau) + 
    cms.ignore(hltL1sBigOrMuXXerIsoTauYYer) + 
    cms.ignore(hltL1sMu22erIsoTau40er) + 
    cms.ignore(hltL1sBigORDoubleTauJet) + 
    cms.ignore(hltL1VBFDiJetIsoTau) + 
    cms.ignore(hltL1sVeryBigORMu18erTauXXer2p1) + 
    cms.ignore(hltL1sDoubleTauBigORWithLowMass) + """

