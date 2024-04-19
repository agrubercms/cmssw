import FWCore.ParameterSet.Config as cms

hltParticleNetONNXJetTags = cms.EDProducer( "BoostedJetONNXJetTagsProducer",
    src = cms.InputTag( "hltParticleNetJetTagInfos" ),
    preprocess_json = cms.string( "RecoBTag/Combined/data/HLT/ParticleNetAK4/V01/preprocess.json" ),
    preprocessParams = cms.PSet(  ),
    model_path = cms.FileInPath( "RecoBTag/Combined/data/HLT/ParticleNetAK4/V01/particle-net.onnx" ),
    flav_names = cms.vstring( 'probtauhp',
      'probtauhm',
      'probb',
      'probc',
      'probuds',
      'probg',
      'ptcorr' ),
    jets = cms.InputTag( "" ),
    produceValueMap = cms.untracked.bool( False ),
    debugMode = cms.untracked.bool( False )
)