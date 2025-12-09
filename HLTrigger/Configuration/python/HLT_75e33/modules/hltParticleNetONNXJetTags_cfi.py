import FWCore.ParameterSet.Config as cms

hltParticleNetONNXJetTags = cms.EDProducer( "BoostedJetONNXJetTagsProducer",
    src = cms.InputTag( "hltParticleNetJetTagInfos" ),
    preprocess_json = cms.string( "preprocess.json" ),
    preprocessParams = cms.PSet(  ),
    model_path = cms.FileInPath( "PNet_notUnified_firstTry.onnx" ),
    flav_names = cms.vstring(
        "probb",     
        "probc",  
        "probuds",   
        "probg",
        "probtaup",
        "probtaum",
    ),
    jets = cms.InputTag( "" ),
    produceValueMap = cms.untracked.bool( False ),
    debugMode = cms.untracked.bool( True )
)
