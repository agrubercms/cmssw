import FWCore.ParameterSet.Config as cms

hltParticleTransformerONNXJetTags = cms.EDProducer( "hltParticleTransformerAK4ONNXJetTagsProducer",
    src = cms.InputTag( "hltParticleTransformerAK4TagInfos" ),
    model_path = cms.FileInPath( "ParticleTransformer_4global.onnx" ),
    flav_names = cms.vstring(
        "probb",     
        "probc",  
        "probuds",   
        "probg",
        "probtaup",
        "probtaum",
    ),
    input_names = cms.vstring(
        "global",
        "cpf",         # Changed from "cpf_features" to "cpf"
        "npf",
        "vtx",
    ),
    output_names = cms.vstring("output"),
    mightGet = cms.optional.untracked.vstring,
    
)

