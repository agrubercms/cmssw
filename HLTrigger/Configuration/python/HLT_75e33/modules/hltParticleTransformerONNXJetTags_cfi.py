import FWCore.ParameterSet.Config as cms

hltParticleTransformerONNXJetTags = cms.EDProducer( "hltParticleTransformerAK4ONNXJetTagsProducer",
    src = cms.InputTag( "hltParticleTransformerAK4TagInfos" ),
    model_path = cms.FileInPath( "particle_transformer_bincl.onnx" ),
    flav_names = cms.vstring(
        "probb",     
        "problepb",  
        "probtau",   
        "probc",     
        "probuds",   
        "probg"      
    ),
    input_names = cms.vstring(
        "global_features",
        "cpf_features",         # Changed from "cpf_features" to "cpf"
        "npf_features",
        "vtx_features",
    ),
    output_names = cms.vstring("output"),
    mightGet = cms.optional.untracked.vstring,
    
)

