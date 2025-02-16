import FWCore.ParameterSet.Config as cms

hltParticleTransformerONNXJetTags = cms.EDProducer( "hltParticleTransformerAK4ONNXJetTagsProducer",
    src = cms.InputTag( "hltParticleTransformerAK4TagInfos" ),
    model_path = cms.FileInPath( "ParticleTransformer.onnx" ),
    flav_names = cms.vstring(
        "probb",       # 'b'
        "probbb",      # 'bb'
        "problepb",    # 'leptonicB'
        "probtau",     # 'tau'
        "probc",       # 'c'
        "probuds"      # 'uds'
    ),
    input_names = cms.vstring(
        "global_features",
        "cpf",         # Changed from "cpf_features" to "cpf"
        "npf",
        "vtx"
    ),
    output_names = cms.vstring("output"),
    mightGet = cms.optional.untracked.vstring,
    
)

