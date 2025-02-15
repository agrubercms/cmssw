import FWCore.ParameterSet.Config as cms

hltParticleTransformerAK4TagInfos = cms.EDProducer("hltParticleTransformerAK4TagInfoProducer",
    candidates = cms.InputTag("hltParticleFlowTmp"),
    fallback_puppi_weight = cms.bool(True),
    fallback_vertex_association = cms.bool(False),
    flip = cms.bool(False),
    is_weighted_jet = cms.bool(False),
    jet_radius = cms.double(0.4),
    jets = cms.InputTag("hltPFJetForBtag"),
    max_jet_eta = cms.double(2.5),
    max_sip3dsig_for_flip = cms.double(99999),
    mightGet = cms.optional.untracked.vstring,
    min_candidate_pt = cms.double(0.95),
    min_jet_pt = cms.double(15),
    puppi_value_map = cms.InputTag(""),
    secondary_vertices = cms.InputTag("hltDeepInclusiveMergedVerticesPF"),
    vertex_associator = cms.InputTag("hltPrimaryVertexAssociation","original"),
    vertices = cms.InputTag("hltPhase2PixelVertices")
)
