import FWCore.ParameterSet.Config as cms

hltParticleNetDiscriminatorsJetTags = cms.EDProducer( "BTagProbabilityToDiscriminator",
  discriminators = cms.VPSet( 
    cms.PSet(  name = cms.string( "BvsAll" ),
      numerator = cms.VInputTag( 'hltParticleNetONNXJetTags:probb' ),
      denominator = cms.VInputTag( 'hltParticleNetONNXJetTags:probb','hltParticleNetONNXJetTags:probc','hltParticleNetONNXJetTags:probuds','hltParticleNetONNXJetTags:probg' )
    ),
    cms.PSet(  name = cms.string( "CvsAll" ),
      numerator = cms.VInputTag( 'hltParticleNetONNXJetTags:probc' ),
      denominator = cms.VInputTag( 'hltParticleNetONNXJetTags:probb','hltParticleNetONNXJetTags:probc','hltParticleNetONNXJetTags:probuds','hltParticleNetONNXJetTags:probg' )
    ),
    cms.PSet(  name = cms.string( "TauhvsAll" ),
      numerator = cms.VInputTag( 'hltParticleNetONNXJetTags:probtauhp','hltParticleNetONNXJetTags:probtauhm' ),
      denominator = cms.VInputTag( 'hltParticleNetONNXJetTags:probb','hltParticleNetONNXJetTags:probc','hltParticleNetONNXJetTags:probuds','hltParticleNetONNXJetTags:probg','hltParticleNetONNXJetTags:probtauhp','hltParticleNetONNXJetTags:probtauhm' )
    ),
    cms.PSet(  name = cms.string( "CvsL" ),
      numerator = cms.VInputTag( 'hltParticleNetONNXJetTags:probc' ),
      denominator = cms.VInputTag( 'hltParticleNetONNXJetTags:probc','hltParticleNetONNXJetTags:probuds','hltParticleNetONNXJetTags:probg' )
    )
  )
)