import FWCore.ParameterSet.Config as cms

hltParticleTransformerDiscriminatorsJetTags = cms.EDProducer("BTagProbabilityToDiscriminator",
    discriminators = cms.VPSet(
      cms.PSet(name = cms.string("BvsAll"),
        numerator = cms.VInputTag(
            'hltParticleTransformerONNXJetTags:probb',
            'hltParticleTransformerONNXJetTags:problepb'
        ),
        denominator = cms.VInputTag(
            'hltParticleTransformerONNXJetTags:probb',
            'hltParticleTransformerONNXJetTags:problepb',
            'hltParticleTransformerONNXJetTags:probc',
            'hltParticleTransformerONNXJetTags:probuds'
        )
      ),
      cms.PSet(name = cms.string("CvsAll"),
        numerator = cms.VInputTag('hltParticleTransformerONNXJetTags:probc'),
        denominator = cms.VInputTag(
            'hltParticleTransformerONNXJetTags:probb',
            'hltParticleTransformerONNXJetTags:problepb',
            'hltParticleTransformerONNXJetTags:probc',
            'hltParticleTransformerONNXJetTags:probuds'
        )
      ),
      cms.PSet(name = cms.string("TauvsAll"),
        numerator = cms.VInputTag('hltParticleTransformerONNXJetTags:probtau'),
        denominator = cms.VInputTag(
            'hltParticleTransformerONNXJetTags:probb',
            'hltParticleTransformerONNXJetTags:problepb',
            'hltParticleTransformerONNXJetTags:probtau',
            'hltParticleTransformerONNXJetTags:probc',
            'hltParticleTransformerONNXJetTags:probuds'
        )
      ),
      cms.PSet(name = cms.string("CvsL"),
        numerator = cms.VInputTag('hltParticleTransformerONNXJetTags:probc'),
        denominator = cms.VInputTag('hltParticleTransformerONNXJetTags:probc',
                                     'hltParticleTransformerONNXJetTags:probuds')
      )
    )
)