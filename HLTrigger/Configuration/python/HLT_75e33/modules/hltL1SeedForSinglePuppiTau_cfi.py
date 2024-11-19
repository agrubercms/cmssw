import FWCore.ParameterSet.Config as cms

hltL1SeedForSinglePuppiTau = cms.EDFilter("PathStatusFilter",
    logicalExpression = cms.string('pSinglePuppiTau52')
)
