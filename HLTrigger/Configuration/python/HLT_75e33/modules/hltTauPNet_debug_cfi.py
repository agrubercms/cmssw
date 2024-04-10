import FWCore.ParameterSet.Config as cms
from RecoTauTag.HLTProducers.applyL2TauTag import *
from HLTrigger.Configuration.dev.P1_BTagProbabilityToDiscriminator_cfi import *

hltDoublePFJets30PNetTauhTagMediumWPL2DoubleTau = cms.EDFilter( "TauTagFilter",
   matchWithSeeds = cms.bool( True ),
   matchingdR = cms.double( 0.5 ),
   maxEta = cms.double( 2.3 ),
   minPt = cms.double( 30 ),
   nExpected = cms.int32( 2 ),
   saveTags = cms.bool( True ),
   seedTypes = cms.vint32( [ -100 ] ),
   seeds = cms.InputTag( "hltHpsDoublePFTau40TrackPt1MediumChargedIsolation" ),
   selection = cms.string( "double t1 = 0.56, t2 = 0.47, t3 = 0.001, t4 = 0, x1 = 30, x2 = 100, x3 = 500, x4 = 1000; if (pt <= x1) return t1; if ((pt > x1) && (pt <= x2)) return (t2 - t1) / (x2 - x1) * (pt - x1) + t1; if ((pt > x2) && (pt <= x3)) return (t3 - t2) / (x3 - x2) * (pt - x2) + t2; if ((pt > x3) && (pt <= x4)) return (t4 - t3) / (x4 - x3) * (pt - x3) + t3; return t4;" ),
   tauPtCorr = cms.InputTag( "hltParticleNetONNXJetTags:ptcorr" ),
   tauTags = cms.InputTag( "hltParticleNetDiscriminatorsJetTags:TauhvsAll" ),
   taus = cms.InputTag( "hltAK4PFJetsForTaus" ),
   usePtCorr = cms.bool( True )
)