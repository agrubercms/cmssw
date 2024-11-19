import FWCore.ParameterSet.Config as cms

from ..modules.hltHpsDoublePFTau35MediumDitauWPDeepTau_cfi import hltHpsDoublePFTau35MediumDitauWPDeepTau

hltHpsSinglePFTau35MediumDitauWPDeepTauNoL1 = hltHpsDoublePFTau35MediumDitauWPDeepTau.clone(
    MinN = 1
)