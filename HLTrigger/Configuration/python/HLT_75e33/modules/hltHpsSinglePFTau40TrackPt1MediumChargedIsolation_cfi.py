import FWCore.ParameterSet.Config as cms

from ..modules.hltHpsDoublePFTau40TrackPt1MediumChargedIsolation_cfi import hltHpsDoublePFTau40TrackPt1MediumChargedIsolation

hltHpsSinglePFTau40TrackPt1MediumChargedIsolation = hltHpsDoublePFTau40TrackPt1MediumChargedIsolation.clone(
    MinN = 1
)