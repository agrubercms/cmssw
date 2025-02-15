#include "RecoBTag/FeatureTools/interface/NeutralCandidateConverter.h"

namespace btagbtvdeep {

  void packedCandidateToFeatures(const pat::PackedCandidate* n_pf,
                                 const pat::Jet& jet,
                                 const bool isWeightedJet,
                                 const float drminpfcandsv,
                                 const float jetR,
                                 const float puppiw,
                                 NeutralCandidateFeatures& n_pf_features) {
    commonCandidateToFeatures(n_pf, jet, isWeightedJet, drminpfcandsv, jetR, puppiw, n_pf_features);

    n_pf_features.hadFrac = n_pf->hcalFraction();
    n_pf_features.puppiw = puppiw;
  }

  void recoCandidateToFeatures(const reco::PFCandidate* n_pf,
                               const reco::Jet& jet,
                               const bool isWeightedJet,
                               const float drminpfcandsv,
                               const float jetR,
                               const float puppiw,
                               NeutralCandidateFeatures& n_pf_features) {
    commonCandidateToFeatures(n_pf, jet, isWeightedJet, drminpfcandsv, jetR, puppiw, n_pf_features);
    n_pf_features.puppiw = puppiw;

    // need to get a value map and more stuff to do properly
    // otherwise will be different than for PackedCandidates
    // https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/packedPFCandidates_cfi.py
    if (abs(n_pf->pdgId()) == 1 || abs(n_pf->pdgId()) == 130) {
      n_pf_features.hadFrac = n_pf->hcalEnergy() / (n_pf->ecalEnergy() + n_pf->hcalEnergy());
    } else {
      n_pf_features.hadFrac = 0;
    }
  }

  // Renamed implementation for packed candidate neutral conversion
  void hltPackedCandidateToFeatures(const pat::PackedCandidate* n_pf,
                                       const pat::Jet& jet,
                                       const bool isWeightedJet,
                                       const float drminpfcandsv,
                                       const float jetR,
                                       const float puppiw,
                                       hltNpfCandidateFeatures& n_pf_features) {
    commonCandidateToFeatures(n_pf, jet, isWeightedJet, drminpfcandsv, jetR, puppiw, reinterpret_cast<NeutralCandidateFeatures&>(n_pf_features));
    n_pf_features.Npfcan_puppiw = puppiw;
    n_pf_features.Npfcan_HadFrac = n_pf->hcalFraction();
  }

  // Renamed implementation for reco candidate neutral conversion
  void hltRecoCandidateToFeatures(const reco::PFCandidate* n_pf,
                                     const reco::Jet& jet,
                                     const bool isWeightedJet,
                                     const float drminpfcandsv,
                                     const float jetR,
                                     const float puppiw,
                                     hltNpfCandidateFeatures& n_pf_features) {
    commonCandidateToFeatures(n_pf, jet, isWeightedJet, drminpfcandsv, jetR, puppiw, reinterpret_cast<NeutralCandidateFeatures&>(n_pf_features));
    n_pf_features.Npfcan_puppiw = puppiw;
    if (abs(n_pf->pdgId()) == 1 || abs(n_pf->pdgId()) == 130) {
      n_pf_features.Npfcan_HadFrac = n_pf->hcalEnergy() / (n_pf->ecalEnergy() + n_pf->hcalEnergy());
    } else {
      n_pf_features.Npfcan_HadFrac = 0;
    }
  }

  // Add wrapper implementations
  void convertNeutralCandidateToFeatures(const pat::PackedCandidate* n_pf,
                                        const pat::Jet& jet,
                                        bool isWeightedJet,
                                        float drminpfcandsv,
                                        float jetR,
                                        float puppiw,
                                        hltNpfCandidateFeatures& n_pf_features) {
    hltPackedCandidateToFeatures(n_pf, jet, isWeightedJet, drminpfcandsv, jetR, puppiw, n_pf_features);
  }

  void convertNeutralCandidateToFeatures(const reco::PFCandidate* n_pf,
                                        const reco::Jet& jet,
                                        bool isWeightedJet,
                                        float drminpfcandsv,
                                        float jetR,
                                        float puppiw,
                                        hltNpfCandidateFeatures& n_pf_features) {
    hltRecoCandidateToFeatures(n_pf, jet, isWeightedJet, drminpfcandsv, jetR, puppiw, n_pf_features);
  }

}  // namespace btagbtvdeep
