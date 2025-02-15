#ifndef RecoBTag_FeatureTools_ChargedCandidateConverter_h
#define RecoBTag_FeatureTools_ChargedCandidateConverter_h

#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "DataFormats/BTauReco/interface/ChargedCandidateFeatures.h"
#include "DataFormats/BTauReco/interface/hltParticleTransformerAK4Features.h"  // for hltCpfCandidateFeatures
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

namespace btagbtvdeep {

  // Existing common candidate conversion (for non-HLT usage)
  template <typename CandidateType>
  void commonCandidateToFeatures(const CandidateType* c_pf,
                                 const reco::Jet& jet,
                                 const TrackInfoBuilder& track_info,
                                 const bool& isWeightedJet,
                                 const float& drminpfcandsv,
                                 const float& jetR,
                                 const float& puppiw,
                                 ChargedCandidateFeatures& c_pf_features,
                                 const bool flip = false,
                                 const float& distminpfcandsv = 0);

  // Add inline definition to resolve linker errors:
  template <typename CandidateType>
  inline void commonCandidateToFeatures(const CandidateType* c_pf,
                                        const reco::Jet& jet,
                                        const TrackInfoBuilder& track_info,
                                        const bool& isWeightedJet,
                                        const float& drminpfcandsv,
                                        const float& jetR,
                                        const float& puppiw,
                                        ChargedCandidateFeatures& c_pf_features,
                                        const bool flip,
                                        const float& distminpfcandsv) {
      // Compute subjet distance features:
      auto drSubjetFeatures = getDRSubjetFeatures(jet, c_pf);
      c_pf_features.drsubjet1 = drSubjetFeatures.first;
      c_pf_features.drsubjet2 = drSubjetFeatures.second;

      // Use weighted candidate kinematics:
      float constituentWeight = isWeightedJet ? puppiw : 1.0f;
      c_pf_features.ptrel = catch_infs_and_bound((c_pf->pt() * constituentWeight) / jet.pt(), 0, -1, 0, -1);
      c_pf_features.ptrel_noclip = (c_pf->pt() * constituentWeight) / jet.pt();
      c_pf_features.erel = (c_pf->energy() * constituentWeight) / jet.energy();

      // Only one deltaR member exists, so assign that
      c_pf_features.deltaR = catch_infs_and_bound(reco::deltaR(*c_pf, jet), 0, -0.6, 0, -0.6);
      
      // Removed assignment for deltaR_noclip:
      // c_pf_features.deltaR_noclip = reco::deltaR(*c_pf, jet);

      // Removed assignment for isGamma:
      // c_pf_features.isGamma = (std::abs(c_pf->pdgId()) == 22 ? 1 : 0);

      // Removed assignment for phirel:
      // c_pf_features.phirel = catch_infs_and_bound(std::fabs(reco::deltaPhi(c_pf->phi(), jet.phi())), 0, -2, 0, -0.5);
      
      c_pf_features.drminsv = catch_infs_and_bound(drminpfcandsv, 0, -jetR, 0, -jetR);
      c_pf_features.etarel = catch_infs_and_bound(std::fabs(c_pf->eta() - jet.eta()), 0, -2, 0, -0.5);
      c_pf_features.pt = c_pf->pt();
      c_pf_features.eta = c_pf->eta();
      c_pf_features.phi = c_pf->phi();
      c_pf_features.e = c_pf->energy();
      
      // Optionally use track_info and flip if needed.
      // ...existing conversion logic...
  }

  void packedCandidateToFeatures(const pat::PackedCandidate* c_pf,
                                 const pat::Jet& jet,
                                 const TrackInfoBuilder& track_info,
                                 const bool isWeightedJet,
                                 const float drminpfcandsv,
                                 const float jetR,
                                 const float puppiw,
                                 ChargedCandidateFeatures& c_pf_features,
                                 const bool flip = false,
                                 const float distminpfcandsv = 0);

  void recoCandidateToFeatures(const reco::PFCandidate* c_pf,
                               const reco::Jet& jet,
                               const TrackInfoBuilder& track_info,
                               const bool isWeightedJet,
                               const float drminpfcandsv,
                               const float jetR,
                               const float puppiw,
                               const int pv_ass_quality,
                               const reco::VertexRef& pv,
                               ChargedCandidateFeatures& c_pf_features,
                               const bool flip = false,
                               const float distminpfcandsv = 0);

  // --- New HLT overloads ---
  void hltPackedCandidateToFeatures(const pat::PackedCandidate* c_pf,
                                    const pat::Jet& jet,
                                    const TrackInfoBuilder& track_info,
                                    const bool isWeightedJet,
                                    const float drminpfcandsv,
                                    const float jetR,
                                    const float puppiw,
                                    hltCpfCandidateFeatures& c_pf_features,
                                    const bool flip,
                                    const float distminpfcandsv);

  void hltRecoCandidateToFeatures(const reco::PFCandidate* c_pf,
                                  const reco::Jet& jet,
                                  const TrackInfoBuilder& track_info,
                                  const bool isWeightedJet,
                                  const float drminpfcandsv,
                                  const float jetR,
                                  const float puppiw,
                                  const int pv_ass_quality,
                                  const reco::VertexRef& pv,
                                  hltCpfCandidateFeatures& c_pf_features,
                                  const bool flip,
                                  const float distminpfcandsv);

}  // namespace btagbtvdeep

#endif  //RecoBTag_FeatureTools_ChargedCandidateConverter_h
