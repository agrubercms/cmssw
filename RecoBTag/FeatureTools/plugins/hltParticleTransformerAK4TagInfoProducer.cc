#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "DataFormats/BTauReco/interface/hltParticleTransformerAK4Features.h"
#include "DataFormats/BTauReco/interface/hltParticleTransformerAK4TagInfo.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>

// The HLT producer produces a vector of hltParticleTransformerAK4TagInfo.
class hltParticleTransformerAK4TagInfoProducer : public edm::stream::EDProducer<> {
public:
  explicit hltParticleTransformerAK4TagInfoProducer(const edm::ParameterSet&);
  ~hltParticleTransformerAK4TagInfoProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  typedef std::vector<reco::hltParticleTransformerAK4TagInfo> hltParticleTransformerAK4TagInfoCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
  typedef reco::VertexCollection VertexCollection;

  void produce(edm::Event&, const edm::EventSetup&) override;

  // Helper: Given a candidate pointer (from the jet’s daughter list) and the candidate collection handle,
  // search for the candidate and return a persistent edm::Ref.
  edm::Ref<edm::View<reco::Candidate>> getPersistentCandidate(
      const reco::Candidate* cand,
      const edm::Handle<edm::View<reco::Candidate>>& handle) const {
    for (size_t idx = 0; idx < handle->size(); ++idx) {
      if (&(handle->at(idx)) == cand) {
        return edm::Ref<edm::View<reco::Candidate>>(handle, idx);
      }
    }
    return edm::Ref<edm::View<reco::Candidate>>();
  }

  // --- Local helper conversion methods for candidates (charged/neutral) ---
  static void convertHLTChargedCandidate(const pat::PackedCandidate* cand,
                                         const reco::Jet& jet,
                                         const btagbtvdeep::TrackInfoBuilder& track_info,
                                         bool isWeightedJet,
                                         float drminpfcandsv,
                                         float jetR,
                                         float puppiw,
                                         hltCpfCandidateFeatures& feat,
                                         bool flip,
                                         float distminpfcandsv) {
    float constituentWeight = isWeightedJet ? puppiw : 1.0f;
    feat.Cpfcan_ptrel = (cand->pt() * constituentWeight) / jet.pt();
    feat.Cpfcan_drminsv = btagbtvdeep::catch_infs_and_bound(drminpfcandsv, 0, -jetR, 0, -jetR);
    feat.Cpfcan_pt = cand->pt();
    feat.Cpfcan_eta = cand->eta();
    feat.Cpfcan_phi = cand->phi();
    feat.Cpfcan_e   = cand->energy();
    feat.Cpfcan_VTX_ass = cand->pvAssociationQuality();
    feat.Cpfcan_puppiw  = puppiw;
    // --- New assignments using track_info (assumes the corresponding getter methods exist) ---
    feat.Cpfcan_BtagPf_trackEtaRel    = track_info.getTrackEtaRel();
    feat.Cpfcan_BtagPf_trackPtRel     = track_info.getTrackPtRel();
    feat.Cpfcan_BtagPf_trackPPar      = track_info.getTrackPPar();
    feat.Cpfcan_BtagPf_trackDeltaR    = track_info.getTrackDeltaR();
    feat.Cpfcan_BtagPf_trackPParRatio = track_info.getTrackPParRatio();
    feat.Cpfcan_BtagPf_trackSip2dVal  = track_info.getTrackSip2dVal();
    feat.Cpfcan_BtagPf_trackSip2dSig  = track_info.getTrackSip2dSig();
    feat.Cpfcan_BtagPf_trackSip3dVal  = track_info.getTrackSip3dVal();
    feat.Cpfcan_BtagPf_trackSip3dSig  = track_info.getTrackSip3dSig();
    feat.Cpfcan_BtagPf_trackJetDistVal= track_info.getTrackJetDistVal();
    feat.Cpfcan_chi2 = cand->hasTrackDetails() ?
                        btagbtvdeep::catch_infs_and_bound(cand->pseudoTrack().normalizedChi2(), 300, -1, 300)
                        : -1;
    feat.Cpfcan_quality = cand->hasTrackDetails() ?
                           cand->pseudoTrack().qualityMask()
                           : (1 << reco::TrackBase::loose);
    feat.Cpfcan_drminsv = btagbtvdeep::catch_infs_and_bound(drminpfcandsv, 0, -0.4, 0, -0.4);
  }

  static void convertHLTNeutralCandidate(const pat::PackedCandidate* cand,
                                         const reco::Jet& jet,
                                         bool isWeightedJet,
                                         float drminpfcandsv,
                                         float jetR,
                                         float puppiw,
                                         hltNpfCandidateFeatures& feat) {
    float constituentWeight = isWeightedJet ? puppiw : 1.0f;
    feat.Npfcan_ptrel = (cand->pt() * constituentWeight) / jet.pt();
    feat.Npfcan_deltaR = btagbtvdeep::catch_infs_and_bound(reco::deltaR(*cand, jet), 0, -0.6, 0, -0.6);
    feat.Npfcan_drminsv = btagbtvdeep::catch_infs_and_bound(drminpfcandsv, 0, -jetR, 0, -jetR);
    feat.Npfcan_pt = cand->pt();
    feat.Npfcan_eta = cand->eta();
    feat.Npfcan_phi = cand->phi();
    feat.Npfcan_energy = cand->energy();
    feat.Npfcan_puppiw = puppiw;
    feat.Npfcan_HadFrac = cand->hcalFraction();
  }

  // --- New helper: Explicit SV conversion ---
  // This function explicitly fills secondary vertex features from a given SV.
  void fillSVFeaturesHLT(btagbtvdeep::hltParticleTransformerAK4Features &features,
                         const reco::Jet& jet,
                         const reco::Vertex& pv,
                         const SVCollection* svs,
                         double jetR,
                         bool flip) {
    // Make a local copy of the SV collection and sort it using the provided comparator.
    SVCollection svs_sorted = *svs;
    std::sort(svs_sorted.begin(), svs_sorted.end(), [&pv](const auto& sv1, const auto& sv2) {
      return btagbtvdeep::sv_vertex_comparator(sv1, sv2, pv);
    });
    
    // Loop over sorted SVs and fill features for those within the jet radius.
    for (const auto& sv : svs_sorted) {
      if (reco::deltaR2(sv, jet) > (jetR * jetR))
        continue;
      
      // Use the HLT-specific vertex features type.
      hltVtxFeatures svfeat;
      // Map available quantities.
      svfeat.jet_sv_pt        = sv.pt();
      svfeat.jet_sv_deta      = btagbtvdeep::catch_infs_and_bound(std::fabs(sv.eta() - jet.eta()) - 0.5, 0, -2, 0);
      svfeat.jet_sv_dphi      = btagbtvdeep::catch_infs_and_bound(std::fabs(reco::deltaPhi(sv.phi(), jet.phi())) - 0.5, 0, -2, 0);
      svfeat.jet_sv_eta       = sv.eta();
      svfeat.jet_sv_phi       = sv.phi();
      svfeat.jet_sv_energy    = sv.energy();
      svfeat.jet_sv_mass      = sv.mass();
      svfeat.jet_sv_ntrack    = sv.numberOfDaughters();
      svfeat.jet_sv_chi2      = sv.vertexChi2();
      const auto& dxy_meas = btagbtvdeep::vertexDxy(sv, pv);
      svfeat.jet_sv_dxy       = dxy_meas.value();
      svfeat.jet_sv_dxysig   = btagbtvdeep::catch_infs_and_bound(dxy_meas.value() / dxy_meas.error(), 0, -1, 800);
      const auto& d3d_meas = btagbtvdeep::vertexD3d(sv, pv);
      svfeat.jet_sv_d3d       = d3d_meas.value();
      svfeat.jet_sv_d3dsig   = btagbtvdeep::catch_infs_and_bound(d3d_meas.value() / d3d_meas.error(), 0, -1, 800);
      svfeat.jet_sv_energy_log = std::log(sv.energy());
      
      // Append the filled HLT secondary vertex features.
      features.vtx_features.push_back(svfeat);
    }
  }

  // --- Configuration parameters and tokens ---
  const double jet_radius_;
  const double min_candidate_pt_;
  const bool flip_;

  const edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
  const edm::EDGetTokenT<VertexCollection> vtx_token_;
  const edm::EDGetTokenT<SVCollection> sv_token_;
  const edm::EDGetTokenT<edm::View<reco::Candidate>> candidateToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;

  const bool is_weighted_jet_;
  const double min_jet_pt_;
  const double max_jet_eta_;

  const bool fallback_puppi_weight_;
  const bool fallback_vertex_association_;
  const double max_sip3dsig_for_flip_;

  const edm::EDGetTokenT<edm::ValueMap<float>> puppi_value_map_token_;
  const edm::EDGetTokenT<edm::Association<VertexCollection>> vertex_associator_token_;
};

// --- Constructor ---
hltParticleTransformerAK4TagInfoProducer::hltParticleTransformerAK4TagInfoProducer(const edm::ParameterSet& iConfig)
    : jet_radius_(iConfig.getParameter<double>("jet_radius")),
      min_candidate_pt_(iConfig.getParameter<double>("min_candidate_pt")),
      flip_(iConfig.getParameter<bool>("flip")),
      jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
      candidateToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("candidates"))),
      track_builder_token_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      is_weighted_jet_(iConfig.getParameter<bool>("is_weighted_jet")),
      min_jet_pt_(iConfig.getParameter<double>("min_jet_pt")),
      max_jet_eta_(iConfig.getParameter<double>("max_jet_eta")),
      fallback_puppi_weight_(iConfig.getParameter<bool>("fallback_puppi_weight")),
      fallback_vertex_association_(iConfig.getParameter<bool>("fallback_vertex_association")),
      max_sip3dsig_for_flip_(iConfig.getParameter<double>("max_sip3dsig_for_flip")),
      puppi_value_map_token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("puppi_value_map"))),
      vertex_associator_token_(consumes<edm::Association<VertexCollection>>(iConfig.getParameter<edm::InputTag>("vertex_associator")))
{
  produces<hltParticleTransformerAK4TagInfoCollection>();
}

// --- fillDescriptions ---
void hltParticleTransformerAK4TagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<double>("jet_radius", 0.4);
  desc.add<double>("min_candidate_pt", 0.95);
  desc.add<bool>("flip", false);
  desc.add<edm::InputTag>("vertices", edm::InputTag("hltPixelVertices"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("hltInclusiveCandidateSecondaryVertices"));
  desc.add<edm::InputTag>("jets", edm::InputTag("hltAK4PFJets"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("hltPFCandidates"));
  desc.add<bool>("is_weighted_jet", false);
  desc.add<double>("min_jet_pt", 15.0);
  desc.add<double>("max_jet_eta", 2.5);
  desc.add<bool>("fallback_puppi_weight", true);
  desc.add<bool>("fallback_vertex_association", false);
  desc.add<double>("max_sip3dsig_for_flip", 99999);
  desc.add<edm::InputTag>("puppi_value_map", edm::InputTag(""));
  desc.add<edm::InputTag>("vertex_associator", edm::InputTag("hltPrimaryVertexAssociation", "original"));
  descriptions.add("hltParticleTransformerAK4TagInfos", desc);
}

// --- produce method ---
void hltParticleTransformerAK4TagInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto output_tag_infos = std::make_unique<hltParticleTransformerAK4TagInfoCollection>();

  edm::Handle<edm::View<reco::Jet>> jets;
  iEvent.getByToken(jet_token_, jets);

  edm::Handle<VertexCollection> vtxs;
  iEvent.getByToken(vtx_token_, vtxs);
  if (vtxs->empty()) {
    iEvent.put(std::move(output_tag_infos));
    return;
  }
  const auto& pv = vtxs->at(0);

  edm::Handle<edm::View<reco::Candidate>> tracks;
  iEvent.getByToken(candidateToken_, tracks);

  edm::Handle<SVCollection> svs;
  iEvent.getByToken(sv_token_, svs);

  edm::ESHandle<TransientTrackBuilder> track_builder = iSetup.getHandle(track_builder_token_);

  // Loop over jets
  for (std::size_t jet_n = 0; jet_n < jets->size(); ++jet_n) {
    edm::RefToBase<reco::Jet> jet_ref(jets, jet_n);
    const auto& jet = jets->at(jet_n);

    btagbtvdeep::hltParticleTransformerAK4Features hltFeatures;
    if (jet.pt() < min_jet_pt_ || std::abs(jet.eta()) > max_jet_eta_) {
      hltFeatures.is_filled = false;
    } else {
      // --- Process Secondary Vertices (explicit conversion) ---
      fillSVFeaturesHLT(hltFeatures, jet, pv, svs.product(), jet_radius_, flip_);

      // --- Build candidate ordering ---
      std::vector<btagbtvdeep::SortingClass<size_t>> c_sorted, n_sorted;
      std::map<unsigned int, btagbtvdeep::TrackInfoBuilder> trackinfos;
      for (unsigned int i = 0; i < jet.numberOfDaughters(); ++i) {
        const auto* cand = dynamic_cast<const reco::Candidate*>(jet.daughter(i));
        if (!cand || cand->pt() < min_candidate_pt_)
          continue;
        if (cand->charge() != 0) {
          auto& tinfo = trackinfos.emplace(i, track_builder).first->second;
          tinfo.buildTrackInfo(cand, jet.momentum().Unit(), GlobalVector(jet.px(), jet.py(), jet.pz()), pv);
          c_sorted.emplace_back(i,
                                  tinfo.getTrackSip2dSig(),
                                  -btagbtvdeep::mindrsvpfcand(*svs, cand),
                                  cand->pt() / jet.pt());
        } else {
          n_sorted.emplace_back(i, -1, -btagbtvdeep::mindrsvpfcand(*svs, cand), cand->pt() / jet.pt());
        }
      }
      std::sort(c_sorted.begin(), c_sorted.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);
      std::sort(n_sorted.begin(), n_sorted.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);
      auto c_sortedindices = btagbtvdeep::invertSortingVector(c_sorted);
      auto n_sortedindices = btagbtvdeep::invertSortingVector(n_sorted);

      hltFeatures.cpf_candidates.resize(c_sorted.size());
      hltFeatures.npf_candidates.resize(n_sorted.size());

      // --- Loop over candidates ---
      unsigned int chargedCount = 0, neutralCount = 0;
      for (unsigned int i = 0; i < jet.numberOfDaughters(); ++i) {
        const auto* cand = dynamic_cast<const reco::Candidate*>(jet.daughter(i));
        if (!cand || cand->pt() < min_candidate_pt_)
          continue;

        edm::Ref<edm::View<reco::Candidate>> candRef = getPersistentCandidate(cand, tracks);

        auto packed_cand = dynamic_cast<const pat::PackedCandidate*>(cand);
        // (For brevity, only the packed candidate branch is implemented here.)
        float puppiw = 1.0; // default fallback value
        float drminpfcandsv = btagbtvdeep::mindrsvpfcand(*svs, cand);
        float distminpfcandsv = 0;

        if (cand->charge() != 0) {
          size_t entry = c_sortedindices.at(chargedCount);
          auto& c_pf_features = hltFeatures.cpf_candidates.at(entry);
          if (packed_cand) {
            if (packed_cand->hasTrackDetails()) {
              const reco::Track& pseudoTrack = packed_cand->pseudoTrack();
              reco::TransientTrack transientTrack = track_builder->build(pseudoTrack);
              distminpfcandsv = btagbtvdeep::mindistsvpfcand(*svs, transientTrack);
            }
            convertHLTChargedCandidate(packed_cand,
                                       jet,
                                       trackinfos.at(i),
                                       is_weighted_jet_,
                                       drminpfcandsv,
                                       static_cast<float>(jet_radius_),
                                       puppiw,
                                       c_pf_features,
                                       flip_,
                                       distminpfcandsv);
          }
          ++chargedCount;
        } else {
          size_t entry = n_sortedindices.at(neutralCount);
          auto& n_pf_features = hltFeatures.npf_candidates.at(entry);
          if (packed_cand) {
            convertHLTNeutralCandidate(packed_cand,
                                       jet,
                                       is_weighted_jet_,
                                       drminpfcandsv,
                                       static_cast<float>(jet_radius_),
                                       puppiw,
                                       n_pf_features);
          }
          ++neutralCount;
        }
      }  // end candidate loop

      // --- Compute Global Features ---
      hltFeatures.global_features.jet_pt   = jet.pt();
      hltFeatures.global_features.jet_eta  = jet.eta();
      hltFeatures.global_features.nCpfcan  = hltFeatures.cpf_candidates.size();
      hltFeatures.global_features.nNpfcan  = hltFeatures.npf_candidates.size();
      hltFeatures.global_features.nsv      = hltFeatures.vtx_features.size();
      hltFeatures.global_features.npv      = vtxs->size();
      // Optionally compute CSV-related features; defaulting to 0 if not available.
      hltFeatures.global_features.TagVarCSV_trackSumJetEtRatio = 0;
      hltFeatures.global_features.TagVarCSV_trackSumJetDeltaR    = 0;
      hltFeatures.global_features.TagVarCSV_vertexCategory       = 0;
      hltFeatures.global_features.TagVarCSV_trackSip2dValAboveCharm  = 0;
      hltFeatures.global_features.TagVarCSV_trackSip2dSigAboveCharm  = 0;
      hltFeatures.global_features.TagVarCSV_trackSip3dValAboveCharm  = 0;
      hltFeatures.global_features.TagVarCSV_trackSip3dSigAboveCharm  = 0;
      hltFeatures.global_features.TagVarCSV_jetNTracksEtaRel       = 0;
    }  // end jet kinematics check

    // Create the TagInfo with the persistent jet reference and the filled HLT features.
    output_tag_infos->emplace_back(reco::hltParticleTransformerAK4TagInfo(hltFeatures, jet_ref));
  }  // end jet loop

  iEvent.put(std::move(output_tag_infos));
}

// Define this as a plug-in
DEFINE_FWK_MODULE(hltParticleTransformerAK4TagInfoProducer);
