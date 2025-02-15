// hltParticleTransformerAK4TagInfoProducer.cc

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
#include "RecoBTag/FeatureTools/interface/SecondaryVertexConverter.h"

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
  // Typedef for the output collection.
  typedef std::vector<reco::hltParticleTransformerAK4TagInfo> hltParticleTransformerAK4TagInfoCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
  typedef reco::VertexCollection VertexCollection;

  void beginStream(edm::StreamID) override {}
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override {}

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

  // --- Local helper conversion methods ---
  // Minimal conversion for charged candidates (packed version)
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
    // Compute a weighted candidate ptrel.
    float constituentWeight = isWeightedJet ? puppiw : 1.0f;
    feat.Cpfcan_ptrel = (cand->pt() * constituentWeight) / jet.pt();
    // (We omit erel since the hlt structure does not have Cpfcan_erel.)
    feat.Cpfcan_drminsv = btagbtvdeep::catch_infs_and_bound(drminpfcandsv, 0, -jetR, 0, -jetR);
    feat.Cpfcan_pt = cand->pt();
    feat.Cpfcan_eta = cand->eta();
    feat.Cpfcan_phi = cand->phi();
    feat.Cpfcan_e   = cand->energy();
    
    // HLT-specific fields.
    feat.Cpfcan_VTX_ass = cand->pvAssociationQuality();
    feat.Cpfcan_puppiw  = puppiw;
    // Additional fields (track chi2, quality, etc.) can be added if needed.
  }

  // Minimal conversion for neutral candidates (packed version)
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
    
    // HLT-specific: puppi weight and hadronic fraction.
    feat.Npfcan_puppiw = puppiw;
    feat.Npfcan_HadFrac = cand->hcalFraction();
  }

  // --- Configuration parameters ---
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

  // Extra parameters.
  const bool fallback_puppi_weight_;
  const bool fallback_vertex_association_;
  const double max_sip3dsig_for_flip_;

  const edm::EDGetTokenT<edm::ValueMap<float>> puppi_value_map_token_;
  const edm::EDGetTokenT<edm::Association<VertexCollection>> vertex_associator_token_;
};

// Constructor
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

// fillDescriptions
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

// produce method
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

  // std::cout << "[DEBUG] Starting hltParticleTransformerAK4TagInfoProducer::produce" << std::endl;
  // std::cout << "[DEBUG] Number of jets: " << jets->size() << std::endl;

  // Loop over jets
  for (std::size_t jet_n = 0; jet_n < jets->size(); ++jet_n) {
    edm::RefToBase<reco::Jet> jet_ref(jets, jet_n);
    const auto& jet = jets->at(jet_n);

    // std::cout << "[DEBUG] Processing jet " << jet_n << " with pt: " << jet.pt() << ", eta: " << jet.eta() << std::endl;

    // Initialize HLT feature container
    btagbtvdeep::hltParticleTransformerAK4Features hltFeatures;
    if (jet.pt() < min_jet_pt_ || std::abs(jet.eta()) > max_jet_eta_) {
      hltFeatures.is_filled = false;
    }

    // --- Process Secondary Vertices ---
    auto svs_sorted = *svs;
    std::sort(svs_sorted.begin(), svs_sorted.end(), [&pv](const auto& sv1, const auto& sv2) {
      return btagbtvdeep::sv_vertex_comparator(sv1, sv2, pv);
    });
    for (const auto& sv : svs_sorted) {
      if (reco::deltaR2(sv, jet) > (jet_radius_ * jet_radius_))
        continue;
      hltFeatures.vtx_features.emplace_back();
      // Use the HLT overload for SV conversion.
      btagbtvdeep::svToFeatures(sv, pv, jet, hltFeatures.vtx_features.back(), flip_);
    }

    // --- Build candidate ordering ---
    std::vector<btagbtvdeep::SortingClass<size_t>> c_sorted, n_sorted;
    std::map<unsigned int, btagbtvdeep::TrackInfoBuilder> trackinfos;
    for (unsigned int i = 0; i < jet.numberOfDaughters(); ++i) {
      const auto* cand = dynamic_cast<const reco::Candidate*>(jet.daughter(i));
      if (!cand || cand->pt() < min_candidate_pt_)
        continue;
      if (cand->charge() != 0) {
        auto& tinfo = trackinfos.emplace(i, track_builder).first->second;
        tinfo.buildTrackInfo(cand, jet.momentum().Unit(),
                             GlobalVector(jet.px(), jet.py(), jet.pz()), pv);
        c_sorted.emplace_back(i,
                                tinfo.getTrackSip2dSig(),
                                -btagbtvdeep::mindrsvpfcand(*svs, cand),
                                cand->pt() / jet.pt());
      } else {
        n_sorted.emplace_back(i, -1,
                              -btagbtvdeep::mindrsvpfcand(*svs, cand),
                              cand->pt() / jet.pt());
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
      auto reco_cand   = dynamic_cast<const reco::PFCandidate*>(cand);

      float puppiw = 1.0; // fallback value
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
          // Use our local helper conversion function for charged candidates.
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
        // (A branch for reco_cand can be added if needed.)
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
        // (A branch for reco_cand can be added if needed.)
        ++neutralCount;
      }
    }  // end candidate loop

    // Create the TagInfo with the persistent jet reference and the HLT features.
    output_tag_infos->emplace_back(reco::hltParticleTransformerAK4TagInfo(hltFeatures, jet_ref));
  }  // end jet loop

  // std::cout << "[DEBUG] Finished jet loop" << std::endl;
  iEvent.put(std::move(output_tag_infos));
  // std::cout << "[DEBUG] Finished put" << std::endl;
}

// Define this as a plug-in
DEFINE_FWK_MODULE(hltParticleTransformerAK4TagInfoProducer);
