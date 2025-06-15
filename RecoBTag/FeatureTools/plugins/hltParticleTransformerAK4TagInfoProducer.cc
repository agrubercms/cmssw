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
#include "DataFormats/JetReco/interface/GenJet.h"

#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>

#define DEBUG

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
      svfeat.jet_sv_chi2      = sv.vertexNormalizedChi2();
      const auto& dxy_meas = btagbtvdeep::vertexDxy(sv, pv);
      svfeat.jet_sv_dxy       = dxy_meas.value();
      svfeat.jet_sv_dxysig   = btagbtvdeep::catch_infs_and_bound(dxy_meas.value() / dxy_meas.error(), 0, -1, 800);
      const auto& d3d_meas = btagbtvdeep::vertexD3d(sv, pv);
      svfeat.jet_sv_d3d       = d3d_meas.value();
      svfeat.jet_sv_d3dsig   = btagbtvdeep::catch_infs_and_bound(d3d_meas.value() / d3d_meas.error(), 0, -1, 800);
      svfeat.jet_sv_pt_log = std::log(sv.pt());
      
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
  bool use_puppi_value_map_ = false;
  bool use_vertex_association_ = true;

  edm::EDGetTokenT<edm::ValueMap<float>> puppi_value_map_token_;
  edm::EDGetTokenT<edm::Association<VertexCollection>> vertex_associator_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> vertex_associator_quality_token_;

  // New GenJet token for input file "ak4GenJets" from process "HLT"
  const edm::EDGetTokenT<std::vector<reco::GenJet>> genjet_token_;
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
      vertex_associator_token_(consumes<edm::Association<VertexCollection>>(iConfig.getParameter<edm::InputTag>("vertex_associator"))),
      vertex_associator_quality_token_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("vertex_associator"))),
      // New token initialization for GenJets from input file "ak4GenJets" from process "HLT"
      genjet_token_(consumes<std::vector<reco::GenJet>>(edm::InputTag(iConfig.getParameter<edm::InputTag>("ak4GenJets"))))
{
  produces<hltParticleTransformerAK4TagInfoCollection>();
}

// --- fillDescriptions ---
void hltParticleTransformerAK4TagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<double>("jet_radius", 0.4);
  desc.add<double>("min_candidate_pt", 0.95);
  desc.add<bool>("flip", false);
  desc.add<edm::InputTag>("vertices", edm::InputTag("hltGoodOfflinePrimaryVertices"));
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
  // New parameter for GenJets input.
  desc.add<edm::InputTag>("ak4GenJets", edm::InputTag("ak4GenJets"));
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

#ifdef DEBUG
  if (tracks.isValid() && !tracks->empty()) {
    const reco::Candidate& cand = tracks->at(0);
    std::cout << "  Candidate type: " << typeid(cand).name() << std::endl;
  } else {
    std::cout << "  Candidate collection is invalid or empty." << std::endl;
  }
#endif

  edm::Handle<SVCollection> svs;
  iEvent.getByToken(sv_token_, svs);
#ifdef DEBUG
  std::cout << "DEBUG: Retrieved SV collection: valid=" << svs.isValid()
            << (svs.isValid() ? ", size=" + std::to_string(svs->size()) : "") << std::endl;
#endif

  // New: Retrieve handle for GenJets
  edm::Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genjet_token_, genJets);
  edm::Handle<edm::ValueMap<float>> puppi_value_map;
  if (use_puppi_value_map_) {
    iEvent.getByToken(puppi_value_map_token_, puppi_value_map);
  }

  edm::Handle<edm::ValueMap<int>> pvasq_value_map;
  edm::Handle<edm::Association<VertexCollection>> pvas;
  if (use_vertex_association_) {
    iEvent.getByToken(vertex_associator_quality_token_, pvasq_value_map);
    iEvent.getByToken(vertex_associator_token_, pvas);
  }

  edm::ESHandle<TransientTrackBuilder> track_builder = iSetup.getHandle(track_builder_token_);

#ifdef DEBUG
  std::cout << "=== Debug: Processing " << jets->size() << " jets ===" << std::endl;
#endif

  // Loop over jets
  for (std::size_t jet_n = 0; jet_n < jets->size(); ++jet_n) {
    edm::RefToBase<reco::Jet> jet_ref(jets, jet_n);
    const auto& jet = jets->at(jet_n);

#ifdef DEBUG
    std::cout << "Processing jet #" << jet_n << ": pt=" << jet.pt() << " eta=" << jet.eta() << " phi=" << jet.phi()
              << std::endl;
#endif

    btagbtvdeep::hltParticleTransformerAK4Features hltFeatures;
    if (jet.pt() < min_jet_pt_ || std::abs(jet.eta()) > max_jet_eta_) {
      hltFeatures.is_filled = false;
#ifdef DEBUG
      std::cout << "  Skipping jet (pt or eta out of range)" << std::endl;
#endif
    } else {
      hltFeatures.is_filled = true;

      // --- Process Secondary Vertices (explicit conversion) ---
      fillSVFeaturesHLT(hltFeatures, jet, pv, svs.product(), jet_radius_, flip_);

#ifdef DEBUG
      std::cout << "  Found " << hltFeatures.vtx_features.size() << " secondary vertices for this jet" << std::endl;
#endif

      // --- Collect and sort PF candidates by pt ---
      std::vector<const reco::PFCandidate*> chargedCandidates;
      std::vector<const reco::PFCandidate*> neutralCandidates;

      for (unsigned int i = 0; i < jet.numberOfDaughters(); ++i) {
        const auto* cand = dynamic_cast<const reco::PFCandidate*>(jet.daughter(i));
        if (!cand || cand->pt() < min_candidate_pt_)
          continue;
        if (cand->charge() != 0)
          chargedCandidates.push_back(cand);
        else
          neutralCandidates.push_back(cand);
      }

      auto sortByPt = [](const reco::PFCandidate* a, const reco::PFCandidate* b) { return a->pt() > b->pt(); };
      std::sort(chargedCandidates.begin(), chargedCandidates.end(), sortByPt);
      std::sort(neutralCandidates.begin(), neutralCandidates.end(), sortByPt);

      hltFeatures.cpf_candidates.reserve(chargedCandidates.size());
      hltFeatures.npf_candidates.reserve(neutralCandidates.size());

      for (const auto* cand : chargedCandidates) {
        float puppiw = 1.0;
        float drminpfcandsv = btagbtvdeep::mindrsvpfcand(*svs, cand);
        if (cand->trackRef().isNonnull()) {
          reco::TransientTrack transientTrack = track_builder->build(*(cand->trackRef()));
        }

        btagbtvdeep::TrackInfoBuilder trackInfo(track_builder);
        // Note: The offline version builds TrackInfo using the jet_dir, jet_ref_track_dir, and pv.
        // Your HLT version uses jet.momentum().Unit(), GlobalVector(jet.px(), jet.py(), jet.pz()), and pv. This seems consistent.
        
        int pv_ass_quality = 0; // Default quality
        reco::VertexRef pv_ass = reco::VertexRef(vtxs, 0); // Default to the leading primary vertex from the main vertex collection

        if (use_vertex_association_) {
          // Get a persistent edm::Ref to the candidate in the original 'tracks' collection
          edm::Ref<edm::View<reco::Candidate>> candRef = getPersistentCandidate(cand, tracks);

          if (candRef.isNonnull()) {
            // Ensure the handles are valid before attempting to access their data
            // (getByToken would throw if the product is not found, but an extra check is safe)
            if (pvas.isValid() && pvasq_value_map.isValid()) {
              // Get the associated vertex using the 'pvas' handle
              const reco::VertexRef& pv_orig = (*pvas)[candRef]; 
              if (pv_orig.isNonnull()) {
                pv_ass = pv_orig; // Update pv_ass to the actual associated vertex
              }
              // Get the association quality
              pv_ass_quality = (*pvasq_value_map)[candRef];
            } else {
              // Optional: Log a warning if maps are expected but not valid
              std::cout << "Warning: Vertex association maps are not valid." << std::endl;
            }
          } else {
            // Optional: Log a warning if candidate not found in original collection
            std::cout << "Warning: Candidate not found in original track collection for PV association." << std::endl;
          }
        }
        
        // Pass the determined pv_ass (dereferenced) to TrackInfoBuilder
        // pv_ass will be the associated PV if found and use_vertex_association_ is true,
        // otherwise it defaults to the leading PV. This is safe because vtxs is checked for non-emptiness earlier.
        trackInfo.buildTrackInfo(cand, jet.momentum().Unit(), GlobalVector(jet.px(), jet.py(), jet.pz()), *pv_ass);

        hltCpfCandidateFeatures feat;
        feat.jet_pfcand_deta = jet.eta() - cand->eta();
        feat.jet_pfcand_dphi = reco::deltaPhi(jet.phi(), cand->phi());
        feat.jet_pfcand_pt_log = (cand->pt() > 0) ? std::log(cand->pt()) : 0;
        feat.jet_pfcand_energy_log = (cand->energy() > 0) ? std::log(cand->energy()) : 0;
        feat.jet_pfcand_charge = static_cast<float>(cand->charge());
        feat.jet_pfcand_frompv = static_cast<float>(pv_ass_quality); // Use the obtained quality
        feat.jet_pfcand_nlostinnerhits = btagbtvdeep::lost_inner_hits_from_pfcand(*cand);
        
        if (cand->bestTrack()) {
          const auto& track = *(cand->bestTrack());
          feat.jet_pfcand_track_chi2 = track.normalizedChi2();
          feat.jet_pfcand_track_qual = static_cast<float>(track.qualityMask());
          // Use the potentially updated pv_ass for dz/dxy calculations
          feat.jet_pfcand_dz = track.dz(pv_ass->position()); 
          feat.jet_pfcand_dzsig =
              fabs(btagbtvdeep::catch_infs_and_bound(track.dz(pv_ass->position()) / track.dzError(), 300, -1, 300));
          feat.jet_pfcand_dxy = track.dxy(pv_ass->position());
          feat.jet_pfcand_dxysig =
              fabs(btagbtvdeep::catch_infs_and_bound(track.dxy(pv_ass->position()) / track.dxyError(), 300, -1, 300));
          feat.jet_pfcand_npixhits = track.hitPattern().numberOfValidPixelHits();
          feat.jet_pfcand_nstriphits = track.hitPattern().stripLayersWithMeasurement();
        } else {
          feat.jet_pfcand_track_chi2 = -1;
          feat.jet_pfcand_track_qual = 0;
          feat.jet_pfcand_dz = 0;
          feat.jet_pfcand_dzsig = 0;
          feat.jet_pfcand_dxy = 0;
          feat.jet_pfcand_dxysig = 0;
          feat.jet_pfcand_npixhits = 0;
          feat.jet_pfcand_nstriphits = 0;
        }
        feat.jet_pfcand_etarel = trackInfo.getTrackEtaRel();
        feat.jet_pfcand_pperp_ratio = trackInfo.getTrackPtRatio(); 
        feat.jet_pfcand_ppara_ratio = trackInfo.getTrackPParRatio();
        feat.jet_pfcand_trackjet_d3d = trackInfo.getTrackSip3dVal();
        feat.jet_pfcand_trackjet_d3dsig = trackInfo.getTrackSip3dSig();
        feat.jet_pfcand_trackjet_dist = -trackInfo.getTrackJetDistVal();
        feat.jet_pfcand_trackjet_decayL = trackInfo.getTrackJetDecayLen();
        feat.jet_pfcand_pt = cand->pt();
        feat.jet_pfcand_eta = cand->eta();
        feat.jet_pfcand_phi = cand->phi();
        feat.jet_pfcand_energy = cand->energy();

        hltFeatures.cpf_candidates.push_back(feat);
      }

      // --- Compute Global Features ---
      hltFeatures.global_features.jet_pt = jet.pt();
      hltFeatures.global_features.jet_eta = jet.eta();
      hltFeatures.global_features.jet_phi = jet.phi();
      hltFeatures.global_features.jet_energy = jet.energy();
    }  // end jet kinematics check

    // Create the TagInfo with the persistent jet reference and the filled HLT features.
    output_tag_infos->emplace_back(reco::hltParticleTransformerAK4TagInfo(hltFeatures, jet_ref));
  }  // end jet loop

  // --- Debug printout for all tag infos individually ---
#ifdef DEBUG
  std::cout << "=== Debug: Individually Printed Tag Infos ===" << std::endl;
  for (std::size_t i = 0; i < output_tag_infos->size(); ++i) {
    const auto& tagInfo = output_tag_infos->at(i);
    const auto& features = tagInfo.features();
    
    std::cout << "Tag Info " << i << ":\n";
    
    // Print global features
    std::cout << "  -- Global Features --\n";
    std::cout << "    jet_pt: " << features.global_features.jet_pt << "\n";
    std::cout << "    jet_eta: " << features.global_features.jet_eta << "\n";
    
    // Print secondary vertex features
    std::cout << "  -- Secondary Vertices (" << features.vtx_features.size() << ") --\n";
    for (size_t j = 0; j < features.vtx_features.size(); ++j) {
      const auto& sv = features.vtx_features[j];
      std::cout << "    SV #" << j << ":\n";
      std::cout << "      jet_sv_ntrack: " << sv.jet_sv_ntrack << "\n";
      std::cout << "      jet_sv_mass: " << sv.jet_sv_mass << "\n";
      std::cout << "      jet_sv_pt_log: " << sv.jet_sv_pt_log << "\n";
      std::cout << "      jet_sv_deta: " << sv.jet_sv_deta << "\n";
      std::cout << "      jet_sv_dphi: " << sv.jet_sv_dphi << "\n";
      std::cout << "      jet_sv_chi2: " << sv.jet_sv_chi2 << "\n";
      std::cout << "      jet_sv_dxy: " << sv.jet_sv_dxy << "\n";
      std::cout << "      jet_sv_dxysig: " << sv.jet_sv_dxysig << "\n";
      std::cout << "      jet_sv_d3d: " << sv.jet_sv_d3d << "\n";
      std::cout << "      jet_sv_d3dsig: " << sv.jet_sv_d3dsig << "\n";
      std::cout << "      jet_sv_pt: " << sv.jet_sv_pt << "\n";
      std::cout << "      jet_sv_eta: " << sv.jet_sv_eta << "\n";
      std::cout << "      jet_sv_phi: " << sv.jet_sv_phi << "\n";
      std::cout << "      jet_sv_energy: " << sv.jet_sv_energy << "\n";
    }
    
    std::cout << std::endl;
  }
#endif

  iEvent.put(std::move(output_tag_infos));
}

// Define this as a plug-in
DEFINE_FWK_MODULE(hltParticleTransformerAK4TagInfoProducer);
