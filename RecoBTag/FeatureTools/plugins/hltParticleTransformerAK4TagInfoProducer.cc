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
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "DataFormats/Math/interface/Vector3D.h" // For GlobalVector

#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>

#include "TVector3.h"

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
    
    GlobalVector jet_dir(jet.px(), jet.py(), jet.pz()); // Jet direction for signing

    // Loop over sorted SVs and fill features for those within the jet radius.
    for (const auto& sv_cand : svs_sorted) { // Renamed sv to sv_cand to avoid conflict with reco_sv
      if (reco::deltaR2(sv_cand, jet) > (jetR * jetR))
        continue;
      
      // Use the HLT-specific vertex features type.
      hltVtxFeatures svfeat;
      // Map available quantities.
      svfeat.jet_sv_pt        = sv_cand.pt();
      svfeat.jet_sv_deta      = sv_cand.eta() - jet.eta();
      svfeat.jet_sv_dphi      = sv_cand.phi() - jet.phi();  // raw subtraction like DeepJetNTupler
      svfeat.jet_sv_eta       = sv_cand.eta();
      svfeat.jet_sv_phi       = sv_cand.phi();
      svfeat.jet_sv_energy    = sv_cand.energy();
      svfeat.jet_sv_mass      = sv_cand.mass();
      svfeat.jet_sv_ntrack    = sv_cand.numberOfDaughters();
      svfeat.jet_sv_chi2      = sv_cand.vertexNormalizedChi2();

      reco::Vertex::CovarianceMatrix csv;
      sv_cand.fillVertexCovariance(csv);
      reco::Vertex svtx(sv_cand.vertex(), csv);

      GlobalVector jet_vec(jet.px(), jet.py(), jet.pz());

      VertexDistanceXY dxy;
      auto dxy_meas = dxy.signedDistance(svtx, pv, jet_vec);
      svfeat.jet_sv_dxy       = dxy_meas.value();
      svfeat.jet_sv_dxysig    = std::fabs(dxy_meas.significance());

      VertexDistance3D d3d;
      auto d3d_meas = d3d.signedDistance(svtx, pv, jet_vec);
      svfeat.jet_sv_d3d       = d3d_meas.value();
      svfeat.jet_sv_d3dsig    = std::fabs(d3d_meas.significance());
      svfeat.jet_sv_pt_log = std::log(sv_cand.pt());

      // costhetasvpv: cosine between SV flight direction (PV->SV) and SV momentum,
      // with optional flip, identical to SecondaryVertexConverter / training
      const float cos_sv_pv = btagbtvdeep::vertexDdotP(sv_cand, pv);
      svfeat.jet_sv_costhetasvpv = (flip ? -1.f : 1.f) * cos_sv_pv;

      // SV energy ratio w.r.t. jet
      svfeat.jet_sv_enratio = (jet.energy() > 0.f ? sv_cand.energy() / jet.energy() : 0.f);

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
  bool use_vertex_association_ = false;  // Will be set dynamically based on fallback_vertex_association_

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

  std::unique_ptr<reco::VertexRefProd> PVRefProd = std::make_unique<reco::VertexRefProd>(vtxs);
  const float min_track_pt_property = 0.5f;
  const int min_valid_pixel_hits = 0;
  const int covarianceVersion = 0;
  const std::vector<int> covariancePackingSchemas = {8, 264, 520, 776, 0};

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
  // Also try to load PUPPI map even if not explicitly enabled, for better matching with ntupler
  if (!use_puppi_value_map_) {
    iEvent.getByToken(puppi_value_map_token_, puppi_value_map);
    use_puppi_value_map_ = puppi_value_map.isValid();
  }

  edm::Handle<edm::ValueMap<int>> pvasq_value_map;
  edm::Handle<edm::Association<VertexCollection>> pvas;
  // Try to use vertex association if not in fallback mode
  use_vertex_association_ = !fallback_vertex_association_;
  if (use_vertex_association_) {
    iEvent.getByToken(vertex_associator_quality_token_, pvasq_value_map);
    iEvent.getByToken(vertex_associator_token_, pvas);
    // If maps aren't valid, switch to fallback mode
    if (!pvasq_value_map.isValid() || !pvas.isValid()) {
      use_vertex_association_ = false;
    }
  }
#ifdef DEBUG
  std::cout << "DEBUG: fallback_vertex_association=" << fallback_vertex_association_
            << " use_vertex_association=" << use_vertex_association_
            << " pvasq valid=" << pvasq_value_map.isValid()
            << " pvas valid=" << pvas.isValid()
            << " puppi_map valid=" << puppi_value_map.isValid()
            << " candidates collection size=" << (tracks.isValid() ? (int)tracks->size() : -1)
            << std::endl;
#endif

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
      std::vector<const reco::PFCandidate*> pfCandidates;

      for (unsigned int i = 0; i < jet.numberOfDaughters(); ++i) {
        const auto* cand = dynamic_cast<const reco::PFCandidate*>(jet.daughter(i));
        if (!cand) continue;
        
        // Apply minimum pT cut (match ntupler's behavior)
        if (cand->pt() < min_candidate_pt_) {
#ifdef DEBUG
          std::cout << "    Skipping cand #" << i << ": pt=" << cand->pt() << " < " << min_candidate_pt_ 
                    << " pdgId=" << cand->pdgId() << std::endl;
#endif
          continue;
        }
        
        // Apply PUPPI weight filtering if available
        bool include_cand = true;
        if (use_puppi_value_map_ && puppi_value_map.isValid()) {
          edm::Ref<edm::View<reco::Candidate>> candRef = getPersistentCandidate(cand, tracks);
          if (candRef.isNonnull()) {
            float puppiw = (*puppi_value_map)[candRef];
            // Skip candidates with too-low PUPPI weights (match ntupler's min_puppi_wgt_ = -1.0)
            const double min_puppi_wgt = -1.0;
            if (puppiw < min_puppi_wgt) {
              include_cand = false;
#ifdef DEBUG
              std::cout << "    Skipping cand #" << i << ": puppiw=" << puppiw << " < " << min_puppi_wgt 
                        << " pt=" << cand->pt() << " pdgId=" << cand->pdgId() << std::endl;
#endif
            }
          }
        }
        
        if (include_cand) {
          pfCandidates.push_back(cand);
#ifdef DEBUG
          std::cout << "    Including cand #" << i << ": pt=" << cand->pt() << " pdgId=" << cand->pdgId() << std::endl;
#endif
        }
      }
      
#ifdef DEBUG
      std::cout << "  Total PF candidates after filtering: " << pfCandidates.size() << std::endl;
#endif

      auto sortByPt = [](const reco::PFCandidate* a, const reco::PFCandidate* b) { return a->pt() > b->pt(); };
      std::sort(pfCandidates.begin(), pfCandidates.end(), sortByPt);

      hltFeatures.cpf_candidates.reserve(pfCandidates.size());

      for (const auto* cand : pfCandidates) {
        float puppiw = 1.0;
        // float drminpfcandsv = btagbtvdeep::mindrsvpfcand(*svs, cand);  // no corresponding HLT feature, drop to avoid warning
        if (cand->trackRef().isNonnull()) {
          reco::TransientTrack transientTrack = track_builder->build(*(cand->trackRef()));
        }
#ifdef DEBUG
        std::cout << "  Processing candidate: pt=" << cand->pt() << " pdgId=" << cand->pdgId() 
                  << " eta=" << cand->eta() << " phi=" << cand->phi() << std::endl;
#endif
        // optionally get Puppi weight from value map
        if (use_puppi_value_map_) {
          edm::Ref<edm::View<reco::Candidate>> candRef = getPersistentCandidate(cand, tracks);
          if (candRef.isNonnull()) {
            if (puppi_value_map.isValid()) {
              puppiw = (*puppi_value_map)[candRef];
            } else if (!fallback_puppi_weight_) {
              puppiw = 0.f;
            }
          }
        }

        btagbtvdeep::TrackInfoBuilder trackInfo(track_builder);
        // Note: The offline version builds TrackInfo using the jet_dir, jet_ref_track_dir, and pv.
        // Your HLT version uses jet.momentum().Unit(), GlobalVector(jet.px(), jet.py(), jet.pz()), and pv. This seems consistent.
        
        int pv_ass_quality = 0; // Default quality
        reco::VertexRef pv_ass = reco::VertexRef(vtxs, 0); // Default to the leading primary vertex from the main vertex collection

        if (use_vertex_association_) {
          // Get a persistent edm::Ref to the candidate in the original 'tracks' collection
          edm::Ref<edm::View<reco::Candidate>> candRef = getPersistentCandidate(cand, tracks);
#ifdef DEBUG
          std::cout << "  DEBUG cand pt=" << cand->pt() << " pdgId=" << cand->pdgId()
                    << " candRef valid=" << candRef.isNonnull();
          if (candRef.isNonnull()) std::cout << " candRef.key=" << candRef.key();
          std::cout << std::endl;
#endif
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
#ifdef DEBUG
              std::cout << "    -> pv_ass_quality=" << pv_ass_quality
                        << " pv_ass.key=" << pv_ass.key() << std::endl;
#endif
            } else {
              std::cout << "Warning: Vertex association maps are not valid." << std::endl;
            }
          } else {
            std::cout << "Warning: Candidate not found in original track collection for PV association." << std::endl;
          }
        }
        
        // Fallback: if no vertex association was found and candidate has a track, find closest PV by dz
        if (!use_vertex_association_ && cand->bestTrack() && pv_ass.key() == 0) {
          const auto& track = *(cand->bestTrack());
          float z_dist = 99999;
          int pv_pos = 0;
          for (size_t iv = 0; iv < vtxs->size(); iv++) {
            float dz = std::abs(track.dz((*vtxs)[iv].position()));
            if (dz < z_dist) {
              z_dist = dz;
              pv_pos = iv;
            }
          }
          pv_ass = reco::VertexRef(vtxs, pv_pos);
        }
        
        const reco::Track* track = cand->bestTrack();

        // Build a PackedCandidate to match the Ntupler/DeepBoostedJetTagInfoProducer logic.
        pat::PackedCandidate packed_candidate;
        math::XYZPoint pv_ass_pos;
        if (not pv_ass.isNonnull()) {
          if (track) {
            float z_dist = 99999;
            int pv_pos = -1;
            for (size_t iv = 0; iv < vtxs->size(); iv++) {
              float dz = std::abs(track->dz(((*vtxs)[iv]).position()));
              if (dz < z_dist) {
                z_dist = dz;
                pv_pos = iv;
              }
            }
            pv_ass = reco::VertexRef(vtxs, pv_pos);
          } else {
            pv_ass = reco::VertexRef(vtxs, 0);
          }
        }
        pv_ass_pos = pv_ass->position();

        if (track) {
          packed_candidate = pat::PackedCandidate(cand->polarP4(),
                                                  track->referencePoint(),
                                                  track->pt(),
                                                  track->eta(),
                                                  track->phi(),
                                                  cand->pdgId(),
                                                  (*PVRefProd),
                                                  pv_ass.key());
          packed_candidate.setAssociationQuality(pat::PackedCandidate::PVAssociationQuality(
              btagbtvdeep::vtx_ass_from_pfcand(*cand, pv_ass_quality, pv_ass)));
          packed_candidate.setCovarianceVersion(covarianceVersion);

          pat::PackedCandidate::LostInnerHits lostHits = pat::PackedCandidate::noLostInnerHits;
          int nlost = track->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
          if (nlost == 0) {
            if (track->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1))
              lostHits = pat::PackedCandidate::validHitInFirstPixelBarrelLayer;
          } else {
            lostHits = (nlost == 1 ? pat::PackedCandidate::oneLostInnerHit : pat::PackedCandidate::moreLostInnerHits);
          }
          packed_candidate.setLostInnerHits(lostHits);
          packed_candidate.setTrkAlgo(static_cast<uint8_t>(track->algo()), static_cast<uint8_t>(track->originalAlgo()));

          const bool use_track_properties = track->pt() > min_track_pt_property;
          if (use_track_properties) {
            packed_candidate.setFirstHit(track->hitPattern().getHitPattern(reco::HitPattern::TRACK_HITS, 0));
            if (std::abs(cand->pdgId()) == 22) {
              packed_candidate.setTrackProperties(*track, covariancePackingSchemas[4], covarianceVersion);
            } else if (track->hitPattern().numberOfValidPixelHits() > min_valid_pixel_hits) {
              packed_candidate.setTrackProperties(*track, covariancePackingSchemas[0], covarianceVersion);
            } else {
              packed_candidate.setTrackProperties(*track, covariancePackingSchemas[1], covarianceVersion);
            }
          } else if (packed_candidate.pt() > min_track_pt_property) {
            if (track->hitPattern().numberOfValidPixelHits() > 0) {
              packed_candidate.setTrackProperties(*track, covariancePackingSchemas[2], covarianceVersion);
            } else {
              packed_candidate.setTrackProperties(*track, covariancePackingSchemas[3], covarianceVersion);
            }
          }
          packed_candidate.setTrackHighPurity(cand->trackRef().isNonnull() &&
                                              cand->trackRef()->quality(reco::Track::highPurity));
        } else {
          packed_candidate = pat::PackedCandidate(cand->polarP4(),
                                                  pv_ass_pos,
                                                  cand->pt(),
                                                  cand->eta(),
                                                  cand->phi(),
                                                  cand->pdgId(),
                                                  (*PVRefProd),
                                                  pv_ass.key());
          packed_candidate.setAssociationQuality(
              pat::PackedCandidate::PVAssociationQuality(pat::PackedCandidate::UsedInFitTight));
        }

        track = packed_candidate.bestTrack();

        // Pass the determined pv_ass (dereferenced) to TrackInfoBuilder, using the packed candidate.
        trackInfo.buildTrackInfo(&packed_candidate,
                                 jet.momentum().Unit(),
                                 GlobalVector(jet.px(), jet.py(), jet.pz()),
                                 *pv_ass);

#ifdef DEBUG
        if (packed_candidate.bestTrack()) {
          std::cout << "    PV association: pv_ass.key=" << pv_ass.key() 
                    << " pv_ass_quality=" << pv_ass_quality
                    << " fromPV()=" << packed_candidate.fromPV() << std::endl;
          std::cout << "    Track: pt=" << packed_candidate.bestTrack()->pt()
                    << " original_track.pt=" << cand->bestTrack()->pt() << std::endl;
          std::cout << "    dz(pos=" << pv_ass_pos.z() << ")=" << packed_candidate.dz(pv_ass_pos)
                    << " dzError=" << packed_candidate.dzError()
                    << " dz/dzError=" << (packed_candidate.dzError() > 0 ? packed_candidate.dz(pv_ass_pos) / packed_candidate.dzError() : 0)
                    << std::endl;
          std::cout << "    dxy=" << packed_candidate.dxy(pv_ass_pos)
                    << " dxyError=" << packed_candidate.dxyError()
                    << " dxy/dxyError=" << (packed_candidate.dxyError() > 0 ? packed_candidate.dxy(pv_ass_pos) / packed_candidate.dxyError() : 0)
                    << std::endl;
        }
#endif

        hltCpfCandidateFeatures feat;
        feat.jet_pfcand_deta = jet.eta() - cand->eta();
        feat.jet_pfcand_dphi = reco::deltaPhi(jet.phi(), cand->phi());
        feat.jet_pfcand_pt_log = (cand->pt() > 0) ? std::log(cand->pt()) : 0;
        feat.jet_pfcand_energy_log = (cand->energy() > 0) ? std::log(cand->energy()) : 0;
        feat.jet_pfcand_charge = static_cast<float>(cand->charge());

        feat.jet_pfcand_frompv = static_cast<float>(packed_candidate.fromPV());
        feat.jet_pfcand_nlostinnerhits = packed_candidate.lostInnerHits();

        bool highPurity = false;

        if (track) {
          feat.jet_pfcand_track_chi2 = track->normalizedChi2();
          feat.jet_pfcand_track_qual = track->qualityMask();
          feat.jet_pfcand_dz = packed_candidate.dz(pv_ass_pos);
          feat.jet_pfcand_dzsig = fabs(packed_candidate.dz(pv_ass_pos) / packed_candidate.dzError());
          feat.jet_pfcand_dxy = packed_candidate.dxy(pv_ass_pos);
          feat.jet_pfcand_dxysig = fabs(packed_candidate.dxy(pv_ass_pos) / packed_candidate.dxyError());
          feat.jet_pfcand_npixhits = packed_candidate.numberOfPixelHits();
          feat.jet_pfcand_nstriphits = packed_candidate.stripLayersWithMeasurement();
          highPurity = packed_candidate.trackHighPurity();
        } else {
          feat.jet_pfcand_track_chi2 = 0;
          feat.jet_pfcand_track_qual = 0;
          feat.jet_pfcand_dz = 0;
          feat.jet_pfcand_dzsig = 0;
          feat.jet_pfcand_dxy = 0;
          feat.jet_pfcand_dxysig = 0;
          feat.jet_pfcand_npixhits = 0;
          feat.jet_pfcand_nstriphits = 0;
        }

        feat.jet_pfcand_etarel = trackInfo.getTrackEtaRel();
        // Compute pperp_ratio and ppara_ratio like DeepJetNTupler:
        // jet_direction.Perp(cand_direction) / cand_direction.Mag()
        // jet_direction.Dot(cand_direction) / cand_direction.Mag()
        TVector3 jet_direction(jet.px(), jet.py(), jet.pz());
        jet_direction = jet_direction.Unit();
        TVector3 cand_direction(cand->px(), cand->py(), cand->pz());
        float cand_mag = cand_direction.Mag();
        feat.jet_pfcand_pperp_ratio = (cand_mag > 0) ? jet_direction.Perp(cand_direction) / cand_mag : 0;
        feat.jet_pfcand_ppara_ratio = (cand_mag > 0) ? jet_direction.Dot(cand_direction) / cand_mag : 0;
        feat.jet_pfcand_trackjet_d3d = trackInfo.getTrackSip3dVal();
        feat.jet_pfcand_trackjet_d3dsig = trackInfo.getTrackSip3dSig();
        feat.jet_pfcand_trackjet_dist = -trackInfo.getTrackJetDistVal();
        feat.jet_pfcand_trackjet_decayL = trackInfo.getTrackJetDecayLen();

        // calorimeter fractions: set to 0 to match DeepJetNTupler behavior
        // The Ntupler creates a PackedCandidate but doesn't call setCaloFraction(),
        // so caloFraction() and hcalFraction() return 0
        feat.jet_pfcand_calofraction = 0.f;
        feat.jet_pfcand_hcalfraction = 0.f;

        // new: Puppi weight, track high-purity flag, and particle ID
        feat.jet_pfcand_puppiw = puppiw;
        feat.jet_pfcand_highpurity = highPurity ? 1.f : 0.f;

        // Use abs(pdgId) to match DeepJetNTupler (e.g., 211 for pion, 13 for muon)
        feat.jet_pfcand_id = static_cast<float>(std::abs(cand->pdgId()));

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
  //=== Dump exactly the four ONNX inputs: global, cpf, npf, sv ===
  constexpr size_t kGlobalFeatures = 0;
  constexpr size_t kCpfCandidates   = 1;
  constexpr size_t kNpfCandidates   = 2;
  constexpr size_t kVtxFeatures     = 3;
  const int n_features_cpf = 31;  // Updated to match YAML: 31 features per charged candidate
  const int n_features_npf = 0;   // we don't store neutral PF in this tagger
  const int n_features_sv  = 16;  // Updated to match YAML: 16 features per SV

  for (std::size_t jet_idx = 0; jet_idx < output_tag_infos->size(); ++jet_idx) {
    const auto& tagInfo = output_tag_infos->at(jet_idx);
    const auto& features = tagInfo.features();
    if (!features.is_filled) continue;

    std::cout << "=== Dumping TagInfoProducer Features ===" << std::endl;

    //--- Global features (size = 4) ---
    std::cout << "  -- Global Features (data_["<<kGlobalFeatures<<"], size=4) --" << std::endl;
    std::cout << "    data_["<<kGlobalFeatures<<"][0]: " << features.global_features.jet_pt     << std::endl;
    std::cout << "    data_["<<kGlobalFeatures<<"][1]: " << features.global_features.jet_eta    << std::endl;
    std::cout << "    data_["<<kGlobalFeatures<<"][2]: " << features.global_features.jet_phi    << std::endl;
    std::cout << "    data_["<<kGlobalFeatures<<"][3]: " << features.global_features.jet_energy << std::endl;
    std::cout << "    ncands: " << features.cpf_candidates.size() << std::endl;

    //--- Charged PF candidates (each has 31 features) ---
    size_t nCpf = features.cpf_candidates.size();
    std::cout << "  -- Charged PF Candidates (data_["<<kCpfCandidates<<"], size=" 
              << (nCpf * n_features_cpf) << ") --" << std::endl;
    for (size_t i = 0; i < nCpf; ++i) {
      const auto& cpf = features.cpf_candidates[i];
      const float feats[n_features_cpf] = {
          cpf.jet_pfcand_deta,
          cpf.jet_pfcand_dphi,
          cpf.jet_pfcand_pt_log,
          cpf.jet_pfcand_energy_log,
          cpf.jet_pfcand_charge,
          cpf.jet_pfcand_frompv,
          static_cast<float>(cpf.jet_pfcand_nlostinnerhits),
          cpf.jet_pfcand_track_chi2,
          cpf.jet_pfcand_track_qual,
          cpf.jet_pfcand_dz,
          cpf.jet_pfcand_dzsig,
          cpf.jet_pfcand_dxy,
          cpf.jet_pfcand_dxysig,
          cpf.jet_pfcand_etarel,
          cpf.jet_pfcand_pperp_ratio,
          cpf.jet_pfcand_ppara_ratio,
          cpf.jet_pfcand_trackjet_d3d,
          cpf.jet_pfcand_trackjet_d3dsig,
          cpf.jet_pfcand_trackjet_dist,
          cpf.jet_pfcand_trackjet_decayL,
          static_cast<float>(cpf.jet_pfcand_npixhits),
          static_cast<float>(cpf.jet_pfcand_nstriphits),
          cpf.jet_pfcand_calofraction,
          cpf.jet_pfcand_hcalfraction,
          cpf.jet_pfcand_puppiw,
          cpf.jet_pfcand_highpurity,
          cpf.jet_pfcand_id,
          cpf.jet_pfcand_pt,
          cpf.jet_pfcand_eta,
          cpf.jet_pfcand_phi,
          cpf.jet_pfcand_energy
      };
      for (int j = 0; j < n_features_cpf; ++j) {
        std::cout << "    data_["<<kCpfCandidates<<"]["<<(i*n_features_cpf + j)<<"]: "
                  << feats[j] << std::endl;
      }
    }

    //--- Neutral PF candidates: none for this producer ---
    std::cout << "  -- Neutral PF Candidates (data_["<<kNpfCandidates<<"], size=0) --" << std::endl;

    //--- Secondary-vertex candidates (each has 16 features) ---
    size_t nSV = std::min<size_t>(features.vtx_features.size(), 5);
    std::cout << "  -- SV Candidates (data_["<<kVtxFeatures<<"], size=" 
              << (nSV * n_features_sv) << ") --" << std::endl;
    for (size_t i = 0; i < nSV; ++i) {
      const auto& sv = features.vtx_features[i];
      const float svfeats[n_features_sv] = {
        sv.jet_sv_deta,
        sv.jet_sv_dphi,
        sv.jet_sv_pt_log,
        sv.jet_sv_mass,
        static_cast<float>(sv.jet_sv_ntrack),
        sv.jet_sv_chi2,
        sv.jet_sv_dxy,
        sv.jet_sv_dxysig,
        sv.jet_sv_d3d,
        sv.jet_sv_d3dsig,
        sv.jet_sv_costhetasvpv,
        sv.jet_sv_enratio,
        sv.jet_sv_pt,
        sv.jet_sv_eta,
        sv.jet_sv_phi,
        sv.jet_sv_energy
      };
      for (int j = 0; j < n_features_sv; ++j) {
        std::cout << "    data_["<<kVtxFeatures<<"]["<<(i*n_features_sv + j)<<"]: "
                  << svfeats[j] << std::endl;
      }
    }

    std::cout << "=== End Feature Dump ===" << std::endl;
  }
#endif

  iEvent.put(std::move(output_tag_infos));
}

// Define this as a plug-in
DEFINE_FWK_MODULE(hltParticleTransformerAK4TagInfoProducer);
