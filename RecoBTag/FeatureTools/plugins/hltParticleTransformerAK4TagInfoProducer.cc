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

//#define DEBUG

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
  // Updated to use reco::PFCandidate instead of pat::PackedCandidate.
  static void convertHLTChargedCandidate(const reco::PFCandidate* cand,
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
    // For PFCandidate, default the association quality.
    feat.Cpfcan_VTX_ass = 0;
    feat.Cpfcan_puppiw  = puppiw;
    
    // --- Track-related features from track_info ---
    feat.Cpfcan_BtagPf_trackEtaRel    = track_info.getTrackEtaRel();
    feat.Cpfcan_BtagPf_trackPtRel     = track_info.getTrackPtRel();
    feat.Cpfcan_BtagPf_trackPPar      = track_info.getTrackPPar();
    feat.Cpfcan_BtagPf_trackDeltaR    = track_info.getTrackDeltaR();
    feat.Cpfcan_BtagPf_trackPParRatio = track_info.getTrackPParRatio();
    feat.Cpfcan_BtagPf_trackSip2dVal  = track_info.getTrackSip2dVal();
    feat.Cpfcan_BtagPf_trackSip2dSig  = track_info.getTrackSip2dSig();
    feat.Cpfcan_BtagPf_trackSip3dVal  = track_info.getTrackSip3dVal();
    feat.Cpfcan_BtagPf_trackSip3dSig  = track_info.getTrackSip3dSig();
    feat.Cpfcan_BtagPf_trackJetDistVal = track_info.getTrackJetDistVal();
    
    if (cand->trackRef().isNonnull()) { // replaced cand->hasTrackDetails()
      const auto& track = *(cand->trackRef()); // replaced cand->pseudoTrack()
      feat.Cpfcan_chi2 = btagbtvdeep::catch_infs_and_bound(track.normalizedChi2(), 300, -1, 300);
      feat.Cpfcan_quality = track.qualityMask();
    } else {
      feat.Cpfcan_chi2 = -1;
      feat.Cpfcan_quality = (1 << reco::TrackBase::loose);
    }

#ifdef DEBUG
    std::cout << "DEBUG Charged Candidate: pt=" << cand->pt() 
              << " eta=" << cand->eta() 
              << " phi=" << cand->phi() 
              << " sip2d=" << feat.Cpfcan_BtagPf_trackSip2dVal 
              << " sip2dSig=" << feat.Cpfcan_BtagPf_trackSip2dSig 
              << " etaRel=" << feat.Cpfcan_BtagPf_trackEtaRel
              << std::endl;
#endif
  }

  static void convertHLTNeutralCandidate(const reco::PFCandidate* cand,
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
    feat.Npfcan_HadFrac = (cand->hcalEnergy())/(cand->energy()); // replaced cand->hcalFraction()
    feat.Npfcan_isGamma = cand->pdgId() == 22;

#ifdef DEBUG
    std::cout << "DEBUG Neutral Candidate: pt=" << cand->pt() 
              << " eta=" << cand->eta() 
              << " phi=" << cand->phi() 
              << " energy=" << cand->energy() 
              << " isGamma=" << feat.Npfcan_isGamma
              << std::endl;
#endif
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
            << (svs.isValid() ? ", size=" + std::to_string(svs->size()) : "") 
            << std::endl;
#endif

  // New: Retrieve handle for GenJets
  edm::Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genjet_token_, genJets);

  edm::ESHandle<TransientTrackBuilder> track_builder = iSetup.getHandle(track_builder_token_);

#ifdef DEBUG
  std::cout << "=== Debug: Processing " << jets->size() << " jets ===" << std::endl;
#endif

  // Loop over jets
  for (std::size_t jet_n = 0; jet_n < jets->size(); ++jet_n) {
    edm::RefToBase<reco::Jet> jet_ref(jets, jet_n);
    const auto& jet = jets->at(jet_n);

#ifdef DEBUG
    std::cout << "Processing jet #" << jet_n << ": pt=" << jet.pt() << " eta=" << jet.eta() << " phi=" << jet.phi() << std::endl;
#endif

    // New: Find and print matching gen jet from the "ak4GenJets" collection.
    /*const reco::GenJet* matchedGenJet = nullptr;
    double minDR2 = 1e6;
    for (const auto& gj : *genJets) {
      double dr2 = reco::deltaR2(jet, gj);
      if (dr2 < minDR2 && dr2 < 0.16) { // using 0.4 threshold squared (0.16)
        minDR2 = dr2;
        matchedGenJet = &gj;
      }
    }
    if (matchedGenJet) {
      std::cout << "  Matched ak4GenJets: pt = " << matchedGenJet->pt() 
                << ", eta = " << matchedGenJet->eta() 
                << ", phi = " << matchedGenJet->phi() 
                << ", energy = " << matchedGenJet->energy()
                << ", pdgId = " << matchedGenJet->pdgId() << std::endl;
    } else {
      std::cout << "  No matching ak4GenJets found" << std::endl;
    }*/

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

      // --- Build candidate ordering ---
      std::vector<btagbtvdeep::SortingClass<size_t>> c_sorted, n_sorted;
      std::map<unsigned int, btagbtvdeep::TrackInfoBuilder> trackinfos;
      
      // First, count valid candidates
      int validChargedCount = 0, validNeutralCount = 0;
      for (unsigned int i = 0; i < jet.numberOfDaughters(); ++i) {
        // Assume all daughters are reco::PFCandidate*
        const reco::PFCandidate* cand = static_cast<const reco::PFCandidate*>(jet.daughter(i));
        if (!cand || cand->pt() < min_candidate_pt_)
          continue;
        
        if (cand->charge() != 0) {
          validChargedCount++;
        } else {
          validNeutralCount++;
        }
      }
      
#ifdef DEBUG
      std::cout << "  Found " << validChargedCount << " charged and " << validNeutralCount 
                << " neutral candidates (above pt threshold " << min_candidate_pt_ << ")" << std::endl;
#endif

      // Pre-allocate vectors for sorting and features
      c_sorted.reserve(validChargedCount);
      n_sorted.reserve(validNeutralCount);
      
      for (unsigned int i = 0; i < jet.numberOfDaughters(); ++i) {
        // Assume all daughters are reco::PFCandidate*
        const reco::PFCandidate* cand = static_cast<const reco::PFCandidate*>(jet.daughter(i));
        if (!cand || cand->pt() < min_candidate_pt_)
          continue;

#ifdef DEBUG
        std::cout << "  Candidate #" << i << " type: " << typeid(*cand).name() << std::endl;
#endif
        
        if (cand->charge() != 0) {
          // For charged candidates, build track info
          auto& tinfo = trackinfos.emplace(i, track_builder).first->second;
          tinfo.buildTrackInfo(cand, jet.momentum().Unit(), GlobalVector(jet.px(), jet.py(), jet.pz()), pv);
          
          // Only include candidates with valid track info
          float sip2dsig = tinfo.getTrackSip2dSig();
          if (std::isfinite(sip2dsig)) {
            c_sorted.emplace_back(i,
                                 sip2dsig,
                                 -btagbtvdeep::mindrsvpfcand(*svs, cand),
                                 cand->pt() / jet.pt());
#ifdef DEBUG
            std::cout << "  Adding charged candidate idx=" << i << " pt=" << cand->pt() 
                      << " sip2dsig=" << sip2dsig << std::endl;
#endif
          }
        } else {
          n_sorted.emplace_back(i, -1, -btagbtvdeep::mindrsvpfcand(*svs, cand), cand->pt() / jet.pt());
#ifdef DEBUG
          std::cout << "  Adding neutral candidate idx=" << i << " pt=" << cand->pt() << std::endl;
#endif
        }
      }

      // Sort candidates by track properties
      std::sort(c_sorted.begin(), c_sorted.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);
      std::sort(n_sorted.begin(), n_sorted.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);
      
      auto c_sortedindices = btagbtvdeep::invertSortingVector(c_sorted);
      auto n_sortedindices = btagbtvdeep::invertSortingVector(n_sorted);

#ifdef DEBUG
      std::cout << "  After sorting: " << c_sorted.size() << " charged, " << n_sorted.size() << " neutral candidates" << std::endl;
#endif

      // Resize feature vectors to match sorted candidate counts
      hltFeatures.cpf_candidates.resize(c_sorted.size());
      hltFeatures.npf_candidates.resize(n_sorted.size());

      // --- Loop over candidates ---
      size_t chargedCount = 0, neutralCount = 0;
      for (unsigned int i = 0; i < jet.numberOfDaughters(); ++i) {
        // Assume all daughters are reco::PFCandidate*
        const reco::PFCandidate* cand = static_cast<const reco::PFCandidate*>(jet.daughter(i));
        if (!cand || cand->pt() < min_candidate_pt_) 
          continue;

        edm::Ref<edm::View<reco::Candidate>> candRef = getPersistentCandidate(cand, tracks);

        // Get PUPPI weight: if not available, default to 1.0.
        float puppiw = 1.0;
        if (fallback_puppi_weight_) {
          // Uncomment below if PFCandidate provides puppiWeight(); otherwise value remains 1.0.
          // puppiw = pf_cand->puppiWeight();
        }

        float drminpfcandsv = btagbtvdeep::mindrsvpfcand(*svs, cand);
        float distminpfcandsv = 0;

        if (cand->charge() != 0) {
          auto it = std::find_if(c_sorted.begin(), c_sorted.end(),
                                 [i](const auto& sc) { return sc.get() == i; });
          if (it == c_sorted.end()) continue;

          size_t entry = c_sortedindices.at(chargedCount);
          ++chargedCount;
          auto& c_pf_features = hltFeatures.cpf_candidates.at(entry);

          if (cand->trackRef().isNonnull()) {
            const reco::Track& pseudoTrack = *(cand->trackRef());
            reco::TransientTrack transientTrack = track_builder->build(pseudoTrack);
            distminpfcandsv = btagbtvdeep::mindistsvpfcand(*svs, transientTrack);
          }

          convertHLTChargedCandidate(cand,
                                     jet,
                                     trackinfos.at(i),
                                     is_weighted_jet_,
                                     drminpfcandsv,
                                     static_cast<float>(jet_radius_),
                                     puppiw,
                                     c_pf_features,
                                     flip_,
                                     distminpfcandsv);
        } else {
          if (neutralCount >= n_sorted.size()) continue;
          size_t entry = n_sortedindices.at(neutralCount);
          ++neutralCount;
          auto& n_pf_features = hltFeatures.npf_candidates.at(entry);

          convertHLTNeutralCandidate(cand,
                                     jet,
                                     is_weighted_jet_,
                                     drminpfcandsv,
                                     static_cast<float>(jet_radius_),
                                     puppiw,
                                     n_pf_features);
        }
      }  // end candidate loop

#ifdef DEBUG
      std::cout << "  Processed " << chargedCount << "/" << c_sorted.size() << " charged and "
                << neutralCount << "/" << n_sorted.size() << " neutral candidates" << std::endl;
#endif

      // --- Compute Global Features ---
      hltFeatures.global_features.jet_pt   = jet.pt();
      hltFeatures.global_features.jet_eta  = jet.eta();
      hltFeatures.global_features.nCpfcan  = hltFeatures.cpf_candidates.size();
      hltFeatures.global_features.nNpfcan  = hltFeatures.npf_candidates.size();
      hltFeatures.global_features.nsv      = hltFeatures.vtx_features.size();
    std::cout << std::endl;
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
    std::cout << "    nCpfcan: " << features.global_features.nCpfcan << "\n";
    std::cout << "    nNpfcan: " << features.global_features.nNpfcan << "\n";
    std::cout << "    nsv: " << features.global_features.nsv << "\n";
    std::cout << "    npv: " << features.global_features.npv << "\n";
    std::cout << "    TagVarCSV_trackSumJetEtRatio: " << features.global_features.TagVarCSV_trackSumJetEtRatio << "\n";
    std::cout << "    TagVarCSV_trackSumJetDeltaR: " << features.global_features.TagVarCSV_trackSumJetDeltaR << "\n";
    std::cout << "    TagVarCSV_vertexCategory: " << features.global_features.TagVarCSV_vertexCategory << "\n";
    std::cout << "    TagVarCSV_trackSip2dValAboveCharm: " << features.global_features.TagVarCSV_trackSip2dValAboveCharm << "\n";
    std::cout << "    TagVarCSV_trackSip2dSigAboveCharm: " << features.global_features.TagVarCSV_trackSip2dSigAboveCharm << "\n";
    std::cout << "    TagVarCSV_trackSip3dValAboveCharm: " << features.global_features.TagVarCSV_trackSip3dValAboveCharm << "\n";
    std::cout << "    TagVarCSV_trackSip3dSigAboveCharm: " << features.global_features.TagVarCSV_trackSip3dSigAboveCharm << "\n";
    std::cout << "    TagVarCSV_jetNTracksEtaRel: " << features.global_features.TagVarCSV_jetNTracksEtaRel << "\n";
    
    // Print charged PF candidates features
    std::cout << "  -- Charged PF Candidates (" << features.cpf_candidates.size() << ") --\n";
    for (size_t j = 0; j < features.cpf_candidates.size(); ++j) {
      const auto& cpf = features.cpf_candidates[j];
      std::cout << "    CPF #" << j << ":\n";
      std::cout << "      Cpfcan_BtagPf_trackEtaRel: " << cpf.Cpfcan_BtagPf_trackEtaRel << "\n";
      std::cout << "      Cpfcan_BtagPf_trackPtRel: " << cpf.Cpfcan_BtagPf_trackPtRel << "\n";
      std::cout << "      Cpfcan_BtagPf_trackPPar: " << cpf.Cpfcan_BtagPf_trackPPar << "\n";
      std::cout << "      Cpfcan_BtagPf_trackDeltaR: " << cpf.Cpfcan_BtagPf_trackDeltaR << "\n";
      std::cout << "      Cpfcan_BtagPf_trackPParRatio: " << cpf.Cpfcan_BtagPf_trackPParRatio << "\n";
      std::cout << "      Cpfcan_BtagPf_trackSip2dVal: " << cpf.Cpfcan_BtagPf_trackSip2dVal << "\n";
      std::cout << "      Cpfcan_BtagPf_trackSip2dSig: " << cpf.Cpfcan_BtagPf_trackSip2dSig << "\n";
      std::cout << "      Cpfcan_BtagPf_trackSip3dVal: " << cpf.Cpfcan_BtagPf_trackSip3dVal << "\n";
      std::cout << "      Cpfcan_BtagPf_trackSip3dSig: " << cpf.Cpfcan_BtagPf_trackSip3dSig << "\n";
      std::cout << "      Cpfcan_BtagPf_trackJetDistVal: " << cpf.Cpfcan_BtagPf_trackJetDistVal << "\n";
      std::cout << "      Cpfcan_ptrel: " << cpf.Cpfcan_ptrel << "\n";
      std::cout << "      Cpfcan_drminsv: " << cpf.Cpfcan_drminsv << "\n";
      std::cout << "      Cpfcan_VTX_ass: " << cpf.Cpfcan_VTX_ass << "\n";
      std::cout << "      Cpfcan_puppiw: " << cpf.Cpfcan_puppiw << "\n";
      std::cout << "      Cpfcan_chi2: " << cpf.Cpfcan_chi2 << "\n";
      std::cout << "      Cpfcan_quality: " << cpf.Cpfcan_quality << "\n";
      std::cout << "      Cpfcan_pt: " << cpf.Cpfcan_pt << "\n";
      std::cout << "      Cpfcan_eta: " << cpf.Cpfcan_eta << "\n";
      std::cout << "      Cpfcan_phi: " << cpf.Cpfcan_phi << "\n";
      std::cout << "      Cpfcan_e: " << cpf.Cpfcan_e << "\n";
    }
    
    // Print neutral PF candidates features
    std::cout << "  -- Neutral PF Candidates (" << features.npf_candidates.size() << ") --\n";
    for (size_t j = 0; j < features.npf_candidates.size(); ++j) {
      const auto& npf = features.npf_candidates[j];
      std::cout << "    NPF #" << j << ":\n";
      std::cout << "      Npfcan_ptrel: " << npf.Npfcan_ptrel << "\n";
      std::cout << "      Npfcan_deltaR: " << npf.Npfcan_deltaR << "\n";
      std::cout << "      Npfcan_isGamma: " << npf.Npfcan_isGamma << "\n";
      std::cout << "      Npfcan_HadFrac: " << npf.Npfcan_HadFrac << "\n";
      std::cout << "      Npfcan_drminsv: " << npf.Npfcan_drminsv << "\n";
      std::cout << "      Npfcan_puppiw: " << npf.Npfcan_puppiw << "\n";
      std::cout << "      Npfcan_pt: " << npf.Npfcan_pt << "\n";
      std::cout << "      Npfcan_eta: " << npf.Npfcan_eta << "\n";
      std::cout << "      Npfcan_phi: " << npf.Npfcan_phi << "\n";
      std::cout << "      Npfcan_energy: " << npf.Npfcan_energy << "\n";
    }
    
    // Print secondary vertex features
    std::cout << "  -- Secondary Vertices (" << features.vtx_features.size() << ") --\n";
    for (size_t j = 0; j < features.vtx_features.size(); ++j) {
      const auto& sv = features.vtx_features[j];
      std::cout << "    SV #" << j << ":\n";
      std::cout << "      jet_sv_ntrack: " << sv.jet_sv_ntrack << "\n";
      std::cout << "      jet_sv_mass: " << sv.jet_sv_mass << "\n";
      std::cout << "      jet_sv_energy_log: " << sv.jet_sv_energy_log << "\n";
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
