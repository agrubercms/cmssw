#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexFeatures.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoBTag/FeatureTools/interface/SecondaryVertexConverter.h"

namespace btagbtvdeep {

  void svToFeatures(const reco::VertexCompositePtrCandidate& sv,
                    const reco::Vertex& pv,
                    const reco::Jet& jet,
                    SecondaryVertexFeatures& sv_features,
                    const bool flip) {
    math::XYZVector jet_dir = jet.momentum().Unit();
    sv_features.pt = sv.pt();
    sv_features.ptrel = sv.pt() / jet.pt();
    sv_features.etarel = catch_infs_and_bound(std::fabs(sv.eta() - jet.eta()) - 0.5, 0, -2, 0);
    sv_features.phirel = catch_infs_and_bound(std::fabs(reco::deltaPhi(sv.phi(), jet.phi())) - 0.5, 0, -2, 0);
    sv_features.eta = sv.eta();
    sv_features.phi = sv.phi();
    sv_features.e = sv.energy();
    sv_features.px = sv.px();
    sv_features.py = sv.py();
    sv_features.pz = sv.pz();
    sv_features.deltaR = catch_infs_and_bound(std::fabs(reco::deltaR(sv, jet_dir)) - 0.5, 0, -2, 0);
    sv_features.mass = sv.mass();
    sv_features.ntracks = sv.numberOfDaughters();
    sv_features.chi2 = sv.vertexChi2();
    sv_features.normchi2 = catch_infs_and_bound(sv_features.chi2 / sv.vertexNdof(), 1000, -1000, 1000);
    const auto& dxy_meas = vertexDxy(sv, pv);
    sv_features.dxy = dxy_meas.value();
    sv_features.dxysig = catch_infs_and_bound(dxy_meas.value() / dxy_meas.error(), 0, -1, 800);
    const auto& d3d_meas = vertexD3d(sv, pv);
    sv_features.d3d = d3d_meas.value();
    sv_features.d3dsig = catch_infs_and_bound(d3d_meas.value() / d3d_meas.error(), 0, -1, 800);
    sv_features.costhetasvpv = (flip ? -1.f : 1.f) * vertexDdotP(sv, pv);
    sv_features.enratio = sv.energy() / jet.energy();
  }

  // New overload for hltVtxFeatures.
  void svToFeatures(const reco::VertexCompositePtrCandidate& sv,
                    const reco::Vertex& pv,
                    const reco::Jet& jet,
                    hltVtxFeatures& sv_features,
                    const bool flip) {
    SecondaryVertexFeatures temp;
    svToFeatures(sv, pv, jet, temp, flip);
    // Map fields from temp to hltVtxFeatures (using the hlt names)
    sv_features.jet_sv_ntrack    = temp.ntracks;
    sv_features.jet_sv_mass      = temp.mass;
    sv_features.jet_sv_energy_log= (temp.e > 0 ? std::log(temp.e) : 0);
    // For deta/dphi, you may decide how to copy or recompute them.
    sv_features.jet_sv_deta      = temp.eta - jet.eta();
    sv_features.jet_sv_dphi      = std::fabs(reco::deltaPhi(temp.phi, jet.phi()));
    sv_features.jet_sv_chi2      = temp.chi2;
    sv_features.jet_sv_dxy       = temp.dxy;
    sv_features.jet_sv_dxysig    = temp.dxysig;
    sv_features.jet_sv_d3d       = temp.d3d;
    sv_features.jet_sv_d3dsig    = temp.d3dsig;
    sv_features.jet_sv_pt        = temp.pt;
    sv_features.jet_sv_eta       = temp.eta;
    sv_features.jet_sv_phi       = temp.phi;
    sv_features.jet_sv_energy    = temp.e;
  }

}  // namespace btagbtvdeep
