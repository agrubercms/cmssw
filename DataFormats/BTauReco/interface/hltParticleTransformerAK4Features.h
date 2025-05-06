#ifndef DataFormats_BTauReco_hltParticleTransformerAK4Features_h
#define DataFormats_BTauReco_hltParticleTransformerAK4Features_h

#include <vector>
#include <string>

class hltGlobalFeatures {
public:
  float jet_pt;
  float jet_eta;
  float jet_phi;
  float jet_energy;
};

class hltCpfCandidateFeatures {
public:
  float jet_pfcand_deta;
  float jet_pfcand_dphi;
  float jet_pfcand_pt_log;
  float jet_pfcand_energy_log;
  float jet_pfcand_charge;
  float jet_pfcand_frompv;
  float jet_pfcand_nlostinnerhits;
  float jet_pfcand_track_chi2;
  float jet_pfcand_track_qual;
  float jet_pfcand_dz;
  float jet_pfcand_dzsig;
  float jet_pfcand_dxy;
  float jet_pfcand_dxysig;
  float jet_pfcand_etarel;
  float jet_pfcand_pperp_ratio;
  float jet_pfcand_ppara_ratio;
  float jet_pfcand_trackjet_d3d;
  float jet_pfcand_trackjet_d3dsig;
  float jet_pfcand_trackjet_dist;
  float jet_pfcand_trackjet_decayL;
  float jet_pfcand_npixhits;
  float jet_pfcand_nstriphits;
  float jet_pfcand_pt;
  float jet_pfcand_eta;
  float jet_pfcand_phi;
  float jet_pfcand_energy;
};

class hltNpfCandidateFeatures {
public:
  float Npfcan_ptrel;
  float Npfcan_deltaR;
  float Npfcan_isGamma;
  float Npfcan_HadFrac;
  float Npfcan_drminsv;
  float Npfcan_puppiw;
  float Npfcan_pt;
  float Npfcan_eta;
  float Npfcan_phi;
  float Npfcan_energy;
};

class hltVtxFeatures {
public:
  float jet_sv_deta;
  float jet_sv_dphi;
  float jet_sv_pt_log;  // updated name from jet_sv_energy_log
  float jet_sv_mass;
  float jet_sv_ntrack;
  float jet_sv_chi2;
  float jet_sv_dxy;
  float jet_sv_dxysig;
  float jet_sv_d3d;
  float jet_sv_d3dsig;
  float jet_sv_pt;
  float jet_sv_eta;
  float jet_sv_phi;
  float jet_sv_energy;
};

namespace btagbtvdeep {

  class hltParticleTransformerAK4Features {
  public:
    bool is_filled = true;
    hltGlobalFeatures global_features;
    std::vector<hltCpfCandidateFeatures> cpf_candidates;
    std::vector<hltNpfCandidateFeatures> npf_candidates;
    std::vector<hltVtxFeatures> vtx_features;
  };

}  // namespace btagbtvdeep

#endif
