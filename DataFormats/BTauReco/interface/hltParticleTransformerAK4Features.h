#ifndef DataFormats_BTauReco_hltParticleTransformerAK4Features_h
#define DataFormats_BTauReco_hltParticleTransformerAK4Features_h

#include <vector>
#include <string>

class hltGlobalFeatures {
public:
  float jet_pt;
  float jet_eta;
  float nCpfcan;
  float nNpfcan;
  float nsv;
  float npv;
  float TagVarCSV_trackSumJetEtRatio;
  float TagVarCSV_trackSumJetDeltaR;
  float TagVarCSV_vertexCategory;
  float TagVarCSV_trackSip2dValAboveCharm;
  float TagVarCSV_trackSip2dSigAboveCharm;
  float TagVarCSV_trackSip3dValAboveCharm;
  float TagVarCSV_trackSip3dSigAboveCharm;
  float TagVarCSV_jetNTracksEtaRel;
};

class hltCpfCandidateFeatures {
public:
  float Cpfcan_BtagPf_trackEtaRel;
  float Cpfcan_BtagPf_trackPtRel;
  float Cpfcan_BtagPf_trackPPar;
  float Cpfcan_BtagPf_trackDeltaR;
  float Cpfcan_BtagPf_trackPParRatio;
  float Cpfcan_BtagPf_trackSip2dVal;
  float Cpfcan_BtagPf_trackSip2dSig;
  float Cpfcan_BtagPf_trackSip3dVal;
  float Cpfcan_BtagPf_trackSip3dSig;
  float Cpfcan_BtagPf_trackJetDistVal;
  float Cpfcan_ptrel;
  float Cpfcan_drminsv;
  float Cpfcan_VTX_ass;
  float Cpfcan_puppiw;
  float Cpfcan_chi2;
  float Cpfcan_quality;
  float Cpfcan_pt;
  float Cpfcan_eta;
  float Cpfcan_phi;
  float Cpfcan_e;
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
  float jet_sv_ntrack;
  float jet_sv_mass;
  float jet_sv_energy_log;
  float jet_sv_deta;
  float jet_sv_dphi;
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
