#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/hltParticleTransformerAK4Features.h"
#include "DataFormats/BTauReco/interface/hltParticleTransformerAK4TagInfo.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoBTag/ONNXRuntime/interface/tensor_fillers.h"
#include "RecoBTag/ONNXRuntime/interface/tensor_configs.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace cms::Ort;

class hltParticleTransformerAK4ONNXJetTagsProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
  public:
    explicit hltParticleTransformerAK4ONNXJetTagsProducer(const edm::ParameterSet&, const ONNXRuntime*);
    ~hltParticleTransformerAK4ONNXJetTagsProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions&);

    static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
    static void globalEndJob(const ONNXRuntime*);

  private:
    typedef std::vector<reco::hltParticleTransformerAK4TagInfo> TagInfoCollection;
    typedef reco::JetTagCollection JetTagCollection;

    void produce(edm::Event&, const edm::EventSetup&) override;

    void make_inputs(const btagbtvdeep::hltParticleTransformerAK4Features& features);
    void get_input_sizes(const reco::hltParticleTransformerAK4TagInfo& taginfo);

    const edm::EDGetTokenT<TagInfoCollection> src_;
    std::vector<std::string> flav_names_;
    std::vector<std::string> input_names_;
    std::vector<std::string> output_names_;

    enum InputIndexes {
      kGlobalFeatures = 0,
      kCpfCandidates = 1,
      kNpfCandidates = 2,
      kVtxFeatures = 3
    };

    // Global features: fixed size remains.
    size_t global_size = 14;

    // For charged PF candidates:
    constexpr static unsigned n_max_cpf_candidates_ = 26; // maximum candidates per jet
    unsigned n_features_cpf_; // number of features per candidate (should be 20)

    // For neutral PF candidates:
    constexpr static unsigned n_max_npf_candidates_ = 25; // maximum candidates per jet
    unsigned n_features_npf_; // per candidate (should be 10)

    // For SV candidates:
    constexpr static unsigned n_max_sv_candidates_ = 5;  // maximum candidates per jet
    unsigned n_features_sv_; // per candidate (should be 14)

    std::vector<unsigned> input_sizes_;
    std::vector<std::vector<int64_t>> input_shapes_;  // shapes of each input group (-1 for dynamic axis)

    // hold the input data
    FloatArrays data_;
  };

  hltParticleTransformerAK4ONNXJetTagsProducer::hltParticleTransformerAK4ONNXJetTagsProducer(const edm::ParameterSet& iConfig,
                                                                                             const ONNXRuntime* cache)
      : src_(consumes<TagInfoCollection>(iConfig.getParameter<edm::InputTag>("src"))),
        flav_names_(iConfig.getParameter<std::vector<std::string>>("flav_names")),
        input_names_(iConfig.getParameter<std::vector<std::string>>("input_names")),
        output_names_(iConfig.getParameter<std::vector<std::string>>("output_names")) {
    // get output names from flav_names
    for (const auto& flav_name : flav_names_) {
      produces<JetTagCollection>(flav_name);
    }
  }

  void hltParticleTransformerAK4ONNXJetTagsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // hltParticleTransformerAK4JetTags
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("hltParticleTransformerAK4TagInfos"));
    desc.add<std::vector<std::string>>("input_names", {"global_features", "cpf_features", "npf_features", "vtx_features"});
    desc.add<edm::FileInPath>("model_path", edm::FileInPath("ParticleTransformer.onnx"));
    desc.add<std::vector<std::string>>("output_names", {"output"});
    desc.add<std::vector<std::string>>("flav_names", {"probb", "probbb", "problepb"});

    descriptions.add("hltParticleTransformerAK4ONNXJetTags", desc);
  }

  std::unique_ptr<ONNXRuntime> hltParticleTransformerAK4ONNXJetTagsProducer::initializeGlobalCache(
      const edm::ParameterSet& iConfig) {
    return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
  }

  void hltParticleTransformerAK4ONNXJetTagsProducer::globalEndJob(const ONNXRuntime* cache) {}

  void hltParticleTransformerAK4ONNXJetTagsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<TagInfoCollection> tag_infos;
    iEvent.getByToken(src_, tag_infos);

    // initialize output collection
    std::vector<std::unique_ptr<JetTagCollection>> output_tags;
    if (!tag_infos->empty()) {
      auto jet_ref = tag_infos->begin()->jet();
      auto ref2prod = edm::makeRefToBaseProdFrom(jet_ref, iEvent);
      for (std::size_t i = 0; i < flav_names_.size(); i++) {
        output_tags.emplace_back(std::make_unique<JetTagCollection>(ref2prod));
      }
    } else {
      for (std::size_t i = 0; i < flav_names_.size(); i++) {
        output_tags.emplace_back(std::make_unique<JetTagCollection>());
      }
    }

    for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
      const auto& taginfo = (*tag_infos)[jet_n];
      std::vector<float> outputs(flav_names_.size(), -1.0);
      if (taginfo.features().is_filled) {
        get_input_sizes(taginfo);

        // Restore proper feature geometry per candidate group:
        input_shapes_ = {
            {(int64_t)1, (int64_t)1, (int64_t)global_size},              // global features
            {(int64_t)1, (int64_t)n_max_cpf_candidates_, (int64_t)n_features_cpf_},         // charged PF candidates
            {(int64_t)1, (int64_t)n_max_npf_candidates_, (int64_t)n_features_npf_},         // neutral PF candidates
            {(int64_t)1, (int64_t)n_max_sv_candidates_,  (int64_t)n_features_sv_}           // SV candidates
        };
        assert(outputs.size() == flav_names_.size());
      }

      const auto& jet_ref = taginfo.jet();
      for (std::size_t flav_n = 0; flav_n < flav_names_.size(); flav_n++) {
        (*(output_tags[flav_n]))[jet_ref] = outputs[flav_n];
      }
    }

    // put into the event
    for (std::size_t flav_n = 0; flav_n < flav_names_.size(); ++flav_n) {
      iEvent.put(std::move(output_tags[flav_n]), flav_names_[flav_n]);
    }
  }

  void hltParticleTransformerAK4ONNXJetTagsProducer::get_input_sizes(
      const reco::hltParticleTransformerAK4TagInfo& taginfo) {
    const auto& features = taginfo.features();
    // Set per-candidate feature dimensions
    n_features_cpf_ = 20; // charged PF: 20 features per candidate
    n_features_npf_ = 10; // neutral PF: 10 features per candidate
    n_features_sv_  = 14; // SV: 14 features per candidate

    // Use maximum candidate counts
    const auto n_cpf_candidates = std::min(features.cpf_candidates.size(), (size_t)n_max_cpf_candidates_);
    const auto n_npf_candidates = std::min(features.npf_candidates.size(), (size_t)n_max_npf_candidates_);
    const auto n_sv_candidates  = std::min(features.vtx_features.size(),  (size_t)n_max_sv_candidates_);

    std::vector<unsigned int> input_sizes = std::vector<unsigned int>{ 
        static_cast<unsigned int>(global_size),
        static_cast<unsigned int>(n_max_cpf_candidates_ * n_features_cpf_),
        static_cast<unsigned int>(n_max_npf_candidates_ * n_features_npf_),
        static_cast<unsigned int>(n_max_sv_candidates_ * n_features_sv_)
    };

    // init data storage
    data_.clear();
    for (const auto& len : input_sizes) {
      data_.emplace_back(1 * len, 0);
    }

    make_inputs(features);
  }

  void hltParticleTransformerAK4ONNXJetTagsProducer::make_inputs(const btagbtvdeep::hltParticleTransformerAK4Features& features) {
    float* ptr = nullptr;
    unsigned offset = 0;

    // Global features: fill in order of declaration
    assert(data_[kGlobalFeatures].size() >= global_size);
    float* start = &data_[kGlobalFeatures][0]; // store start pointer
    ptr = start;
    *ptr++ = features.global_features.jet_pt;
    *ptr++ = features.global_features.jet_eta;
    *ptr++ = features.global_features.nCpfcan;
    *ptr++ = features.global_features.nNpfcan;
    *ptr++ = features.global_features.nsv;
    *ptr++ = features.global_features.npv;
    *ptr++ = features.global_features.TagVarCSV_trackSumJetEtRatio;
    *ptr++ = features.global_features.TagVarCSV_trackSumJetDeltaR;
    *ptr++ = features.global_features.TagVarCSV_vertexCategory;
    *ptr++ = features.global_features.TagVarCSV_trackSip2dValAboveCharm;
    *ptr++ = features.global_features.TagVarCSV_trackSip2dSigAboveCharm;
    *ptr++ = features.global_features.TagVarCSV_trackSip3dValAboveCharm;
    *ptr++ = features.global_features.TagVarCSV_trackSip3dSigAboveCharm;
    *ptr++ = features.global_features.TagVarCSV_jetNTracksEtaRel;
    // Assert that exactly 14 values were written
    assert(ptr == start + global_size);
    
    // Charged PF candidates
    assert(data_[kCpfCandidates].size() >= n_max_cpf_candidates_ * n_features_cpf_);
    offset = 0;
    for (std::size_t c_pf_n = 0; c_pf_n < std::min(features.cpf_candidates.size(), (std::size_t)n_max_cpf_candidates_); c_pf_n++) {
      ptr = &data_[kCpfCandidates][offset + c_pf_n * n_features_cpf_];
      const auto& cpf = features.cpf_candidates[c_pf_n];
      float* start_cpf = ptr; // store pointer start for this candidate
      *ptr++ = cpf.Cpfcan_BtagPf_trackEtaRel;
      *ptr++ = cpf.Cpfcan_BtagPf_trackPtRel;
      *ptr++ = cpf.Cpfcan_BtagPf_trackPPar;
      *ptr++ = cpf.Cpfcan_BtagPf_trackDeltaR;
      *ptr++ = cpf.Cpfcan_BtagPf_trackPParRatio;
      *ptr++ = cpf.Cpfcan_BtagPf_trackSip2dVal;
      *ptr++ = cpf.Cpfcan_BtagPf_trackSip2dSig;
      *ptr++ = cpf.Cpfcan_BtagPf_trackSip3dVal;
      *ptr++ = cpf.Cpfcan_BtagPf_trackSip3dSig;
      *ptr++ = cpf.Cpfcan_BtagPf_trackJetDistVal;
      *ptr++ = cpf.Cpfcan_ptrel;
      *ptr++ = cpf.Cpfcan_drminsv;
      *ptr++ = cpf.Cpfcan_VTX_ass;
      *ptr++ = cpf.Cpfcan_puppiw;
      *ptr++ = cpf.Cpfcan_chi2;
      *ptr++ = cpf.Cpfcan_quality;
      *ptr++ = cpf.Cpfcan_pt;
      *ptr++ = cpf.Cpfcan_eta;
      *ptr++ = cpf.Cpfcan_phi;
      *ptr++ = cpf.Cpfcan_e;
      // Now expect 20 features.
      int written = ptr - start_cpf;
      if (written != static_cast<int>(n_features_cpf_)) {
          std::cout << "Charged candidate " << c_pf_n << ": wrote " << written 
                    << " features, expected " << n_features_cpf_ << std::endl;
      }
      assert(written == static_cast<int>(n_features_cpf_));
    }
    
    // Neutral PF candidates
    assert(data_[kNpfCandidates].size() >= n_max_npf_candidates_ * n_features_npf_);
    offset = 0;
    for (std::size_t n_pf_n = 0; n_pf_n < std::min(features.npf_candidates.size(), (std::size_t)n_max_npf_candidates_); n_pf_n++) {
      ptr = &data_[kNpfCandidates][offset + n_pf_n * n_features_npf_];
      const auto& npf = features.npf_candidates[n_pf_n];
      float* start_npf = ptr; // store pointer start for this candidate
      *ptr++ = npf.Npfcan_ptrel;
      *ptr++ = npf.Npfcan_deltaR;
      *ptr++ = npf.Npfcan_isGamma;
      *ptr++ = npf.Npfcan_HadFrac;
      *ptr++ = npf.Npfcan_drminsv;
      *ptr++ = npf.Npfcan_puppiw;
      *ptr++ = npf.Npfcan_pt;
      *ptr++ = npf.Npfcan_eta;
      *ptr++ = npf.Npfcan_phi;
      *ptr++ = npf.Npfcan_energy;
      int writtenNpf = ptr - start_npf;
      if (writtenNpf != static_cast<int>(n_features_npf_)) {
          std::cout << "Neutral candidate " << n_pf_n << ": wrote " << writtenNpf
                    << " features, expected " << n_features_npf_ << std::endl;
      }
      assert(writtenNpf == static_cast<int>(n_features_npf_));
    }
    
    // SV candidates
    assert(data_[kVtxFeatures].size() >= n_max_sv_candidates_ * n_features_sv_);
    offset = 0;
    for (std::size_t sv_n = 0; sv_n < std::min(features.vtx_features.size(), (std::size_t)n_max_sv_candidates_); sv_n++) {
      ptr = &data_[kVtxFeatures][offset + sv_n * n_features_sv_];
      const auto& sv = features.vtx_features[sv_n];
      float* start_sv = ptr; // store pointer start for this candidate
      *ptr++ = sv.jet_sv_ntrack;
      *ptr++ = sv.jet_sv_mass;
      *ptr++ = sv.jet_sv_energy_log;
      *ptr++ = sv.jet_sv_deta;
      *ptr++ = sv.jet_sv_dphi;
      *ptr++ = sv.jet_sv_chi2;
      *ptr++ = sv.jet_sv_dxy;
      *ptr++ = sv.jet_sv_dxysig;
      *ptr++ = sv.jet_sv_d3d;
      *ptr++ = sv.jet_sv_d3dsig;
      *ptr++ = sv.jet_sv_pt;
      *ptr++ = sv.jet_sv_eta;
      *ptr++ = sv.jet_sv_phi;
      *ptr++ = sv.jet_sv_energy;
      int writtenSv = ptr - start_sv;
      if (writtenSv != static_cast<int>(n_features_sv_)) {
          std::cout << "SV candidate " << sv_n << ": wrote " << writtenSv
                    << " features, expected " << n_features_sv_ << std::endl;
      }
      assert(writtenSv == static_cast<int>(n_features_sv_));
    }

  }

//define this as a plug-in
DEFINE_FWK_MODULE(hltParticleTransformerAK4ONNXJetTagsProducer);
