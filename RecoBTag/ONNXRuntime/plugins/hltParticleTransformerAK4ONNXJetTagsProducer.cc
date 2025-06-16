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

    // Global features: fixed size now 4.
    size_t global_size = 4;

    // For charged PF candidates:
    constexpr static unsigned n_max_cpf_candidates_ = 50; // updated maximum candidates per jet
    unsigned n_features_cpf_; // now 22 features per candidate

    // For neutral PF candidates:
    constexpr static unsigned n_max_npf_candidates_ = 0; // no neutral candidates
    unsigned n_features_npf_; // set to 0

    // For SV candidates:
    constexpr static unsigned n_max_sv_candidates_ = 5;  // remains
    unsigned n_features_sv_; // updated to 14

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
            {(int64_t)1, (int64_t)global_size},              // global features: shape [1,4]
            {(int64_t)1, (int64_t)n_max_cpf_candidates_, (int64_t)n_features_cpf_},         // cpf features: shape [1,50,22]
            {(int64_t)1, (int64_t)n_max_npf_candidates_, (int64_t)n_features_npf_},         // npf features: shape [1,0,0]
            {(int64_t)1, (int64_t)n_max_sv_candidates_,  (int64_t)n_features_sv_}           // vtx features: shape [1,5,14]
        };
        
        // Run inference - directly assign outputs using the same pattern as UnifiedParticleTransformer
        outputs = globalCache()->run(input_names_, data_, input_shapes_, output_names_, 1)[0];
        
        assert(outputs.size() == flav_names_.size());
        // Debug output with clearer labels
        std::cout << "Jet " << jet_n << " probabilities:" << std::endl;
        for (size_t i = 0; i < outputs.size(); ++i) {
          std::cout << "  " << flav_names_[i] << ": " << outputs[i] << std::endl;
        }
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
    n_features_cpf_ = 22; // updated: 22 features per charged candidate (removed frompv, track_chi2, track_qual, dzsig)
    n_features_npf_ = 0;  // updated: 0 features per neutral candidate
    n_features_sv_  = 14; // updated: 14 features per vertex candidate

    std::vector<unsigned int> input_sizes = {
        static_cast<unsigned int>(global_size),
        static_cast<unsigned int>(n_max_cpf_candidates_ * n_features_cpf_),
        static_cast<unsigned int>(n_max_npf_candidates_ * n_features_npf_),
        static_cast<unsigned int>(n_max_sv_candidates_ * n_features_sv_)
    };

    data_.clear();
    for (const auto& len : input_sizes) {
      data_.emplace_back(1 * len, 0);
    }

    make_inputs(features);
  }

  void hltParticleTransformerAK4ONNXJetTagsProducer::make_inputs(const btagbtvdeep::hltParticleTransformerAK4Features& features) {
    float* ptr = nullptr;
    // Global features: new order: jet_pt, jet_eta, jet_phi, jet_energy
    {
      assert(data_[kGlobalFeatures].size() >= global_size);
      float* start = &data_[kGlobalFeatures][0];
      ptr = start;
      *ptr++ = features.global_features.jet_pt;
      *ptr++ = features.global_features.jet_eta;
      *ptr++ = features.global_features.jet_phi;
      *ptr++ = features.global_features.jet_energy;
      assert(ptr == start + global_size);
    }
    
    // Charged PF candidates (new order, 22 features per candidate):
    {
      assert(data_[kCpfCandidates].size() >= n_max_cpf_candidates_ * n_features_cpf_);
      unsigned offset = 0;
      for (std::size_t c_pf_n = 0; c_pf_n < std::min(features.cpf_candidates.size(), (std::size_t)n_max_cpf_candidates_); c_pf_n++) {
        ptr = &data_[kCpfCandidates][offset + c_pf_n * n_features_cpf_];
        const auto& cpf = features.cpf_candidates[c_pf_n];
        float* start_cpf = ptr;
        *ptr++ = cpf.jet_pfcand_deta;
        *ptr++ = cpf.jet_pfcand_dphi;
        *ptr++ = cpf.jet_pfcand_pt_log;
        *ptr++ = cpf.jet_pfcand_energy_log;
        *ptr++ = cpf.jet_pfcand_charge;
        *ptr++ = cpf.jet_pfcand_nlostinnerhits;
        *ptr++ = cpf.jet_pfcand_dz;
        *ptr++ = cpf.jet_pfcand_dxy;
        *ptr++ = cpf.jet_pfcand_dxysig;
        *ptr++ = cpf.jet_pfcand_etarel;
        *ptr++ = cpf.jet_pfcand_pperp_ratio;
        *ptr++ = cpf.jet_pfcand_ppara_ratio;
        *ptr++ = cpf.jet_pfcand_trackjet_d3d;
        *ptr++ = cpf.jet_pfcand_trackjet_d3dsig;
        *ptr++ = cpf.jet_pfcand_trackjet_dist;
        *ptr++ = cpf.jet_pfcand_trackjet_decayL;
        *ptr++ = cpf.jet_pfcand_npixhits;
        *ptr++ = cpf.jet_pfcand_nstriphits;
        *ptr++ = cpf.jet_pfcand_pt;
        *ptr++ = cpf.jet_pfcand_eta;
        *ptr++ = cpf.jet_pfcand_phi;
        *ptr++ = cpf.jet_pfcand_energy;
        int written = ptr - start_cpf;
        if (written != static_cast<int>(n_features_cpf_)) {
          std::cout << "Charged candidate " << c_pf_n << ": wrote " << written 
                    << " features, expected " << n_features_cpf_ << std::endl;
        }
        assert(written == static_cast<int>(n_features_cpf_));
      }
    }
    
    // Neutral PF candidates: none expected.
    {
      assert(data_[kNpfCandidates].size() == 0 || data_[kNpfCandidates].empty());
    }
    
    // SV candidates (new order: jet_sv_deta, jet_sv_dphi, jet_sv_pt_log, jet_sv_mass,
    // jet_sv_ntrack, jet_sv_chi2, jet_sv_dxy, jet_sv_dxysig, jet_sv_d3d, jet_sv_d3dsig,
    // jet_sv_pt, jet_sv_eta, jet_sv_phi, jet_sv_energy):
    {
      assert(data_[kVtxFeatures].size() >= n_max_sv_candidates_ * n_features_sv_);
      unsigned offset = 0;
      for (std::size_t sv_n = 0; sv_n < std::min(features.vtx_features.size(), (std::size_t)n_max_sv_candidates_); sv_n++) {
        ptr = &data_[kVtxFeatures][offset + sv_n * n_features_sv_];
        const auto& sv = features.vtx_features[sv_n];
        float* start_sv = ptr;
        *ptr++ = sv.jet_sv_deta;
        *ptr++ = sv.jet_sv_dphi;
        *ptr++ = sv.jet_sv_pt_log;
        *ptr++ = sv.jet_sv_mass;
        *ptr++ = sv.jet_sv_ntrack;
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
        // --- Printout of the filled tensors ---
    std::cout << "=== Dumping Tensors for ONNX Runtime ===" << std::endl;
    // Print Global Features
    std::cout << "  -- Global Features (data_[" << kGlobalFeatures << "], size: " << data_[kGlobalFeatures].size() << ") --" << std::endl;
    for (size_t i = 0; i < data_[kGlobalFeatures].size(); ++i) {
      std::cout << "    data_[" << kGlobalFeatures << "][" << i << "]: " << data_[kGlobalFeatures][i] << std::endl;
    }

    // Print Charged PF Candidates Features
    std::cout << "  -- Charged PF Candidates (data_[" << kCpfCandidates << "], size: " << data_[kCpfCandidates].size() << ") --" << std::endl;
    std::cout << "    (Formatted as [cand_idx * n_features + feature_idx])" << std::endl;
    for (size_t i = 0; i < data_[kCpfCandidates].size(); ++i) {
      // Print only a few elements for brevity if the tensor is too large, or print all
      // For now, printing all. Can be adjusted if output is too verbose.
      std::cout << "    data_[" << kCpfCandidates << "][" << i << "]: " << data_[kCpfCandidates][i] << std::endl;
    }

    // Print Neutral PF Candidates Features
    std::cout << "  -- Neutral PF Candidates (data_[" << kNpfCandidates << "], size: " << data_[kNpfCandidates].size() << ") --" << std::endl;
    for (size_t i = 0; i < data_[kNpfCandidates].size(); ++i) {
      std::cout << "    data_[" << kNpfCandidates << "][" << i << "]: " << data_[kNpfCandidates][i] << std::endl;
    }

    // Print Vertex Features
    std::cout << "  -- Vertex Features (data_[" << kVtxFeatures << "], size: " << data_[kVtxFeatures].size() << ") --" << std::endl;
    std::cout << "    (Formatted as [vtx_idx * n_features + feature_idx])" << std::endl;
    for (size_t i = 0; i < data_[kVtxFeatures].size(); ++i) {
      std::cout << "    data_[" << kVtxFeatures << "][" << i << "]: " << data_[kVtxFeatures][i] << std::endl;
    }
    std::cout << "=== End Tensor Dump ===" << std::endl;
    // --- End Printout ---
  }

//define this as a plug-in
DEFINE_FWK_MODULE(hltParticleTransformerAK4ONNXJetTagsProducer);
