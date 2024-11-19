#ifndef HLTP2GTTauFilter_h
#define HLTP2GTTauFilter_h

#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/L1Trigger/interface/P2GTAlgoBlock.h"

namespace edm {
  class ConfigurationDescriptions;
}

class HLTP2GTTauFilter : public HLTFilter {
public:
  HLTP2GTTauFilter(const edm::ParameterSet&);
  ~HLTP2GTTauFilter() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  bool hltFilter(edm::Event&,
                 const edm::EventSetup&,
                 trigger::TriggerFilterObjectWithRefs& filterproduct) const override;

private:
  edm::InputTag m_l1GTAlgoBlockTag;
  edm::EDGetTokenT<l1t::P2GTAlgoBlockMap> m_algoBlockToken;
  std::vector<std::string> m_l1GTAlgoNames;
  double m_minPt;
  unsigned int m_minN;
  double m_maxAbsEta;
  bool m_saveTags;
};

#endif  //HLTP2GTTauFilter_h