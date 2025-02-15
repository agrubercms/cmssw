#ifndef RecoBTag_FeatureTools_SecondaryVertexConverter_h
#define RecoBTag_FeatureTools_SecondaryVertexConverter_h

#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexFeatures.h"
#include "DataFormats/BTauReco/interface/hltParticleTransformerAK4Features.h" // for hltVtxFeatures
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

namespace btagbtvdeep {

  // Overload for standard SecondaryVertexFeatures.
  void svToFeatures(const reco::VertexCompositePtrCandidate& sv,
                    const reco::Vertex& pv,
                    const reco::Jet& jet,
                    SecondaryVertexFeatures& sv_features,
                    const bool flip = false);

  // New overload for hltVtxFeatures.
  void svToFeatures(const reco::VertexCompositePtrCandidate& sv,
                    const reco::Vertex& pv,
                    const reco::Jet& jet,
                    hltVtxFeatures& sv_features,
                    const bool flip);
                    
}

#endif  //RecoBTag_FeatureTools_SecondaryVertexConverter_h
