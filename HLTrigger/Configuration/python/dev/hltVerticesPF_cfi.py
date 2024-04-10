import FWCore.ParameterSet.Config as cms

hltVerticesPF = cms.EDProducer( "PrimaryVertexProducer",
    vertexCollections = cms.VPSet( 
      cms.PSet(  chi2cutoff = cms.double( 3.0 ),
        label = cms.string( "" ),
        useBeamConstraint = cms.bool( False ),
        minNdof = cms.double( 0.0 ),
        maxDistanceToBeam = cms.double( 1.0 ),
        algorithm = cms.string( "AdaptiveVertexFitter" )
      ),
      cms.PSet(  chi2cutoff = cms.double( 3.0 ),
        label = cms.string( "WithBS" ),
        useBeamConstraint = cms.bool( True ),
        minNdof = cms.double( 0.0 ),
        maxDistanceToBeam = cms.double( 1.0 ),
        algorithm = cms.string( "AdaptiveVertexFitter" )
      )
    ),
    verbose = cms.untracked.bool( False ),
    TkFilterParameters = cms.PSet( 
      maxEta = cms.double( 100.0 ),
      minPt = cms.double( 0.0 ),
      minSiliconLayersWithHits = cms.int32( 5 ),
      minPixelLayersWithHits = cms.int32( 2 ),
      maxNormalizedChi2 = cms.double( 20.0 ),
      trackQuality = cms.string( "any" ),
      algorithm = cms.string( "filter" ),
      maxD0Significance = cms.double( 999.0 )
    ),
    beamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
    TrackLabel = cms.InputTag( "hltPFMuonMerging" ),
    TrackTimeResosLabel = cms.InputTag( "dummy_default" ),
    TrackTimesLabel = cms.InputTag( "dummy_default" ),
    trackMTDTimeQualityVMapTag = cms.InputTag( "dummy_default" ),
    TkClusParameters = cms.PSet( 
      TkDAClusParameters = cms.PSet( 
        zmerge = cms.double( 0.01 ),
        Tstop = cms.double( 0.5 ),
        d0CutOff = cms.double( 999.0 ),
        dzCutOff = cms.double( 4.0 ),
        vertexSize = cms.double( 0.15 ),
        coolingFactor = cms.double( 0.6 ),
        Tpurge = cms.double( 2.0 ),
        Tmin = cms.double( 2.4 ),
        uniquetrkweight = cms.double( 0.9 )
      ),
      algorithm = cms.string( "DA_vect" )
    ),
    isRecoveryIteration = cms.bool( False ),
    recoveryVtxCollection = cms.InputTag( "" ),
    useMVACut = cms.bool( False ),
    minTrackTimeQuality = cms.double( 0.8 )
)
hltVerticesPFSelector = cms.EDFilter( "PrimaryVertexObjectFilter",
    filterParams = cms.PSet( 
      maxZ = cms.double( 24.0 ),
      minNdof = cms.double( 4.0 ),
      maxRho = cms.double( 2.0 ),
      pvSrc = cms.InputTag( "hltVerticesPF" )
    ),
    src = cms.InputTag( "hltVerticesPF" )
)
hltVerticesPFFilter = cms.EDFilter( "VertexSelector",
    src = cms.InputTag( "hltVerticesPFSelector" ),
    cut = cms.string( "!isFake" ),
    filter = cms.bool( True )
)
hltPFJetForBtagSelector = cms.EDFilter( "HLT1PFJet",
    saveTags = cms.bool( True ),
    inputTag = cms.InputTag( "hltAK4PFJetsCorrected" ),
    triggerType = cms.int32( 86 ),
    MinE = cms.double( -1.0 ),
    MinPt = cms.double( 30.0 ),
    MinMass = cms.double( -1.0 ),
    MaxMass = cms.double( -1.0 ),
    MinEta = cms.double( -1.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 1 )
)
hltPFJetForBtag = cms.EDProducer( "HLTPFJetCollectionProducer",
    HLTObject = cms.InputTag( "hltPFJetForBtagSelector" ),
    TriggerTypes = cms.vint32( 86 )
)
hltDeepBLifetimeTagInfosPF = cms.EDProducer( "CandIPProducer",
    primaryVertex = cms.InputTag( "hltVerticesPFFilter" ),
    computeProbabilities = cms.bool( True ),
    computeGhostTrack = cms.bool( True ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 3 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 5.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False ),
    jets = cms.InputTag( "hltPFJetForBtag" ),
    candidates = cms.InputTag( "hltParticleFlow" ),
    maxDeltaR = cms.double( 0.4 )
)


hltPrimaryVertexAssociation = cms.EDProducer( "PFCandidatePrimaryVertexSorter",
    sorting = cms.PSet(  ),
    assignment = cms.PSet( 
      maxDxyForJetAxisAssigment = cms.double( 0.1 ),
      maxDzForJetAxisAssigment = cms.double( 0.1 ),
      useTiming = cms.bool( False ),
      preferHighRanked = cms.bool( False ),
      EtaMinUseDz = cms.double( -1.0 ),
      maxDistanceToJetAxis = cms.double( 0.07 ),
      PtMaxCharged = cms.double( -1.0 ),
      minJetPt = cms.double( 25.0 ),
      maxDzSigForPrimaryAssignment = cms.double( 5.0 ),
      OnlyUseFirstDz = cms.bool( False ),
      maxDzErrorForPrimaryAssignment = cms.double( 0.05 ),
      maxDzForPrimaryAssignment = cms.double( 0.1 ),
      maxJetDeltaR = cms.double( 0.5 ),
      maxDxySigForNotReconstructedPrimary = cms.double( 2.0 ),
      DzCutForChargedFromPUVtxs = cms.double( 0.2 ),
      maxDtSigForPrimaryAssignment = cms.double( 3.0 ),
      maxDxyForNotReconstructedPrimary = cms.double( 0.01 ),
      useVertexFit = cms.bool( True ),
      NumOfPUVtxsForCharged = cms.uint32( 0 )
    ),
    qualityForPrimary = cms.int32( 2 ),
    usePVMET = cms.bool( True ),
    particles = cms.InputTag( "hltParticleFlow" ),
    vertices = cms.InputTag( "hltVerticesPFFilter" ),
    jets = cms.InputTag( "hltPFJetForBtag" ),
    produceAssociationToOriginalVertices = cms.bool( True ),
    produceSortedVertices = cms.bool( False ),
    producePileUpCollection = cms.bool( False ),
    produceNoPileUpCollection = cms.bool( False )
)
hltParticleNetJetTagInfos = cms.EDProducer( "DeepBoostedJetTagInfoProducer",
    jet_radius = cms.double( 0.4 ),
    min_jet_pt = cms.double( 30.0 ),
    max_jet_eta = cms.double( 2.5 ),
    min_pt_for_track_properties = cms.double( 0.95 ),
    min_pt_for_pfcandidates = cms.double( 0.1 ),
    use_puppiP4 = cms.bool( False ),
    include_neutrals = cms.bool( True ),
    sort_by_sip2dsig = cms.bool( False ),
    min_puppi_wgt = cms.double( -1.0 ),
    flip_ip_sign = cms.bool( False ),
    sip3dSigMax = cms.double( -1.0 ),
    use_hlt_features = cms.bool( True ),
    vertices = cms.InputTag( "hltVerticesPFFilter" ),
    secondary_vertices = cms.InputTag( "hltDeepInclusiveMergedVerticesPF" ),
    pf_candidates = cms.InputTag( "hltParticleFlow" ),
    jets = cms.InputTag( "hltPFJetForBtag" ),
    puppi_value_map = cms.InputTag( "" ),
    vertex_associator = cms.InputTag( 'hltPrimaryVertexAssociation','original' ),
    use_scouting_features = cms.bool( False ),
    normchi2_value_map = cms.InputTag( "" ),
    dz_value_map = cms.InputTag( "" ),
    dxy_value_map = cms.InputTag( "" ),
    dzsig_value_map = cms.InputTag( "" ),
    dxysig_value_map = cms.InputTag( "" ),
    lostInnerHits_value_map = cms.InputTag( "" ),
    quality_value_map = cms.InputTag( "" ),
    trkPt_value_map = cms.InputTag( "" ),
    trkEta_value_map = cms.InputTag( "" ),
    trkPhi_value_map = cms.InputTag( "" )
)
hltParticleNetONNXJetTags = cms.EDProducer( "BoostedJetONNXJetTagsProducer",
    src = cms.InputTag( "hltParticleNetJetTagInfos" ),
    preprocess_json = cms.string( "RecoBTag/Combined/data/HLT/ParticleNetAK4/V01/preprocess.json" ),
    preprocessParams = cms.PSet(  ),
    model_path = cms.FileInPath( "RecoBTag/Combined/data/HLT/ParticleNetAK4/V01/particle-net.onnx" ),
    flav_names = cms.vstring( 'probtauhp',
      'probtauhm',
      'probb',
      'probc',
      'probuds',
      'probg',
      'ptcorr' ),
    jets = cms.InputTag( "" ),
    produceValueMap = cms.untracked.bool( False ),
    debugMode = cms.untracked.bool( False )
)
hltParticleNetDiscriminatorsJetTags = cms.EDProducer( "BTagProbabilityToDiscriminator",
    discriminators = cms.VPSet( 
      cms.PSet(  name = cms.string( "BvsAll" ),
        numerator = cms.VInputTag( 'hltParticleNetONNXJetTags:probb' ),
        denominator = cms.VInputTag( 'hltParticleNetONNXJetTags:probb','hltParticleNetONNXJetTags:probc','hltParticleNetONNXJetTags:probuds','hltParticleNetONNXJetTags:probg' )
      ),
      cms.PSet(  name = cms.string( "CvsAll" ),
        numerator = cms.VInputTag( 'hltParticleNetONNXJetTags:probc' ),
        denominator = cms.VInputTag( 'hltParticleNetONNXJetTags:probb','hltParticleNetONNXJetTags:probc','hltParticleNetONNXJetTags:probuds','hltParticleNetONNXJetTags:probg' )
      ),
      cms.PSet(  name = cms.string( "TauhvsAll" ),
        numerator = cms.VInputTag( 'hltParticleNetONNXJetTags:probtauhp','hltParticleNetONNXJetTags:probtauhm' ),
        denominator = cms.VInputTag( 'hltParticleNetONNXJetTags:probb','hltParticleNetONNXJetTags:probc','hltParticleNetONNXJetTags:probuds','hltParticleNetONNXJetTags:probg','hltParticleNetONNXJetTags:probtauhp','hltParticleNetONNXJetTags:probtauhm' )
      ),
      cms.PSet(  name = cms.string( "CvsL" ),
        numerator = cms.VInputTag( 'hltParticleNetONNXJetTags:probc' ),
        denominator = cms.VInputTag( 'hltParticleNetONNXJetTags:probc','hltParticleNetONNXJetTags:probuds','hltParticleNetONNXJetTags:probg' )
      )
    )
)
hltDoublePFJets30PNetTauhTagMediumWPL2DoubleTau = cms.EDFilter( "TauTagFilter",
    saveTags = cms.bool( True ),
    nExpected = cms.int32( 2 ),
    taus = cms.InputTag( "hltPFJetForBtag" ),
    tauTags = cms.InputTag( 'hltParticleNetDiscriminatorsJetTags','TauhvsAll' ),
    tauPtCorr = cms.InputTag( 'hltParticleNetONNXJetTags','ptcorr' ),
    seeds = cms.InputTag( "hltL2DoubleTauTagNNFilter" ),
    seedTypes = cms.vint32( -100 ),
    selection = cms.string( "double t1 = 0.56, t2 = 0.47, t3 = 0.001, t4 = 0, x1 = 30, x2 = 100, x3 = 500, x4 = 1000; if (pt <= x1) return t1; if ((pt > x1) && (pt <= x2)) return (t2 - t1) / (x2 - x1) * (pt - x1) + t1; if ((pt > x2) && (pt <= x3)) return (t3 - t2) / (x3 - x2) * (pt - x2) + t2; if ((pt > x3) && (pt <= x4)) return (t4 - t3) / (x4 - x3) * (pt - x3) + t3; return t4;" ),
    minPt = cms.double( 30.0 ),
    maxEta = cms.double( 2.3 ),
    usePtCorr = cms.bool( True ),
    matchWithSeeds = cms.bool( True ),
    matchingdR = cms.double( 0.5 )
)