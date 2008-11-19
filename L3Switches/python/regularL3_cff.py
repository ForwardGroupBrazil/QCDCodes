# /dev/CMSSW_2_1_5/HLT/V11 (CMSSW_2_1_5_HLT1)

import FWCore.ParameterSet.Config as cms


HLTConfigVersion = cms.PSet(
  tableName = cms.string('/dev/CMSSW_2_1_5/HLT/V11')
)

KFFittingSmoother = cms.ESProducer(
    "KFFittingSmootherESProducer",
    ComponentName = cms.string( "KFFittingSmoother" ),
    Fitter = cms.string( "KFFitter" ),
    Smoother = cms.string( "KFSmoother" ),
    EstimateCut = cms.double( -1.0 ),
    MinNumberOfHits = cms.int32( 5 ),
    RejectTracks = cms.bool( True ),
    BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
    NoInvalidHitsBeginEnd = cms.bool( False ),
    appendToDataLabel = cms.string( "" )
    )
KFFitter = cms.ESProducer(
    "KFTrajectoryFitterESProducer",
    ComponentName = cms.string( "KFFitter" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    Updator = cms.string( "KFUpdator" ),
    Estimator = cms.string( "Chi2" ),
    minHits = cms.int32( 3 ),
    appendToDataLabel = cms.string( "" )
    )
KFSmoother = cms.ESProducer(
    "KFTrajectorySmootherESProducer",
    ComponentName = cms.string( "KFSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    Updator = cms.string( "KFUpdator" ),
    Estimator = cms.string( "Chi2" ),
    errorRescaling = cms.double( 100.0 ),
    minHits = cms.int32( 3 ),
    appendToDataLabel = cms.string( "" )
    )

hltL3TkTracksFromL2 = cms.EDProducer(
     "TrackProducer",
     TrajectoryInEvent = cms.bool( True ),
     useHitsSplitting = cms.bool( False ),
     clusterRemovalInfo = cms.InputTag( "" ),
     alias = cms.untracked.string( "" ),
     Fitter = cms.string( "KFFittingSmoother" ),
     Propagator = cms.string( "PropagatorWithMaterial" ),
     src = cms.InputTag( "hltL3TrackCandidateFromL2" ),
     beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
     TTRHBuilder = cms.string( "WithTrackAngle" ),
     AlgorithmName = cms.string( "undefAlgorithm" )
  )
 

hltL3Muons = cms.EDProducer( "L3MuonProducer",
    MuonCollectionLabel = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    L3TrajBuilderParameters = cms.PSet( 
      TrackTransformer = cms.PSet( 
        Fitter = cms.string( "KFFitterForRefitInsideOut" ),
        Smoother = cms.string( "KFSmootherForRefitInsideOut" ),
        TrackerRecHitBuilder = cms.string( "WithTrackAngle" ),
        MuonRecHitBuilder = cms.string( "MuonRecHitBuilder" ),
        RefitDirection = cms.string( "insideOut" ),
        RefitRPCHits = cms.bool( True ),
        DoPredictionsOnly = cms.bool( False )
      ),
      GlobalMuonTrackMatcher = cms.PSet( 
        Chi2Cut = cms.double( 50.0 ),
        MinP = cms.double( 2.5 ),
        MinPt = cms.double( 1.0 ),
        DeltaRCut = cms.double( 0.2 ),
        DeltaDCut = cms.double( 10.0 )
      ),
      TrackerPropagator = cms.string( "SteppingHelixPropagatorAny" ),
      Chi2CutRPC = cms.double( 1.0 ),
      Chi2CutCSC = cms.double( 150.0 ),
      Chi2CutDT = cms.double( 10.0 ),
      HitThreshold = cms.int32( 1 ),
      Chi2ProbabilityCut = cms.double( 30.0 ),
      PtCut = cms.double( 1.0 ),
      Direction = cms.int32( 0 ),
      MuonHitsOption = cms.int32( 1 ),
      ScaleTECxFactor = cms.double(-1.0),
      ScaleTECyFactor = cms.double(-1.0),
      TrackRecHitBuilder = cms.string( "WithTrackAngle" ),
      RPCRecSegmentLabel = cms.InputTag( "hltRpcRecHits" ),
      CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
      DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
      MuonTrackingRegionBuilder = cms.PSet( 
        beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
        UseVertex = cms.bool( False ),
        Rescale_eta = cms.double( 3.0 ),
        Rescale_phi = cms.double( 3.0 ),
        Rescale_Dz = cms.double( 3.0 ),
        DeltaZ_Region = cms.double( 15.9 ),
        DeltaR = cms.double( 0.2 ),
        EscapePt = cms.double( 1.5 ),
        Phi_min = cms.double( 0.05 ),
        Eta_min = cms.double( 0.05 ),
        UseFixedRegion = cms.bool( False ),
        EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
        EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
        PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
        PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
        vertexCollection = cms.InputTag( "pixelVertices" ),
        Eta_fixed = cms.double( 0.2 ),
        Phi_fixed = cms.double( 0.2 ),
        OnDemand = cms.double( -1.0 )
      ),
      StateOnTrackerBoundOutPropagator = cms.string( "SmartPropagatorAny" ),
#xx2      l3SeedLabel = cms.InputTag( "" ),
      tkTrajLabel = cms.InputTag( "hltL3TkTracksFromL2" ),
      TkTrackBuilder = cms.string( "muonCkfTrajectoryBuilder" ),
#xx1      SeedGeneratorParameters = cms.PSet( 
#xx1        ComponentName = cms.string( "TSGFromOrderedHits" ),
#xx1        TTRHBuilder = cms.string( "WithTrackAngle" ),
#xx1        OrderedHitsFactoryPSet = cms.PSet( 
#xx1          ComponentName = cms.string( "StandardHitPairGenerator" ),
#xx1          SeedingLayers = cms.string( "PixelLayerPairs" )
#xx1        )
#xx1      ),
      KFFitter = cms.string( "L3MuKFFitter" ),
      TransformerOutPropagator = cms.string( "SmartPropagatorAny" ),
      MatcherOutPropagator = cms.string( "SmartPropagator" ),
      TrackerRecHitBuilder = cms.string( "WithTrackAngle" ),
      MuonRecHitBuilder = cms.string( "MuonRecHitBuilder" ),
      RefitRPCHits = cms.bool( True )
    ),
    ServiceParameters = cms.PSet( 
      UseMuonNavigation = cms.untracked.bool( True ),
      RPCLayers = cms.bool( True ),
      Propagators = cms.untracked.vstring( 'SteppingHelixPropagatorAny',
        'SteppingHelixPropagatorAlong',
        'SteppingHelixPropagatorOpposite',
        'PropagatorWithMaterial',
        'PropagatorWithMaterialOpposite',
        'SmartPropagator',
        'SmartPropagatorOpposite',
        'SmartPropagatorAnyOpposite',
        'SmartPropagatorAny',
        'SmartPropagatorRK',
        'SmartPropagatorAnyRK' )
    ),
    TrackLoaderParameters = cms.PSet( 
      MuonUpdatorAtVertexParameters = cms.PSet( 
        Propagator = cms.string( "SteppingHelixPropagatorOpposite" ),
        BeamSpotPosition = cms.vdouble( 0.0, 0.0, 0.0 ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 ),
        MaxChi2 = cms.double( 1000000.0 )
      ),
      VertexConstraint = cms.bool( False ),
      PutTkTrackIntoEvent = cms.untracked.bool( True ),
      MuonSeededTracksInstance = cms.untracked.string( "L2Seeded" ),
      DoSmoothing = cms.bool( True ),
      SmoothTkTrack = cms.untracked.bool( False ),
      Smoother = cms.string( "KFSmootherForMuonTrackLoader" )
    )
)

