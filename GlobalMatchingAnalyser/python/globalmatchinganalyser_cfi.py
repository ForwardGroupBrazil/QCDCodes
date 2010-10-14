import FWCore.ParameterSet.Config as cms

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.GlobalTrackingTools.MuonTrackingRegionCommon_cff import *
from RecoMuon.GlobalTrackingTools.GlobalMuonTrackMatcher_cff import *

globalMatchingAnalyser = cms.EDAnalyzer(
    'GlobalMatchingAnalyser',
    MuonServiceProxy,
    MuonTrackingRegionCommon,
    GlobalMuonTrackMatcher,
    trackLabel = cms.InputTag("generalTracks"),
    muonLabel = cms.InputTag("muons"),
    staLabel = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
    useAll = cms.int32(1),
    classification = cms.InputTag("classByHitsGlb"),
    TrackerRecHitBuilder = cms.string('WithTrackAngle'),
    MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
    TrackerPropagator = cms.string('SteppingHelixPropagatorAny'),
    RefitRPCHits = cms.bool(True),
    GlbRefitterParameters = cms.PSet(
       DTRecSegmentLabel = cms.InputTag("dt4DSegments"),
       CSCRecSegmentLabel = cms.InputTag("cscSegments"),
       
       MuonHitsOption = cms.int32(1),
       PtCut = cms.double(1.0),
       Chi2ProbabilityCut = cms.double(30.0),
       Chi2CutCSC = cms.double(150.0),
       Chi2CutDT = cms.double(10.0),
       Chi2CutRPC = cms.double(1.0),
       HitThreshold = cms.int32(1),
       
       Fitter = cms.string('GlbMuKFFitter'),
       Propagator = cms.string('SmartPropagatorAnyRK'),
       TrackerRecHitBuilder = cms.string('WithTrackAngle'),
       MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
       DoPredictionsOnly = cms.bool(False),
       RefitDirection = cms.string('insideOut'),
       PropDirForCosmics = cms.bool(False),
       RefitRPCHits = cms.bool(True),
       
       # muon station to be skipped
       SkipStation		= cms.int32(-1),
       
       # PXB = 1, PXF = 2, TIB = 3, TID = 4, TOB = 5, TEC = 6
       TrackerSkipSystem	= cms.int32(-1),
       
       # layer, wheel, or disk depending on the system
       TrackerSkipSection	= cms.int32(-1)
       ),
    TrackTransformer = cms.PSet(
       Fitter = cms.string('KFFitterForRefitInsideOut'),
       TrackerRecHitBuilder = cms.string('WithTrackAngle'),
       Smoother = cms.string('KFSmootherForRefitInsideOut'),
       MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
       RefitDirection = cms.string('alongMomentum'),
       RefitRPCHits = cms.bool(True),
       DoPredictionsOnly = cms.bool(False),
       Propagator = cms.string('SmartPropagatorAnyRKOpposite') #SmartPropagatorAnyRK
       ),

    
    )
