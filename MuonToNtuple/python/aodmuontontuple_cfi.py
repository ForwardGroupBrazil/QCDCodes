import FWCore.ParameterSet.Config as cms
from RecoMuon.GlobalTrackingTools.GlobalMuonRefitter_cff import *
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from SimTracker.TrackHistory.TrackClassifier_cff import *

aodmuontontuple = cms.EDAnalyzer("AODMuonToNtuple",
    MuonServiceProxy,

    trackClassifier,

    RefitterParameters = cms.PSet(
    GlobalMuonRefitter,
    ),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(100000.0),

                                 TrackTransformer = cms.PSet(DoPredictionsOnly = cms.bool(False),
                                                             Fitter = cms.string('KFFitterForRefitInsideOut'),
                                                             #        TrackerRecHitBuilder = cms.string('WithTrackAngleAndTemplate'),
                                                             TrackerRecHitBuilder = cms.string('WithTrackAngle'),
                                                             Smoother = cms.string('KFSmootherForRefitInsideOut'),
                                                             MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
                                                             RefitDirection = cms.string('insideOut'),
                                                             RefitRPCHits = cms.bool(True),
                                                             Propagator = cms.string('SmartPropagatorAnyRKOpposite')
                                                             ),

                                 
        out = cms.string("aodtest.root"),
        Muon = cms.InputTag("muons"),
        Jet = cms.InputTag("selectedLayer1Jets"),
        MET = cms.InputTag("layer1METs"),
        isMC = cms.bool(True),
        CrossSection = cms.double(1233),
        FilterEfficiency = cms.double(1.0),
        TotalNevents = cms.double(10100),
        BDiscriminant = cms.double(0.7)
)
