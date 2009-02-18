import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy

missingtracker = cms.EDAnalyzer("MissingTracker",
                      MuonServiceProxy,
                      glbLabel = cms.untracked.InputTag("globalMuons"),
                      muLabel = cms.untracked.InputTag("muons"),
                      out = cms.string('validationPlots.root'),
                      dirName = cms.string('RecoMuonV/MultiTrack/'),
                      nintHit = cms.int32(201),
                      min = cms.double(-2.5),
                      maxHit = cms.double(100.5),
                      minHit = cms.double(-100.5),
                      max = cms.double(2.5),
                      nint = cms.int32(50),

                      RefitRPCHits = cms.bool(True),
                      TrackerRecHitBuilder = cms.string('WithTrackAngle'),
                      MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
                      TrackerPropagator = cms.string('SteppingHelixPropagatorAny'),
                      
                      )                      
