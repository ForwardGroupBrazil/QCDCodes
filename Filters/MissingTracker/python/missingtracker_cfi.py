import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer("MissingTracker",
                      glbLabel = cms.untracked.InputTag("globalMuons"),
                      muLabel = cms.untracked.InputTag("muons"),
                      out = cms.string('validationPlots.root'),
                      dirName = cms.string('RecoMuonV/MultiTrack/'),
                      nintHit = cms.int32(75),
                      min = cms.double(-2.5),
                      maxHit = cms.double(75.0),
                      minHit = cms.double(0.0),
                      max = cms.double(2.5),
                      nint = cms.int32(50)
                      )                      
