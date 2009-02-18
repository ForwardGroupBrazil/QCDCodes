import FWCore.ParameterSet.Config as cms

paramswitch = cms.EDAnalyzer('ParamSwitch',
                             muonsTag = cms.InputTag("muons"),
                             selectionTag = cms.InputTag("AllGlobalMuons"),
                             chi2Cut = cms.double(3.0),
                             hitCut = cms.int32(3),
                             minPtCut = cms.double(80.0),
                             maxPtCut = cms.double(120.0),
                             chiFlag = cms.bool(True),
                             hitFlag = cms.bool(True),
                             ptFlag = cms.bool(False),
                             
)                    
