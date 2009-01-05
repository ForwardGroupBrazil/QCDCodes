import FWCore.ParameterSet.Config as cms

paramswitch = cms.EDAnalyzer('ParamSwitch',
    muonsTag = cms.InputTag("muons"),
                              selectionTag = cms.InputTag("AllGlobalMuons"),
)                    
