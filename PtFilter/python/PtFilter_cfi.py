import FWCore.ParameterSet.Config as cms

ptFilter = cms.EDFilter("PtFilter",
                                 muLabel = cms.InputTag("muons"),
                                 MinPtCut = cms.double(0.0),
                                 MaxPtCut = cms.double(500.0),
                                 )
