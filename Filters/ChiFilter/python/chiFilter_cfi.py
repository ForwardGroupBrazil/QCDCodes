import FWCore.ParameterSet.Config as cms

chiFilter = cms.EDFilter(
    "ChiFilter",
    muLabel = cms.InputTag("muons"),
    chi2Cut = cms.double(1.0),
    )
