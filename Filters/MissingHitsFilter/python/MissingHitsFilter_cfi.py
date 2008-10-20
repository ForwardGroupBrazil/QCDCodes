import FWCore.ParameterSet.Config as cms

missingHitsFilter = cms.EDFilter("MissingHitsFilter",
                                   muLabel = cms.InputTag("muons"),
                                   hitCut = cms.int32(3),
                                   hitFraction = cms.double(0.20)

    )
