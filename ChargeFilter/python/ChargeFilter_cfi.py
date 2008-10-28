import FWCore.ParameterSet.Config as cms

process.chargeFltr = cms.EDFilter(
    "ChargeFilter",
    trkMuAssocLabel =cms.InputTag("tpToGlbTrackAssociation"),
    charge = cms.int32(1),
    simLabel = cms.InputTag("muonTP"),
    trkLabel = cms.InputTag("globalMuons")
    )
