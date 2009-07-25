import FWCore.ParameterSet.Config as cms

aodmuontontuple = cms.EDAnalyzer("AODMuonToNtuple",
        out = cms.string("aodtest.root"),
        Muon = cms.InputTag("muons"),
        Jet = cms.InputTag("selectedLayer1Jets"),
        MET = cms.InputTag("selectedLayer1METs"),
        isMC = cms.bool(True),
        CrossSection = cms.double(1233),
        FilterEfficiency = cms.double(1.0),
        TotalNevents = cms.double(10100),
        BDiscriminant = cms.double(0.7)
)
