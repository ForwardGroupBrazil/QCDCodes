import FWCore.ParameterSet.Config as cms

muontontuple = cms.EDAnalyzer("MuonToNtuple",
        out = cms.string("test.root"),
        Muon = cms.InputTag("selectedLayer1Muons"),
        Jet = cms.InputTag("selectedLayer1Jets"),
        MET = cms.InputTag("selectedLayer1METs"),
        isMC = cms.bool(True),
        CrossSection = cms.double(1233),
        FilterEfficiency = cms.double(1.0),
        TotalNevents = cms.double(10100),
        BDiscriminant = cms.double(0.7)
)
