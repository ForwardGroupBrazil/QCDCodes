import FWCore.ParameterSet.Config as cms

l2l3analyzer = cms.EDAnalyzer(
    'L2L3PtAnalyzer',
    L2Label = cms.InputTag('hltL2Muons:UpdatedAtVtx'),
    L3TkLabel = cms.InputTag('hltL3TkTracksFromL2'),
)
