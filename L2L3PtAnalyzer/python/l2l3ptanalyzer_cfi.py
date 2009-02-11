import FWCore.ParameterSet.Config as cms

l2l3analyzer = cms.EDAnalyzer(
    'L2L3PtAnalyzer',
    out = cms.string(''),
    L2Label = cms.InputTag('hltL2Muons:UpdatedAtVtx'),
    L3TkLabel = cms.InputTag('hltL3TkTracksFromL2'),
    TPLabel = cms.InputTag('mergedtruth:MergedTrackTruth'),
    trkMuAssocLabel = cms.InputTag('TrackAssociatorByHits'),
    ignoremissingl2collection=cms.untracked.bool(True),
    ignoremissingl3tkcollection=cms.untracked.bool(True),
)
