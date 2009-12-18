import FWCore.ParameterSet.Config as cms

globalMatchingAnalyser = cms.EDAnalyzer(
    'GlobalMatchingAnalyser',
    outputFileName = cms.untracked.string("matchAnalyser.root"),
    trackLabel = cms.InputTag("generalTracks"),
    muonLabel = cms.InputTag("muons"),
)
