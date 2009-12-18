import FWCore.ParameterSet.Config as cms

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.GlobalTrackingTools.MuonTrackingRegionCommon_cff import *


globalMatchingAnalyser = cms.EDAnalyzer(
    'GlobalMatchingAnalyser',
    MuonServiceProxy,
    MuonTrackingRegionCommon,
    outputFileName = cms.untracked.string("matchAnalyser.root"),
    trackLabel = cms.InputTag("generalTracks"),
    muonLabel = cms.InputTag("muons"),
)
