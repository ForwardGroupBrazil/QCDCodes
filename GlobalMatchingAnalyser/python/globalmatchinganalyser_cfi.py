import FWCore.ParameterSet.Config as cms

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.GlobalTrackingTools.MuonTrackingRegionCommon_cff import *
from RecoMuon.GlobalTrackingTools.GlobalMuonTrackMatcher_cff import *

globalMatchingAnalyser = cms.EDAnalyzer(
    'GlobalMatchingAnalyser',
    MuonServiceProxy,
    MuonTrackingRegionCommon,
    GlobalMuonTrackMatcher,
    trackLabel = cms.InputTag("generalTracks"),
    muonLabel = cms.InputTag("muons"),
    useAll = cms.int32(1)
)
