import FWCore.ParameterSet.Config as cms

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *
from RecoMuon.GlobalTrackingTools.GlobalTrajectoryBuilderCommon_cff import *
from UserCode.GlobalTruncator.GlobalTruncRefitter_cff import *
globalTrunc = cms.EDProducer(
    'GlobalTruncator',
    MuonTrackLoaderForGLB,
    #    InputTag MuonCollectionLabel = standAloneMuons:UpdatedAtVtx
    MuonServiceProxy,
    RefitIndex    = cms.vint32( 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
    RefitSubIndex = cms.vint32(-1,-1, 1, 2, 3, 4, 5,-1, 1, 2, 3, 4, 5),
    Refits = cms.vstring('tracker',
                         'global',
                         'station1',
                         'station2',
                         'station3',
                         'station4',
                         'station5',
                         'first',
                         'first1',
                         'first2',
                         'first3',
                         'first4',
                         'first5',
                         ), 
    MuonCollectionLabel = cms.InputTag("globalMuons"),
    RefitterParameters = cms.PSet(
    GlobalTruncRefitter
    )
    
    )
