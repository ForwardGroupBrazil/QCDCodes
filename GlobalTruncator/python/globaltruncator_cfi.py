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
    RefitIndex    = cms.vint32( 1, 1, 1, 1, 1,
                                2,  2, 2, 2, 2, 2,
                                3,  3, 3, 3, 3, 3,
                                5,  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5
                                ),
    RefitSubIndex = cms.vint32( 1 , 2, 3, 4, 5,
                                -1, 1, 2, 3, 4, 5,
                                -1, 1, 2, 3, 4, 5,
                                -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30
                                ),
    Refits = cms.vstring(
                         
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
                         'picky',
                         'picky1',
                         'picky2',
                         'picky3',
                         'picky4',
                         'picky5',
                         'hit',
                         'hit0',
                         'hit1',
                         'hit2',
                         'hit3',
                         'hit4',
                         'hit5',
                         'hit6',
                         'hit7',
                         'hit8',
                         'hit9',
                         'hit10',
                         'hit11',
                         'hit12',
                         'hit13',
                         'hit14',
                         'hit15',
                         'hit16',
                         'hit17',
                         'hit18',
                         'hit19',
                         'hit20',
                         'hit21',
                         'hit22',
                         'hit23',
                         'hit24',
                         'hit25',
                         'hit26',
                         'hit27',
                         'hit28',
                         'hit29',
                         'hit30',                        
                         ), 
    MuonCollectionLabel = cms.InputTag("globalMuons"),
    RefitterParameters = cms.PSet(
    GlobalTruncRefitter
    )
    
    )
