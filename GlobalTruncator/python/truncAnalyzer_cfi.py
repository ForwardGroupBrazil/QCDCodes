import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy

from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import *

truncAnalyzer = cms.EDAnalyzer(
    "TruncAnalyzer",
    MuonServiceProxy,
    #TrackingParticleSelectionForEfficiency,
    src = cms.InputTag("mergedtruth", "MergedTrackTruth"),
    pdgIdTP = cms.vint32(13, -13),
    tipTP = cms.double(3.5),
    lipTP = cms.double(30.0),
    minHitTP = cms.int32(0),
    ptMinTP = cms.double(0.9),
    minRapidityTP = cms.double(-2.4),
    maxRapidityTP = cms.double(2.4),
    signalOnlyTP = cms.bool(True),
    chargedOnlyTP = cms.bool(True),

    
    #simLabel = cms.InputTag("mergedtruth","MergedTrackTruth"),
    simLabel = cms.InputTag("muonTP"),
    glbMuLabel = cms.InputTag("globalMuons"),
    
    glbMuAssocLabel = cms.InputTag("tpToGlbTrackAssociation"),
    
    outputFileName = cms.untracked.string(''),
    subDir = cms.untracked.string('RecoMuonV/Trunc/'),

    label = cms.VInputTag(
    cms.InputTag("tracker"),
    cms.InputTag("globalMuons"),
    cms.InputTag("global"),
    cms.InputTag("globalTrunc:station1"),
    cms.InputTag("globalTrunc:station2"),
    cms.InputTag("globalTrunc:station3"),
    cms.InputTag("globalTrunc:station4"),
    cms.InputTag("globalTrunc:station5"),
    cms.InputTag("globalTrunc:first"),
    cms.InputTag("globalTrunc:first1"),
    cms.InputTag("globalTrunc:first2"),
    cms.InputTag("globalTrunc:first3"),
    cms.InputTag("globalTrunc:first4"),
    cms.InputTag("globalTrunc:first5"),
    cms.InputTag("globalTrunc:picky"),
    cms.InputTag("globalTrunc:picky1"),
    cms.InputTag("globalTrunc:picky2"),
    cms.InputTag("globalTrunc:picky3"),
    cms.InputTag("globalTrunc:picky4"),
    cms.InputTag("globalTrunc:picky5"),
    cms.InputTag("globalTrunc:picky5"),
    cms.InputTag("globalTrunc:hit"),
    cms.InputTag("globalTrunc:hit0"),
    cms.InputTag("globalTrunc:hit1"),
    cms.InputTag("globalTrunc:hit2"),
    cms.InputTag("globalTrunc:hit3"),
    cms.InputTag("globalTrunc:hit4"),
    cms.InputTag("globalTrunc:hit5"),    
    ),

    minStations = cms.untracked.uint32(4),
    
    #
    # Histogram dimensions     #
    #
    nBinRes = cms.untracked.uint32(100),
    
    # pT resolution     #
    minResPt = cms.untracked.double(-0.2),
    maxResPt = cms.untracked.double(0.2),
    
    
    )
