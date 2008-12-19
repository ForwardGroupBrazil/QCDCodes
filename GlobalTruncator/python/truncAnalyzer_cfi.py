import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy

from Validation.RecoTrack.TrackingParticleSelectionForEfficiency_cfi import *
truncAnalyzer = cms.EDAnalyzer(
    "TruncAnalyzer",
    MuonServiceProxy,
    TrackingParticleSelectionForEfficiency,

    
    simLabel = cms.InputTag("mergedtruth","MergedTrackTruth"),
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
    ),
    
    
    #
    # Histogram dimensions     #
    #
    nBinRes = cms.untracked.uint32(25),
    
    # pT resolution     #
    minResPt = cms.untracked.double(-0.2),
    maxResPt = cms.untracked.double(0.2),
    
    
    )
