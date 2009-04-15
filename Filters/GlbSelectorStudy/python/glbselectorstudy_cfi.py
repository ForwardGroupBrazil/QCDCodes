import FWCore.ParameterSet.Config as cms
from TrackingTools.TrackRefitter.TracksToTrajectories_cff import *
from SimTracker.TrackHistory.TrackClassifier_cff import *

glbSelStudy = cms.EDAnalyzer(
    'GlbSelectorStudy',

    trackClassifier,

    simLabel = cms.InputTag("mergedtruth","MergedTrackTruth"),
    trkMuLabel = cms.InputTag("generalTracks"),
    staMuLabel = cms.InputTag("standAloneMuons:UpdatedAtVtx"),
    glbMuLabel = cms.InputTag("globalMuons"),
    muonLabel = cms.InputTag("muons"),

    trkMuAssocLabel = cms.InputTag("TrackAssociatorByHits"),
    staMuAssocLabel = cms.InputTag("TrackAssociatorByDeltaR"),
    glbMuAssocLabel = cms.InputTag("TrackAssociatorByDeltaR"),
    doAssoc = cms.untracked.bool(True),

    TrackTransformer = cms.PSet(DoPredictionsOnly = cms.bool(False),
                                Fitter = cms.string('KFFitterForRefitInsideOut'),
                                #        TrackerRecHitBuilder = cms.string('WithTrackAngleAndTemplate'),
                                TrackerRecHitBuilder = cms.string('WithTrackAngle'),
                                Smoother = cms.string('KFSmootherForRefitInsideOut'),
                                MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
                                RefitDirection = cms.string('alongMomentum'),
                                RefitRPCHits = cms.bool(True),
                                Propagator = cms.string('SmartPropagatorAnyRKOpposite')
                                ),
    
    tpSelector = cms.PSet(
       src = cms.InputTag("mergedtruth", "MergedTrackTruth"),
       pdgId = cms.vint32(13, -13),
       tip = cms.double(3.5),
       lip = cms.double(30.0),
       minHit = cms.int32(0),
       ptMin = cms.double(0.9),
       minRapidity = cms.double(-2.4),
       maxRapidity = cms.double(2.4),
       signalOnly = cms.bool(True),
       chargedOnly = cms.bool(True)
       ),
    
    IDconverttoBinNum = cms.PSet(
    ranges = cms.VPSet(cms.PSet(
    pdgIDs = cms.vint32(0),
    label = cms.string('ID=0')
    ), 
        cms.PSet(
            pdgIDs = cms.vint32(211, -211),
            label = cms.string('#pi+/-')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(321, -321, 130, -130),
            label = cms.string('K')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(411, -411, 421, -421, 431, 
                -431),
            label = cms.string('D')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(521, -521, 511, -511, 531, 
                -531),
            label = cms.string('B')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(5122, -5122),
            label = cms.string('#Lambda_{b}')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(443, -443),
            label = cms.string('J/#Psi')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(553, -553, 100553, -100553, 200553, 
                -200553, 300553, -300553),
            label = cms.string('#Upsilon(nS)')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(24, -24),
            label = cms.string('W^{+/-}')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(23),
            label = cms.string('Z^{0}')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(15, -15),
            label = cms.string('#tau^{+/-}')
        ), 
        cms.PSet(
            pdgIDs = cms.vint32(13, -13),
            label = cms.string('#mu^{+/-}')
        )),
    title = cms.string('a set of PDGids')
)
                             
)
