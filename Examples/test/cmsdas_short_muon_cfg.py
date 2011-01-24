import FWCore.ParameterSet.Config as cms

process = cms.Process("PATMuons")

### standard includes
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.MessageLogger.cerr.FwkReport.reportEvery = 50 
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

### global tag
process.GlobalTag.globaltag = 'START38_V13::All'

### source
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource",
    fileNames = readFiles,
    secondaryFileNames = secFiles,
    )
readFiles.extend( [
    '/store/mc/Fall10/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/START38_V12-v1/0003/00F38C7F-49CB-DF11-9110-00237DF1FFB0.root'
    #'/store/relval/CMSSW_3_8_7/RelValTTbar_Tauola/GEN-SIM-RECO/START38_V13_PU_E7TeV_AVE_2_BX156-v1/0018/F0A1D986-4CFD-DF11-89B5-001A92971B62.root',
    ])

#--- For including HLT DEBUG, usually not available ---
#secFiles.extend( [
#    '/store/relval/CMSSW_3_8_7/RelValTTbar_Tauola/GEN-SIM-DIGI-RAW-HLTDEBUG/START38_V13_PU_E7TeV_AVE_2_BX156-v1/0018/F6D53708-4CFD-DF11-86C9-001A92811716.root',
#    '/store/relval/CMSSW_3_8_7/RelValTTbar_Tauola/GEN-SIM-DIGI-RAW-HLTDEBUG/START38_V13_PU_E7TeV_AVE_2_BX156-v1/0018/E6075587-4CFD-DF11-AF23-0018F3D09680.root',
#    '/store/relval/CMSSW_3_8_7/RelValTTbar_Tauola/GEN-SIM-DIGI-RAW-HLTDEBUG/START38_V13_PU_E7TeV_AVE_2_BX156-v1/0018/868C1C87-4CFD-DF11-90A8-001A928116F8.root',
#    '/store/relval/CMSSW_3_8_7/RelValTTbar_Tauola/GEN-SIM-DIGI-RAW-HLTDEBUG/START38_V13_PU_E7TeV_AVE_2_BX156-v1/0018/066FB883-4CFD-DF11-A27B-002618943984.root',
#    ] )

### number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )


### Adding Trigger Info from TriggerResultsSummary to the PATMuon
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import changeTriggerProcessName, switchOffAmbiguityResolution;

## ==== Change trigger process name ====
#changeTriggerProcessName(process, "HLT2")   # Custom re-run HLT
#changeTriggerProcessName(process, "REDIGI") # Spring10 ReDigi

## ==== Change which track is used for L1 matching ====
#process.muonL1Info.useTrack = 'global'
#process.muonL1Info.useState = 'outermost'
#process.muonMatchHLTL1.useTrack = 'global'
#process.muonMatchHLTL1.useState = 'outermost'

## === Switch off ambiguity resolution in trigger matching ====
## can be a good idea if you're working with events with a single true muon + ghosts
## and you're applying the quality cuts that kill the ghosts only afterwards
switchOffAmbiguityResolution(process)

### Adding Info about the Muon Station involved to the PATMuon
process.load("MuonAnalysis.Examples.muonStations_cfi")
from MuonAnalysis.Examples.muonStations_cfi import addUserData as addStations
addStations(process.patMuonsWithoutTrigger)

### Adding Info about the Muon Hits involved to the PATMuon
# Requires MuonAnalysis/Examples V00-03-00+
process.load("UserCode.Examples.muonHitCount_cfi")
from UserCode.Examples.muonHitCount_cfi import addUserData as addHitCount
addHitCount(process.patMuonsWithoutTrigger)

process.load("UserCode.Examples.muonCosmicMaps_cfi")
from UserCode.Examples.muonCosmicMaps_cfi import addUserData as addCosmicData
addCosmicData(process.patMuonsWithoutTrigger)

### FILTERS for GoodCollision
process.oneGoodVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof >= 4 && abs(z) <= 15 && position.Rho <= 2"),
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)
process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.globalMuFilter = cms.EDFilter("TrackCountFilter",    src = cms.InputTag("globalMuons"), minNumber = cms.uint32(1))
process.recoMuFilter   = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("muons"),       minNumber = cms.uint32(1))

process.preFilter = cms.Sequence(process.oneGoodVertexFilter * process.noScraping)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltL1MuOpen  = hltHighLevel.clone(HLTPaths = ['HLT_L1MuOpen'], )

process.skim_RecoMu            = cms.Path(process.preFilter * process.recoMuFilter  )
process.skim_GlbMu             = cms.Path(process.preFilter * process.globalMuFilter)
process.skim_RecoMuHLTL1MuOpen = cms.Path(process.preFilter * process.recoMuFilter   * process.hltL1MuOpen)

### Add MC classification by hits
process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")
from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
addClassByHits(process.patMuonsWithoutTrigger, extraInfo=True)
## add two other classifications with preselection, to get different ghosts
process.classByHitsTMA = process.classByHitsTM.clone(
    muonPreselection = cms.string("muonID('TrackerMuonArbitrated')")
)
process.classByHitsGlbPT = process.classByHitsGlb.clone(
    muonPreselection = cms.string("muonID('GlobalMuonPromptTight')")
)
process.muonClassificationByHits.replace(process.classByHitsTM,  process.classByHitsTM  + process.classByHitsTMA)
process.muonClassificationByHits.replace(process.classByHitsGlb, process.classByHitsGlb + process.classByHitsGlbPT)
addClassByHits(process.patMuonsWithoutTrigger, labels=["classByHitsTMA","classByHitsGlbPT"], extraInfo=False)

#--- For including HLT DEBUG, usually not available ---
## process.load("MuonAnalysis.MuonAssociators.triggerMatcherToHLTDebug_cfi")
## from MuonAnalysis.MuonAssociators.triggerMatcherToHLTDebug_cfi import addUserData as addTriggerMatchToHLTDebug
## addTriggerMatchToHLTDebug(process.patMuonsWithoutTrigger) 

#--- For including the Muon Shower Block, requires reReco in 3_8_7---
## from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
## process.muonShower = cms.EDProducer("MuonShowerInformationProducer",
##     MuonServiceProxy,
##     muonCollection = cms.InputTag("muons"),
##     trackCollection = cms.InputTag("generalTracks"),
##     ShowerInformationFillerParameters = cms.PSet(
##       MuonServiceProxy,
##
##       InputTrackCollection = cms.InputTag("generalTracks"),
##       InputMuonCollections = cms.VInputTag(cms.InputTag("globalMuons"), cms.InputTag("muons")),
##       DTRecSegmentLabel = cms.InputTag("dt1DRecHits"),
##       CSCRecSegmentLabel = cms.InputTag("csc2DRecHits"),
##       RPCRecSegmentLabel = cms.InputTag("rpcRecHits"),
##       DT4DRecSegmentLabel = cms.InputTag("dt4DSegments"),
##       CSCSegmentLabel = cms.InputTag("cscSegments"),
##
##       TrackerRecHitBuilder = cms.string('WithTrackAngle'),
##       MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
##      )
##
## )
## process.fillShower = cms.EDAnalyzer('TestShower',
##    inputMuonShowerInformationValueMap = cms.untracked.InputTag("muonShower")
## )
## process.load("UserCode.Examples.muonShowerMaps_cfi")
## from UserCode.Examples.muonShowerMaps_cfi import addUserData as addShowerInfo
## addShowerInfo(process.patMuonsWithoutTrigger)
## process.TFileService = cms.Service("TFileService",
##     fileName = cms.string('mc_zmumu.root')
## )

process.p = cms.Path(
    process.preFilter  +
    process.recoMuFilter  +
    ##     (    process.muonShower*
    ##          process.fillShower+
    ##          process.muonShowerInfoLoader ) +
    process.muonClassificationByHits +
    process.muonStations +
    process.muonHitCounts +
    process.cosmicCompatibilityLoader  +
    ##    process.triggerMatcherToHLTDebug + 
    process.patMuonsWithTriggerSequence 
    )

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('tupleMC.root'),
    outputCommands = cms.untracked.vstring('drop *',
        'keep patMuons_patMuonsWithTrigger__*',                    ## All PAT muons
        'keep edmTriggerResults_TriggerResults_*_*',               ## HLT info, per path (cheap)
        'keep l1extraL1MuonParticles_l1extraParticles_*_*',        ## L1 info (cheap)
        #--- Other things you might want ---
        'keep *_offlinePrimaryVertices__*',                  ## 
        'keep *_offlineBeamSpot__*',                         ##
        'keep recoTrackExtras_standAloneMuons_*_*',          ## track states at the muon system, to which patMuons sta tracks point (useful if you want variables of the innermost or outermost state)
        'keep TrackingRecHitsOwned_standAloneMuons_*_*',     ## muon rechits, to compute things like number of stations
        #--- Other things which can be expensive but you might still want them ---
        'keep triggerTriggerEvent_hltTriggerSummaryAOD__*',  ## HLT 4-vectors
        'keep *_patTrigger__*',                              ## As above, in PAT format
        #--- HLT DEBUG, usually not available ---
        #'keep recoTracks_hltL2Muons_*_*',
        #'keep recoTracks_hltL3Muons_*_*',
        #'keep recoTrackExtras_hltL2Muons_*_*',
        #'keep recoTrackExtras_hltL3Muons_*_*',
        #'keep recoRecoChargedCandidates_hltL2MuonCandidates__*',
        #'keep recoRecoChargedCandidates_hltL3MuonCandidates__*',
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring("skim_RecoMu") ),
)
process.e = cms.EndPath(process.out)

