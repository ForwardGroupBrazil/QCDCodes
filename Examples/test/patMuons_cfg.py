import FWCore.ParameterSet.Config as cms

process = cms.Process("PATMuon")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'START3X_V26A::All'
             
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

#process.load("SimGeneral.MixingModule.mixNoPU_cfi")
#process.load("SimGeneral.TrackingAnalysis.trackingParticlesNoSimHits_cfi")    # On RECO
#process.load("SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi")  # On RECO
##process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")            # On RAW+RECO
##process.load("SimMuon.MCTruth.MuonAssociatorByHitsESProducer_cfi")           # On RAW+RECO


from Configuration.EventContent.EventContent_cff import *
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/user/aeverett/CMSSW_3_4_1/SingleMuPt0_500/aeverett/SingleMuPt0_500_CMSSW_3_4_1_step1/SingleMuPt0_500_CMSSW_3_4_1_step2/bd7cda0b32f61291da1aab637c754773//step2_99.root'
    #'/store/user/aeverett/CMSSW_3_4_1/SinglePiPt2_200/aeverett/SinglePiPt2_200_CMSSW_3_4_1_step1/SinglePiPt2_200_CMSSW_3_4_1_step2/e90b8558827cdd944d7d7455b9e9feec//step2_99_1.root'
    ),
    # secondaryFileNames = cms.untracked.vstring(
    # '/store/user/aeverett/CMSSW_3_4_1/SinglePiPt2_200/aeverett/SinglePiPt2_200_CMSSW_3_4_1_step1/SinglePiPt2_200_CMSSW_3_4_1_step1/85a90b0814831c7eda376193ed517e95/SinglePiPt2_200_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_3_1.root'
    # ),
    inputCommands = RECOSIMEventContent.outputCommands,            # keep only RECO out of RAW+RECO, for tests
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),  # keep only RECO out of RAW+RECO, for tests
    )

process.maxEvents =cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   =cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.bscFilter = hltLevel1GTSeed.clone(L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)'))
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
process.countCollisionEvents = cms.EDProducer("EventCountProducer")
process.preFilter = cms.Sequence(process.bscFilter * process.oneGoodVertexFilter * process.noScraping * process.countCollisionEvents )

process.mergedTruth = cms.EDProducer(
    "GenPlusSimParticleProducer",
    src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
    setStatus     = cms.int32(5),             # set status = 8 for GEANT GPs
    # particleTypes = cms.vstring("mu+"),
    filter        = cms.vstring("pt > 0.0"),  # just for testing (optional)
    genParticles   = cms.InputTag("genParticles") # original genParticle list
    )
process.genMuons = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *  ",                     # this is the default
    "++keep abs(pdgId) = 13",        # keep muons and their parents
    "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
    )
    )

process.load("PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi")

process.filter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("patMuons"),
    minNumber = cms.uint32(1),
    )

import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.patMuons = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone(
    muonSource = 'muons',
    # embed the tracks, so we don't have to carry them around
    embedTrack          = True,
    embedCombinedMuon   = True,
    embedStandAloneMuon = True,
    # then switch off some features we don't need
    #addTeVRefits = False, ## <<--- this doesn't work. PAT bug ??
    embedPickyMuon = False,
    embedTpfmsMuon = False, 
    userIsolation = cms.PSet(),   # no extra isolation beyond what's in reco::Muon itself
    isoDeposits = cms.PSet(), # no heavy isodeposits
    addGenMatch = True,       
    embedGenMatch = True,
)

### Adding MCtruth Info to the PATMuon
## from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_8E29_cff import addMCinfo
addMCinfo(process)

### Adding Info about the MC truth matching
# Requires MuonAnalysis/MuonAssociators
# can run on SIM-RECO as well as HLTDEBUG
process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")
from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
addClassByHits(process.patMuons, extraInfo=True)

### Adding Info about the Muon Station involved to the PATMuon
# Requires MuonAnalysis/Examples V00-03-00+
process.load("MuonAnalysis.Examples.muonStations_cfi")
from MuonAnalysis.Examples.muonStations_cfi import addUserData as addStations
addStations(process.patMuons)

### Adding Info about the Muon Hits involved to the PATMuon
# Requires UserCode/Examples and MuonAnalysis/Examples
process.load("UserCode.Examples.muonHitCount_cfi")
from UserCode.Examples.muonHitCount_cfi import addUserData as addHitCount
addHitCount(process.patMuons)




process.go = cms.Path(
#    process.preFilter +
    ( process.mergedTruth *
      process.genMuons    *
      process.muonMatch   +
      process.muonClassificationByHits +
      process.muonStations +
      process.muonHitCounts
      ) *
    process.patMuons
    * process.filter
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("pat_MuPat.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_mergedTruth_*_*",
        "keep *_patMuons_*_*",
        "keep *_genMuons_*_*",
        "keep *_genParticles_*_*",
        "keep *_countCollisionEvents_*_*",
        "keep recoTrackExtras_standAloneMuons_*_*",          ## track states at the muon system, used both by patMuons and standAloneMuons
        "keep recoTracks_standAloneMuons__*",                ## bare standalone muon tracks, using standalone muon momentum (without BS constraint)
        "keep edmTriggerResults_*_*_HLT",                    ## 
        "keep l1extraL1MuonParticles_l1extraParticles_*_*",  ## 
        "keep *_offlinePrimaryVertices__*",                  ## 
        "keep *_offlineBeamSpot__*",                         ##
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring("go")),
)
process.end = cms.EndPath(process.out)

#process.mergedtruthNoSimHits.vertexDistanceCut = 100.0
