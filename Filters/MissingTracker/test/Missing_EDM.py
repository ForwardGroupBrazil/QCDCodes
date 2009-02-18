import FWCore.ParameterSet.Config as cms

process = cms.Process("EDMtoMEConvert")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories = ['MissingTracker',]
process.MessageLogger.debugModules = ['*']
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    TrackAssociator = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    MissingTracker = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    )
)
process.MessageLogger.cerr = cms.untracked.PSet(
    placeholder = cms.untracked.bool(True)
)


process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.load("DQMServices.Components.DQMEnvironment_cfi")
#process.load("Configuration.StandardSequences.FakeConditions_cff")
#process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "IDEAL_V11::All"
process.load("Validation.Configuration.postValidation_cff")
process.load("Validation.RecoMuon.PostProcessorHLT_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    "PoolSource",
#    fileNames = cms.untracked.vstring("file:reVal-step2-SingleMuPt0_500_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_47.root",)

    fileNames = cms.untracked.vstring("file:reVal-step2-SingleMuPt100_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1.root",)
    #fileNames = cms.untracked.vstring("file:outputFileName.root",)
    
#    fileNames = cms.untracked.vstring("file:reVal-step2-SingleMuPt200_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1.root",)
#    fileNames = cms.untracked.vstring("file:reVal-step2-SingleMuPt500_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_44.root",)
    )

process.DQMStore.referenceFileName = ""
process.DQMStore.collateHistograms = False

process.dqmSaver.convention = "Offline"
#Settings equivalent to 'RelVal' convention:
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
#End of 'RelVal convention settings
process.dqmSaver.workflow = "/GlobalValidation/Test/RECO"

process.load("UserCode.MissingTracker.missingtracker_cfi")
process.missingtracker.nintHit=201
process.missingtracker.maxHit = 100.5
process.missingtracker.minHit=-100.5

process.p1 = cms.Path(process.EDMtoMEConverter*
                      process.missingtracker*
                      process.postValidation*
#                      process.recoMuonPostProcessorsHLT*
                      process.dqmSaver)
