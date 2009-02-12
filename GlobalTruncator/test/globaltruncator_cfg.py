import FWCore.ParameterSet.Config as cms

processName = "reTrunc"
process = cms.Process(processName)

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( (
    #'file:myTestCopy.root',
#    '/store/user/aeverett/truncStudy220/SingleMuPt1000-step2/step2-SingleMuPt1000_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_9.root',
#'/store/user/aeverett/noteReco220/SingleMuPt100-step2/step2-SingleMuPt100_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_2.root',
#'/store/user/aeverett/truncStudy220FullEvt/SingleMuPt1000-step2//step2-SingleMuPt1000_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_4.root',
    #'file:validationEDM.root',
    'file:step2-SingleMuPt100_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_11.root',
    ))
secFiles.extend((
#'/store/user/aeverett/note220/SingleMuPt100/SingleMuPt100_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_2.root',
    ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.EventContent.EventContent_cff")


process.out = cms.OutputModule(
    "PoolOutputModule",
    process.RECOEventContent,
    outputCommands = cms.untracked.vstring('drop *', "keep *_MEtoEDMConverter_*_"+processName),
    fileName = cms.untracked.string('validationEDM.root')
    )
#process.out.outputCommands.append('keep *_globalTrunc_*_*')
#process.out.outputCommands.append('keep *')
process.outpath = cms.EndPath(process.out)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.categories = ['TruncAnalyzer', 'TrackValidator','RecoMuonValidator']
process.MessageLogger.debugModules = ['truncAnalyzer','recoMuonVTrackAssoc']
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    ),
    TruncAnalyzer = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    ),
    RecoMuonValidator = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    ),
)
process.MessageLogger.cerr = cms.untracked.PSet(
    placeholder = cms.untracked.bool(True)
)

#process.load("Configuration/StandardSequences/SimulationRandomNumberGeneratorSeeds_cff")

process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.MEtoEDMConverter_step = cms.Path(process.MEtoEDMConverter)

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration/StandardSequences/RawToDigi_cff')
#process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.load("UserCode.L3Switches.TPLink_cfi")
#process.simParticle_step = cms.Path(process.TPLink)

#process.load("RecoMuon.Configuration.RecoMuon_cff")
process.load("UserCode.GlobalTruncator.globaltruncator_cfi")

#process.p1 = cms.Path(process.RawToDigi*process.reconstruction)

#process.p1 = cms.Path(process.RawToDigi*process.localreco*process.globalreco) #*process.globalTrunc)

#process.p1 = cms.Path(process.RawToDigi*process.tevMuons*process.globalTrunc)

process.raw2digi_step = cms.Path(process.RawToDigi)

#process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
#process.raw2digi_step = cms.Path(process.RawToDigi)


process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.GlobalTag.globaltag = "IDEAL_V11::All"

#---- Validation stuffs ----#
## Default validation modules
process.load("Configuration.StandardSequences.Validation_cff")
process.validation_step = cms.Path(process.validation)
## Load muon validation modules
#process.recoMuonVMuAssoc.outputFileName = 'validationME.root'
process.load("Validation.RecoMuon.muonValidation_cff")
#process.muonValidation_step = cms.Path(cms.SequencePlaceholder("mix")+process.recoMuonValidation)
process.muonValidation_step = cms.Path(process.recoMuonValidation)

from Validation.RecoMuon.selectors_cff import *
from Validation.RecoMuon.associators_cff import *
process.load("UserCode.GlobalTruncator.truncAnalyzer_cfi")
process.truncAnalyzer.outputFileName = 'validationME2.root'
#process.truncationV_step = cms.Path(cms.SequencePlaceholder("mix")+process.truncAnalyzer)
process.truncationV_step = cms.Path(process.truncAnalyzer)

process.schedule = cms.Schedule(
#    process.p1,
#    process.simParticle_step,
    process.raw2digi_step,

#    process.validation_step,
    process.muonValidation_step,
    process.truncationV_step,
    process.MEtoEDMConverter_step,
    #process.outpath
    )

