import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_R_42_V12::All'
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.Geometry_cff')
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Summer11/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6/AODSIM/PU_S3_START42_V11-v2/0004/FA6FEF7F-8E7E-E011-AC35-001A92971B78.root')
)

process.filter = cms.EDFilter('MultiGenJetFilter',
    jets     = cms.InputTag('ak5GenJets'), 
    minPt    = cms.double(30),
    minNjets = cms.int32(4)
)

process.skimPath = cms.Path(process.filter)

#############   output module ########################
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_ak5GenJets_*_*',
        'keep GenEventInfoProduct_generator_*_*'
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('skimPath')),
    fileName = cms.untracked.string('testFilter.root')
)

process.p = cms.EndPath(process.out)


