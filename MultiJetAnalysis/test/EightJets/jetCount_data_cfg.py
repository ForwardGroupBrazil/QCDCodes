import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('jetCount.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
#process.source = cms.Source("PoolSource",
#        fileNames = cms.untracked.vstring(
#       'file://./patTuple.root'
#        )
#)
process.load('KKousour.MultiJetAnalysis.HT_Run2011A_Aug05_cfi')
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.jetCount = cms.EDAnalyzer('PatMultijetCounter',
    jets    = cms.InputTag('jetExtender','extendedPatJets'),
    etaMAX  = cms.double(2.5),
    ptMIN   = cms.double(30),
    htMIN   = cms.double(750),
    betaMIN = cms.double(0)
)

process.hlt  = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_HT550_v*','HLT_HT600_v*','HLT_HT650_v*','HLT_HT700_v*','HLT_HT750_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.p = cms.Path(process.hlt * process.jetCount)

