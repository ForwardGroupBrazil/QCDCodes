import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
       'file://./test_vbf_patuple.root'
        )
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1
##-------------------- User analyzer  --------------------------------
process.test = cms.EDAnalyzer('PatTest',
    jets    = cms.InputTag('jetExtender','extendedPatJets'),
    met     = cms.InputTag('pfMet'),
    rho     = cms.InputTag('kt6PFJets','rho')
)

process.testCHS = cms.EDAnalyzer('PatTest',
    jets    = cms.InputTag('jetExtenderCHS','extendedPatJetsCHS'),
    met     = cms.InputTag('pfMet'),
    rho     = cms.InputTag('kt6PFJets','rho')
)

process.hlt  = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_HT550_v*','HLT_HT600_v*','HLT_HT650_v*','HLT_HT700_v*','HLT_HT750_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.p = cms.Path(process.test * process.testCHS)

