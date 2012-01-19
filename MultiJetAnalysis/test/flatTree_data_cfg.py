import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        '/store/user/kkousour/HT/HT_Run2011A_Aug05-Multijets-PAT/708474a592ff9c61b51f3f4524309977/patuple_multijets_100_1_1Ke.root'
        )
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.multijets = cms.EDAnalyzer('PatMultijetSearchTree',
    jets    = cms.InputTag('selectedPatJets'),
    met     = cms.InputTag('pfMet'),
    beta    = cms.string('betaAK5PF'),
    etaMAX  = cms.double(2.5),
    ptMIN   = cms.double(30),
    betaMAX = cms.double(1)
)

process.hlt  = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_HT550_v*','HLT_HT600_v*','HLT_HT650_v*','HLT_HT700_v*','HLT_HT750_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.p = cms.Path(process.hlt * process.multijets)

