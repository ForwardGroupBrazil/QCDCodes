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
        'file://./patTuple.root'
        )
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('PatVBFTree',
    jets    = cms.InputTag('jetExtender','extendedPatJets'),
    met     = cms.InputTag('pfMet'),
    rho     = cms.InputTag('kt6PFJets','rho'),
    btagger = cms.string('combinedSecondaryVertexBJetTags')
)

process.hlt  = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_QuadJet40_v*','HLT_QuadJet70_v*','HLT_QuadJet80_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.p = cms.Path(process.hlt * process.Hbb)

