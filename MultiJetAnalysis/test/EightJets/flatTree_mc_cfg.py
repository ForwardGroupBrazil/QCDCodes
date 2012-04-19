import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

##-------------------- Define the source  ----------------------------
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        'file://./patTuple.root'
        )
)
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.multijets = cms.EDAnalyzer('PatMultijetSearchTree',
    jets    = cms.InputTag('jetExtender','extendedPatJets'),
    met     = cms.InputTag('pfMet'),
    rho     = cms.InputTag('kt6PFJets','rho'),
    etaMAX  = cms.double(2.5),
    ptMIN   = cms.double(30),
    betaMIN = cms.double(0),
    isMC    = cms.untracked.bool(True),
    pu      = cms.untracked.string('addPileupInfo')
)

process.p = cms.Path(process.multijets)

