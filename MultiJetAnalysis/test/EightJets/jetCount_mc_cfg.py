import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('jetCount.root'))

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
process.jetCount = cms.EDAnalyzer('PatMultijetCounter',
    jets    = cms.InputTag('jetExtender','extendedPatJets'),
    etaMAX  = cms.double(2.5),
    ptMIN   = cms.double(30),
    htMIN   = cms.double(750),
    betaMIN = cms.double(0)
)

process.p = cms.Path(process.jetCount)

