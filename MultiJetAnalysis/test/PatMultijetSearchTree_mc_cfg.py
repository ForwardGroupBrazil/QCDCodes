import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('PatMultijetSearchTree_mc.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        'file://./patuple_madgraph.root'
        )
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 100
##-------------------- User analyzer  --------------------------------
process.multijets = cms.EDAnalyzer('PatMultijetSearchTree',
    jets    = cms.InputTag('selectedPatJets'),
    met     = cms.InputTag('pfMet'),
    beta    = cms.string('betaAK5PF'),
    etaMAX  = cms.double(2.5),
    ptMIN   = cms.double(30),
    betaMAX = cms.double(1)
)

process.p = cms.Path(process.multijets)

