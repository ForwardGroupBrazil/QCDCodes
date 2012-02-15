import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree_mc.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        'file://./patTuple_mc.root'
        )
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('PatVBFTree',
    jets    = cms.InputTag('jetExtender','extendedPatJets'),
    met     = cms.InputTag('pfMet'),
    rho     = cms.InputTag('kt6PFJets','rho'),
    pu      = cms.untracked.string('addPileupInfo'),
    isBKG   = cms.untracked.bool(True),
    btagger = cms.string('combinedSecondaryVertexBJetTags')
)

process.p = cms.Path(process.Hbb)

