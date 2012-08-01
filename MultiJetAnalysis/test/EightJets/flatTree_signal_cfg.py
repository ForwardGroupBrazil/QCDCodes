import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree_signal.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
       'file:///afs/cern.ch/work/k/kkousour/private/data/pat_RECO_Pythia6Z2_8Jets_7TeV_GEN_1000_333_100_1_1_1_D6a.root'
        )
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 100
##-------------------- User analyzer  --------------------------------
process.multijets = cms.EDAnalyzer('PatMultijetSearchTree',
    jets    = cms.InputTag('jetExtender','extendedPatJets'),
    met     = cms.InputTag('pfMet'),
    rho     = cms.InputTag('kt6PFJets','rho'),
    etaMAX  = cms.double(2.5),
    ptMIN   = cms.double(20),
    betaMIN = cms.double(0),
    isMC    = cms.untracked.bool(True)
)

process.p = cms.Path(process.multijets)

