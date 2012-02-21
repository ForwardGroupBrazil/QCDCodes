import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('test_flat_signal.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
#process.load('KKousour.MultiJetAnalysis.Axi900sig300_cfi')
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
       '/store/user/kkousour/HT/HT_Run2011A_Aug05-Multijets-PAT/708474a592ff9c61b51f3f4524309977/patuple_multijets_9_1_3EX.root'
        )
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.multijets = cms.EDAnalyzer('PatMultijetSearchTree',
    jets    = cms.InputTag('selectedPatJets'),
    met     = cms.InputTag('pfMet'),
    rho     = cms.InputTag('kt6PFJets','rho'),
    beta    = cms.string('betaAK5PF'),
    etaMAX  = cms.double(2.5),
    ptMIN   = cms.double(30),
    betaMAX = cms.double(1)
)

process.p = cms.Path(process.multijets)

