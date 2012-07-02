import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')
from CMGTools.Production.datasetToSource import *

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(100000)
        )
process.source = datasetToSource('cmgtools','/DoubleElectron/Run2012A-PromptReco-v1/RECO/PAT_CMG_V5_4_0_runrange_190605-194076/', 'cmgTuple_.*root')

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.zjet = cms.EDAnalyzer('FlatTreeProducer',
    jets        = cms.InputTag('cmgPFJetSel'),
    met         = cms.InputTag('nopuMet'),
    muons       = cms.InputTag('cmgMuonSel'),
    electrons   = cms.InputTag('cmgElectronSel'),
    diMuons     = cms.InputTag('cmgDiMuonSel'),
    diElectrons = cms.InputTag('cmgDiElectronSel'),
    photons     = cms.InputTag('cmgPhotonSel'),
    rho         = cms.InputTag('kt6PFJets','rho'),
    minJetPt    = cms.double(50.0),
    maxJetEta   = cms.double(3.0),
    minPhotonPt = cms.double(50.0),
)

############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.p = cms.Path(process.hltFilter * process.zjet)
