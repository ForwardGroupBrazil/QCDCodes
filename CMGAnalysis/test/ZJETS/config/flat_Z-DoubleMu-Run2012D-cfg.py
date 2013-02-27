import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

from CMGTools.Production.datasetToSource import *
from CMGTools.Common.Tools.applyJSON_cff import *

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.source = datasetToSource('cmgtools','/DoubleMu/Run2012D-PromptReco-v1/AOD/PAT_CMG_V5_13_0','cmgTuple_.*.root')
applyJSON(process,'Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt')
############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring(
       'HLT_Mu17_Mu8_v*'
    ),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.zjets = cms.EDAnalyzer('ZJetsFlatTreeProducer',
    jets          = cms.InputTag('cmgPFJetSel'),
    met           = cms.InputTag('cmgPFMETRaw'),
    muons         = cms.InputTag('cmgMuonSel'),
    electrons     = cms.InputTag('cmgElectronSel'),
    rho           = cms.InputTag('kt6PFJets','rho'),
    minJetPt      = cms.double(10.0),
    maxJetEta     = cms.double(5.0),
    minMuonPt     = cms.double(20.0),
    minElectronPt = cms.double(20.0),
    btagger       = cms.string('combinedSecondaryVertexBJetTags'),
    ## trigger ###################################
    triggerAlias     = cms.vstring('HLT_Mu17_Mu8'),
    triggerSelection = cms.vstring(
      'HLT_Mu17_Mu8_v*'
    ),
    triggerConfiguration = cms.PSet(
      hltResults            = cms.InputTag('TriggerResults','','HLT'),
      l1tResults            = cms.InputTag(''),
      daqPartitions         = cms.uint32(1),
      l1tIgnoreMask         = cms.bool(False),
      l1techIgnorePrescales = cms.bool(False),
      throw                 = cms.bool(False)
    )
)

process.p = cms.Path(process.hltFilter * process.zjets)
