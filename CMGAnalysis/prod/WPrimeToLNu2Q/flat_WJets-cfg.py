import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

from CMGTools.Production.datasetToSource import *

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.source = datasetToSource('cmgtools','/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_13_0','cmgTuple_.*.root')

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.wprime = cms.EDAnalyzer('WPrimeFlatTreeProducer',
    jets          = cms.InputTag('cmgPFJetSel'),
    met           = cms.InputTag('cmgPFMETRaw'),
    muons         = cms.InputTag('cmgMuonSel'),
    electrons     = cms.InputTag('cmgElectronSel'),
    rho           = cms.InputTag('kt6PFJets','rho'),
    minJetPt      = cms.double(20.0),
    maxJetEta     = cms.double(2.5),
    minMuonPt     = cms.double(20.0),
    minElectronPt = cms.double(20.0),
    btagger       = cms.string('combinedSecondaryVertexBJetTags'),
    pu            = cms.untracked.string('addPileupInfo'),
    ## trigger ###################################
    triggerAlias     = cms.vstring('HLT_Mu17_Mu8','HLT_Ele17_Ele8'),
    triggerSelection = cms.vstring(
      'HLT_Mu17_Mu8_v*',
      'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*'
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

process.p = cms.Path(process.wprime)
