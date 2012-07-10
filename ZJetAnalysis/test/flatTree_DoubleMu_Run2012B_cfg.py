import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')
from CMGTools.Production.datasetToSource import *
from CMGTools.Common.Tools.applyJSON_cff import *

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
process.source = datasetToSource('cmgtools','/DoubleMu/Run2012B-PromptReco-v1/AOD/PAT_CMG_V5_4_0_runrange_194480-195016', 'cmgTuple_.*root')
data = createDataset('cmgtools','/DoubleMu/Run2012B-PromptReco-v1/AOD/PAT_CMG_V5_4_0_runrange_195017-195396','cmgTuple_.*root', False)
process.source.fileNames.extend(data.listOfGoodFiles())
data = createDataset('cmgtools','/DoubleMu/Run2012B-PromptReco-v1/AOD/PAT_CMG_V5_4_0_runrange_195397-195775','cmgTuple_.*root', False)
process.source.fileNames.extend(data.listOfGoodFiles())
data = createDataset('cmgtools','/DoubleMu/Run2012B-PromptReco-v1/AOD/PAT_CMG_V5_4_0_runrange_195776-195947','cmgTuple_.*root', False)
process.source.fileNames.extend(data.listOfGoodFiles())
data = createDataset('cmgtools','/DoubleMu/Run2012B-PromptReco-v1/AOD/PAT_CMG_V5_4_0_runrange_195948-196509','cmgTuple_.*root', False)
process.source.fileNames.extend(data.listOfGoodFiles())

applyJSON(process,'Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt')

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.zjet = cms.EDAnalyzer('FlatTreeProducer',
    jets          = cms.InputTag('cmgPFJetSel'),
    met           = cms.InputTag('nopuMet'),
    muons         = cms.InputTag('cmgMuonSel'),
    electrons     = cms.InputTag('cmgElectronSel'),
    rho           = cms.InputTag('kt6PFJets','rho'),
    minJetPt      = cms.double(50.0),
    maxJetEta     = cms.double(3.0),
    minMuonPt     = cms.double(20.0),
    minElectronPt = cms.double(20.0)
)

############# hlt filter #########################
process.hltElFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.hltMuFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_Mu17_Mu8_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.p = cms.Path(process.hltMuFilter * process.zjet)
