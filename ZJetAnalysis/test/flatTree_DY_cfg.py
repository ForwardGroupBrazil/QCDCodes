import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')
from CMGTools.Production.datasetToSource import *

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
process.source = datasetToSource('cmgtools','/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v2/AODSIM/V5/PAT_CMG_V5_4_0_NewType1MET/','cmgTuple_.*root')

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

process.p = cms.Path(process.zjet)
