import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_311_V2::All'
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK7PF'),
            label  = cms.untracked.string('AK7PF')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK7Calo'),
            label  = cms.untracked.string('AK7Calo')
            )
      ), 
      connect = cms.string('sqlite:Jec10V3.db')
)
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/161/311/AE2177F9-DF57-E011-BD15-001D09F26509.root')
)
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_data.root'))

process.ak7 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.InputTag('ak7PFJets'),
    calojets        = cms.InputTag('ak7CaloJets'),
    ## database entry for the uncertainties ######
    PFPayloadName   = cms.string('AK7PF'),
    CaloPayloadName = cms.string('AK7Calo'),
    ## calojet ID and extender for the JTA #######
    calojetID       = cms.InputTag('ak7JetID'),
    calojetExtender = cms.InputTag('ak7JetExtender'),
    ## set the conditions for bood Vtx counting ##
    goodVtxNdof     = cms.double(4), 
    goodVtxZ        = cms.double(24),
    ## preselection cuts ###############
    maxY            = cms.double(3.5), 
    minPFPt         = cms.double(20),
    minCaloPt       = cms.double(20),
    minNPFJets      = cms.int32(1),
    minNCaloJets    = cms.int32(1), 
    ## trigger ##############################
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring('HLT_Jet30_v1','HLT_Jet30_v2',
                                  'HLT_Jet60_v1','HLT_Jet60_v2',
                                  'HLT_Jet80_v1','HLT_Jet80_v2',
                                  'HLT_Jet110_v1','HLT_Jet110_v2',
                                  'HLT_Jet150_v1','HLT_Jet150_v2',
                                  'HLT_Jet190_v1','HLT_Jet190_v2',
                                  'HLT_Jet240_v1','HLT_Jet240_v2',
                                  'HLT_Jet300_v1',
                                  'HLT_Jet370_v1','HLT_Jet370_v2'),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jec services ##############################
    pfjecService    = cms.string('ak7PFL1FastL2L3Residual'),
    calojecService  = cms.string('ak7CaloL1L2L3Residual')
)
############# turn-on the fastjet area calculation needed for the L1Fastjet ##############
############# applied only to PFJets because if CaloJets are re-recoed the JetID map will be lost #####
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(5.0)
process.ak7PFJets.doAreaFastjet = True
process.ak7PFJets.Rho_EtaMax = cms.double(5.0)

process.path = cms.Path(process.kt6PFJets * process.ak7PFJets * process.ak7)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10


