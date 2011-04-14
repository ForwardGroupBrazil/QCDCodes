import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_311_V2::All'
process.load("CondCore.DBCommon.CondDBCommon_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Jec11V0_AK7PF'),
            label  = cms.untracked.string('AK7PF')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Jec11V0_AK7Calo'),
            label  = cms.untracked.string('AK7Calo')
            )
      ), 
      connect = cms.string('sqlite:Jec11V0.db')
)
##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
#############   Define the source file ###############
process.load('KKousour.QCDAnalysis.SkimFileNames_cfi')
#############   Import the HLT filters ###############
process.load('KKousour.QCDAnalysis.hltFilters_cff')
#############   Define the source file ###############
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/data/kkousour/JetAOD_9_1_Bh3.root')
#)
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_data.root'))

process.ak7 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.string('ak7PFJets'),
    calojets        = cms.string('ak7CaloJets'),
    ## database entry for the uncertainties ######
    PFPayloadName   = cms.string(''),
    CaloPayloadName = cms.string(''),
    ## calojet ID and extender for the JTA #######
    calojetID       = cms.string('ak7JetID'),
    calojetExtender = cms.string('ak7JetExtender'),
    ## set the conditions for bood Vtx counting ##
    goodVtxNdof     = cms.double(4), 
    goodVtxZ        = cms.double(24),
    ## number of jets to be stored ###############
    minPFPt         = cms.double(20),
    minCaloPt       = cms.double(20),
    minNPFJets      = cms.int32(1),
    minNCaloJets    = cms.int32(1), 
    ## trigger ##############################
    processName     = cms.string('HLT'),
    triggerName     = cms.string('HLT_Jet110_v1'),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jec services ##############################
    pfjecService    = cms.string('ak7PFL1L2L3'),
    calojecService  = cms.string('ak7CaloL1L2L3')
)

process.path = cms.Path(process.hlt110v1 * process.ak7)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 100


