import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_311_V2::All'
process.load("CondCore.DBCommon.CondDBCommon_cfi")
##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
#############   Define the source file ###############
#process.load('KKousour.QCDAnalysis.MCFileNames_cfi')
#############   Import the HLT filters ###############
process.load('KKousour.QCDAnalysis.hltFilters_cff')
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Spring11/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/GEN-SIM-RECODEBUG/E7TeV_FlatDist10_2011EarlyData_50ns_START311_V1G1-v1/0002/FEA7702A-033E-E011-989F-00215E21DD26.root')
)
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_mc.root'))

process.ak7 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.string('ak7PFJets'),
    calojets        = cms.string('ak7CaloJets'),
    ## monte carlo specific #####################
    genjets         = cms.untracked.string('ak7GenJets'),
    isMCarlo        = cms.untracked.bool(True),
    srcPU           = cms.untracked.InputTag('addPileupInfo'),
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
    processName     = cms.string('REDIGI311X'),
    triggerName     = cms.string('HLT_Jet100U_v3'),
    triggerResults  = cms.InputTag("TriggerResults","","REDIGI311X"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","REDIGI311X"),
    ## jec services ##############################
    pfjecService    = cms.string('ak7PFL1L2L3'),
    calojecService  = cms.string('ak7CaloL1L2L3')
)

process.path = cms.Path(process.ak7)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


