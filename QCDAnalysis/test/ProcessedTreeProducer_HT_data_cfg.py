import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_42_V19::All'
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/data/Run2011A/HT/AOD/PromptReco-v4/000/168/229/14A1D37B-C9A4-E011-A918-003048F1110E.root'
    )
)
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_HT_data.root'))

process.ak5 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.InputTag('ak5PFJets'),
    calojets        = cms.InputTag('ak5CaloJets'),
    ## database entry for the uncertainties ######
    PFPayloadName   = cms.string('AK5PF'),
    CaloPayloadName = cms.string('AK5Calo'),
    ## calojet ID and extender for the JTA #######
    calojetID       = cms.InputTag('ak5JetID'),
    calojetExtender = cms.InputTag('ak5JetExtender'),
    ## set the conditions for good Vtx counting ##
    offlineVertices = cms.InputTag('offlinePrimaryVertices'),
    goodVtxNdof     = cms.double(4), 
    goodVtxZ        = cms.double(24),
    ## rho #######################################
    srcCaloRho      = cms.InputTag('kt6CaloJets','rho'),
    srcPFRho        = cms.InputTag('kt6PFJets','rho'),
    ## preselection cuts #########################
    maxY            = cms.double(3.5), 
    minPFPt         = cms.double(20),
    minPFFatPt      = cms.double(10),
    maxPFFatEta     = cms.double(2.5),
    minCaloPt       = cms.double(20),
    minNPFJets      = cms.int32(4),
    minNCaloJets    = cms.int32(4), 
    minJJMass       = cms.double(-1),
    ## trigger ###################################
    printTriggerMenu = cms.untracked.bool(True),
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring(
      'HLT_HT450_v1','HLT_HT450_v2','HLT_HT450_v3','HLT_HT450_v4','HLT_HT450_v5','HLT_HT450_v6','HLT_HT450_v7','HLT_HT450_v8',
      'HLT_HT500_v1','HLT_HT500_v2','HLT_HT500_v3','HLT_HT500_v4','HLT_HT500_v5','HLT_HT500_v6','HLT_HT500_v7','HLT_HT500_v8',
      'HLT_HT550_v1','HLT_HT550_v2','HLT_HT550_v3','HLT_HT550_v4','HLT_HT550_v5','HLT_HT550_v6','HLT_HT550_v7','HLT_HT550_v8'
    ),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jec services ##############################
    pfjecService    = cms.string('ak5PFL1FastL2L3Residual'),
    calojecService  = cms.string('ak5CaloL1L2L3Residual')
)

############# turn-on the fastjet area calculation needed for the L1Fastjet ##############
############# applied only to PFJets because if CaloJets are re-recoed the JetID map will be lost #####
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(5.0)
process.kt6CaloJets.doRhoFastjet = True
process.kt6CaloJets.Rho_EtaMax = cms.double(5.0)
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(5.0)
process.ak5PFJets.jetPtMin = cms.double(5.0)
############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring(
       'HLT_HT450_v1','HLT_HT450_v2','HLT_HT450_v3','HLT_HT450_v4','HLT_HT450_v5','HLT_HT450_v6','HLT_HT450_v7','HLT_HT450_v8',
       'HLT_HT500_v1','HLT_HT500_v2','HLT_HT500_v3','HLT_HT500_v4','HLT_HT500_v5','HLT_HT500_v6','HLT_HT500_v7','HLT_HT500_v8',
       'HLT_HT550_v1','HLT_HT550_v2','HLT_HT550_v3','HLT_HT550_v4','HLT_HT550_v5','HLT_HT550_v6','HLT_HT550_v7','HLT_HT550_v8'
    ),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.path = cms.Path(process.hltFilter * process.kt6PFJets * process.kt6CaloJets * process.ak5PFJets * process.ak5)


