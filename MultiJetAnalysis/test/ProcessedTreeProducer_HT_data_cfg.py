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
    input = cms.untracked.int32(-1)
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/data/Run2011B/Jet/AOD/PromptReco-v1/000/178/079/1CDC5E5C-F4F2-E011-BE17-0019B9F70607.root'
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
    minPFPt         = cms.double(25),
    minPFFatPt      = cms.double(10),
    maxPFFatEta     = cms.double(2.5),
    minCaloPt       = cms.double(25),
    minNPFJets      = cms.int32(8),
    minNCaloJets    = cms.int32(8), 
    minJJMass       = cms.double(-1),
    ## trigger ###################################
    printTriggerMenu = cms.untracked.bool(True),
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring(
      'HLT_HT450_v1','HLT_HT450_v2','HLT_HT450_v3','HLT_HT450_v4','HLT_HT450_v5','HLT_HT450_v6','HLT_HT450_v7','HLT_HT450_v8','HLT_HT450_v9','HLT_HT450_v10','HLT_HT450_v11',
      'HLT_HT500_v1','HLT_HT500_v2','HLT_HT500_v3','HLT_HT500_v4','HLT_HT500_v5','HLT_HT500_v6','HLT_HT500_v7','HLT_HT500_v8','HLT_HT500_v9','HLT_HT500_v10','HLT_HT500_v11',
      'HLT_HT550_v1','HLT_HT550_v2','HLT_HT550_v3','HLT_HT550_v4','HLT_HT550_v5','HLT_HT550_v6','HLT_HT550_v7','HLT_HT550_v8','HLT_HT550_v9','HLT_HT550_v10','HLT_HT550_v11',
      'HLT_HT600_v1','HLT_HT600_v2','HLT_HT600_v3','HLT_HT600_v4',
      'HLT_HT650_v1','HLT_HT650_v2','HLT_HT650_v3','HLT_HT650_v4',
      'HLT_HT700_v1','HLT_HT700_v2',
      'HLT_QuadJet40_v1','HLT_QuadJet40_v2','HLT_QuadJet40_v3','HLT_QuadJet40_v4','HLT_QuadJet40_v5','HLT_QuadJet40_v6','HLT_QuadJet40_v7','HLT_QuadJet40_v8','HLT_QuadJet40_v9','HLT_QuadJet40_v10','HLT_QuadJet40_v11',
      'HLT_QuadJet70_v1','HLT_QuadJet70_v2','HLT_QuadJet70_v3','HLT_QuadJet70_v4','HLT_QuadJet70_v5','HLT_QuadJet70_v6','HLT_QuadJet70_v7','HLT_QuadJet70_v8','HLT_QuadJet70_v9','HLT_QuadJet70_v10',
      'HLT_QuadJet80_v1','HLT_QuadJet80_v2','HLT_QuadJet80_v3','HLT_QuadJet80_v4','HLT_QuadJet80_v5',
      'HLT_EightJet35_v1','HLT_EightJet35_v2','HLT_EightJet35_v3'
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
    HLTPaths           = cms.vstring('HLT_HT450_v*','HLT_HT500_v*','HLT_HT550_v*','HLT_HT600_v*','HLT_HT650_v*','HLT_HT700_v*',
                                     'HLT_QuadJet40_v*','HLT_QuadJet70_v*','HLT_QuadJet80_v*', 
                                     'HLT_EightJet35_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.path = cms.Path(process.hltFilter * process.kt6PFJets * process.kt6CaloJets * process.ak5PFJets * process.ak5)


