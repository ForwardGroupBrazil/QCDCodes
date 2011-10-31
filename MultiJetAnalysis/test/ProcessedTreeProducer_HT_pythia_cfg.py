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
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(4300)
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#############   Define the source file ###############
process.load('KKousour.QCDAnalysis.axiaxi800_hulubar_ussbar_hulbaru_ubarddbar_cfi')
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('/uscms_data/d2/kkousour/ProcessedTree_axiaxi800.root'))

process.ak5 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.InputTag('ak5PFJets'),
    calojets        = cms.InputTag('ak5CaloJets'),
    genjets         = cms.untracked.InputTag('ak5GenJets'),
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
    ## MC $ Generator flags ######################
    isMCarlo        = cms.untracked.bool(True),
    useGenInfo      = cms.untracked.bool(False),
    ## simulated PU ##############################
    srcPU           = cms.untracked.InputTag('addPileupInfo'),
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
    triggerName     = cms.vstring(),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jec services ##############################
    pfjecService    = cms.string('ak5PFL1FastL2L3'),
    calojecService  = cms.string('ak5CaloL1L2L3')
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

process.path = cms.Path(process.kt6PFJets * process.kt6CaloJets * process.ak5PFJets * process.ak5)


