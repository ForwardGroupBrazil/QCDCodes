import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_42_V12::All'
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
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Summer11/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6/AODSIM/PU_S3_START42_V11-v2/0004/FA6FEF7F-8E7E-E011-AC35-001A92971B78.root')
)
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_mc.root'))

process.ak7 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.InputTag('ak7PFJets'),
    calojets        = cms.InputTag('ak7CaloJets'),
    genjets         = cms.untracked.InputTag('ak7GenJets'),
    isMCarlo        = cms.untracked.bool(True),
    ## database entry for the uncertainties ######
    PFPayloadName   = cms.string(''),
    CaloPayloadName = cms.string(''),
    ## calojet ID and extender for the JTA #######
    calojetID       = cms.InputTag('ak7JetID'),
    calojetExtender = cms.InputTag('ak7JetExtender'),
    ## set the conditions for good Vtx counting ##
    offlineVertices = cms.InputTag('offlinePrimaryVertices'),
    goodVtxNdof     = cms.double(4), 
    goodVtxZ        = cms.double(24),
    ## number of jets to be stored ###############
    maxY            = cms.double(5.0), 
    minPFPt         = cms.double(20),
    minPFFatPt      = cms.double(10),
    minCaloPt       = cms.double(20),
    minGenPt        = cms.untracked.double(20),
    minNPFJets      = cms.int32(1),
    minNCaloJets    = cms.int32(1), 
    minJJMass       = cms.double(-1),
    ## trigger ##############################
    printTriggerMenu = cms.untracked.bool(True),
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring('HLT_Jet110'),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jec services ##############################
    pfjecService    = cms.string('ak7PFL1FastL2L3'),
    calojecService  = cms.string('ak7CaloL1L2L3')
)

process.ak5 = process.ak7.clone(
    pfjets           = 'ak5PFJets',
    calojets         = 'ak5CaloJets',
    genjets          = 'ak7GenJets',
    PFPayloadName    = '',
    CaloPayloadName  = '',
    calojetID        = 'ak5JetID',
    calojetExtender  = 'ak5JetExtender',
    pfjecService     = 'ak5PFL1FastL2L3',
    calojecService   = 'ak5CaloL1L2L3',
    printTriggerMenu = False 
)

process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(5.0)
process.ak7PFJets.doAreaFastjet = True
process.ak7PFJets.Rho_EtaMax = cms.double(5.0)
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(5.0)

process.path = cms.Path(process.kt6PFJets * process.ak5PFJets * process.ak7PFJets * process.ak5 * process.ak7)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


