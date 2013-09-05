import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START53_V27::All'
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##---- Jet-Flavor Matching ------------------------------------------
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi") 
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 50

#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
 '/store/mc/Summer12_DR53X/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v2/00000/D8C12B63-760E-E211-AC09-002618943867.root'
)
)
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_mc.root'))

process.ak7 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.InputTag('ak7PFJets'),
    calojets        = cms.InputTag('ak7CaloJets'),
    genjets         = cms.untracked.InputTag('ak7GenJets'),
    ## database entry for the uncertainties ######
    PFPayloadName   = cms.string('AK7PF'),
    CaloPayloadName = cms.string('AK7Calo'),   
    jecUncSrc       = cms.string('Summer12_V2_DATA_AK7PF_UncertaintySources.txt'),
    jecUncSrcNames  = cms.vstring('Absolute','HighPtExtra','SinglePion','Flavor','Time',
                                  'RelativeJEREC1','RelativeJEREC2','RelativeJERHF',
                                  'RelativeStatEC2','RelativeStatHF','RelativeFSR',
                                  'PileUpDataMC','PileUpOOT','PileUpPt','PileUpBias','PileUpJetRate',
                                  'SubTotalPileUp','SubTotalRelative','SubTotalPt','SubTotalDataMC','Total'),
    ## calojet ID and extender for the JTA #######
    calojetID       = cms.InputTag('ak7JetID'),
    calojetExtender = cms.InputTag('ak7JetExtender'),
    ## set the conditions for good Vtx counting ##
    offlineVertices = cms.InputTag('offlinePrimaryVertices'),
    goodVtxNdof     = cms.double(4),
    goodVtxZ        = cms.double(24),
    ## rho #######################################
    srcCaloRho      = cms.InputTag('kt6CaloJets','rho'),
    srcPFRho        = cms.InputTag('kt6PFJets','rho'),
    ## MC $ Generator flags ######################
    isMCarlo        = cms.untracked.bool(True),
    useGenInfo      = cms.untracked.bool(True),
    ## simulated PU ##############################
    srcPU           = cms.untracked.InputTag('addPileupInfo'),
    ## preselection cuts #########################
    maxY            = cms.double(5.0),
    minPFPt         = cms.double(20),
    minPFFatPt      = cms.double(30),
    maxPFFatEta     = cms.double(2.5),
    minCaloPt       = cms.double(20),
    minGenPt        = cms.untracked.double(20),
    minNPFJets      = cms.int32(1),
    minNCaloJets    = cms.int32(1),
    minJJMass       = cms.double(-1),
    ## trigger ###################################
    printTriggerMenu = cms.untracked.bool(True),
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring('HLT_PFJet40_v2'),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jec services ##############################
    pfjecService    = cms.string('ak7PFL1FastL2L3'),
    calojecService  = cms.string('ak7CaloL1FastL2L3'),
    ## Jet-Flavor Matching #######################
    jetFlavourMatching     = cms.untracked.string('AK7byValAlgo'),
    Xsec           = cms.untracked.double(1033680.0)
)

process.ak5 = process.ak7.clone(
    pfjets           = 'ak5PFJets',
    calojets         = 'ak5CaloJets',
    genjets          = 'ak5GenJets',
    PFPayloadName    = 'AK5PF',
    CaloPayloadName  = 'AK5Calo',
    jecUncSrc        = 'Summer12_V2_DATA_AK5PF_UncertaintySources.txt',
    calojetID        = 'ak5JetID',
    calojetExtender  = 'ak5JetExtender',
    pfjecService     = 'ak5PFL1FastL2L3',
    calojecService   = 'ak5CaloL1FastL2L3',
    jetFlavourMatching = 'AK5byValAlgo',
    printTriggerMenu = False
)

process.AK7PFbyRef = process.AK7byRef.clone(jets = cms.InputTag("ak7PFJets"))
process.AK7PFbyValPhys = process.AK7byValPhys.clone(srcByReference = cms.InputTag("AK7PFbyRef"))
process.AK7PFbyValAlgo = process.AK7byValAlgo.clone(srcByReference = cms.InputTag("AK7PFbyRef"))
process.AK7PFFlavour = cms.Sequence(process.AK7PFbyRef*process.AK7PFbyValPhys*process.AK7PFbyValAlgo)

process.AK5PFbyRef = process.AK5byRef.clone(jets = cms.InputTag("ak5PFJets"))
process.AK5PFbyValPhys = process.AK5byValPhys.clone(srcByReference = cms.InputTag("AK5PFbyRef"))
process.AK5PFbyValAlgo = process.AK5byValAlgo.clone(srcByReference = cms.InputTag("AK5PFbyRef"))
process.AK5PFFlavour = cms.Sequence(process.AK5PFbyRef*process.AK5PFbyValPhys*process.AK5PFbyValAlgo)

process.path = cms.Path(process.myPartons * (process.AK7Flavour + process.AK7PFFlavour) * (process.AK5Flavour + process.AK5PFFlavour) * process.ak7 * process.ak5)


