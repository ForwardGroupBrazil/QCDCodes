import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

###############################
####### Parameters ############
###############################

options = VarParsing ('python')

options.register('runOnMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on MC"
)
options.register('reportEvery', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1)"
)
options.register('wantSummary', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
options.register('mcGlobalTag', 'START53_V27',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag"
)
options.register('dataGlobalTag', 'FT53_V21A_AN6',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Data global tag"
)

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()

print "Running on MC: %s"%('True' if options.runOnMC else 'False')

## Global tag
globalTag = options.mcGlobalTag
if not options.runOnMC:
    globalTag = options.dataGlobalTag

## Jet energy corrections
jetCorrectionsAK5PFchs = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
jetCorrectionsAK7PFchs = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
jetCorrectionsAK5Calo = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
jetCorrectionsAK7Calo = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if not options.runOnMC:
    jetCorrectionsAK5PFchs[1].append('L2L3Residual')
    jetCorrectionsAK7PFchs[1].append('L2L3Residual')
    jetCorrectionsAK5Calo[1].append('L2L3Residual')
    jetCorrectionsAK7Calo[1].append('L2L3Residual')

triggerNames = cms.vstring('HLT_PFJet40_v2')
if not options.runOnMC:
    triggerNames = cms.vstring('HLT_PFJet40_v3','HLT_PFJet40_v4','HLT_PFJet40_v5',
                               'HLT_PFJet80_v3','HLT_PFJet80_v4','HLT_PFJet80_v5',
                               'HLT_PFJet140_v3','HLT_PFJet140_v4','HLT_PFJet140_v5',
                               'HLT_PFJet200_v3','HLT_PFJet200_v4','HLT_PFJet200_v5',
                               'HLT_PFJet260_v3','HLT_PFJet260_v4','HLT_PFJet260_v5',
                               'HLT_PFJet320_v3','HLT_PFJet320_v4','HLT_PFJet320_v5',
                               'HLT_PFJet400_v3','HLT_PFJet400_v4','HLT_PFJet400_v5',
                               'HLT_HT200_v1','HLT_HT200_v2','HLT_HT200_v3','HLT_HT200_v4','HLT_HT200_v5',
                               'HLT_HT250_v1','HLT_HT250_v2','HLT_HT250_v3','HLT_HT250_v4','HLT_HT250_v5',
                               'HLT_HT300_v1','HLT_HT300_v2','HLT_HT300_v3','HLT_HT300_v4','HLT_HT300_v5',
                               'HLT_HT350_v1','HLT_HT350_v2','HLT_HT350_v3','HLT_HT350_v4','HLT_HT350_v5',
                               'HLT_HT400_v1','HLT_HT400_v2','HLT_HT400_v3','HLT_HT400_v4','HLT_HT400_v5',
                               'HLT_HT450_v1','HLT_HT450_v2','HLT_HT450_v3','HLT_HT450_v4','HLT_HT450_v5',
                               'HLT_HT500_v1','HLT_HT500_v2','HLT_HT500_v3','HLT_HT500_v4','HLT_HT500_v5','HLT_HT500_v6','HLT_HT500_v7','HLT_HT500_v8','HLT_HT500_v9','HLT_HT500_v10',
                               'HLT_HT550_v1','HLT_HT550_v2','HLT_HT550_v3','HLT_HT550_v4','HLT_HT550_v5','HLT_HT550_v6','HLT_HT550_v7','HLT_HT550_v8','HLT_HT550_v9','HLT_HT550_v10',
                               'HLT_HT650_v1','HLT_HT650_v2','HLT_HT650_v3','HLT_HT650_v4','HLT_HT650_v5','HLT_HT650_v6','HLT_HT650_v7','HLT_HT650_v8','HLT_HT650_v9','HLT_HT650_v10',
                               'HLT_HT750_v1','HLT_HT750_v2','HLT_HT750_v3','HLT_HT750_v4','HLT_HT750_v5','HLT_HT750_v6','HLT_HT750_v7','HLT_HT750_v8','HLT_HT750_v9','HLT_HT750_v10',
                               'HLT_PFHT650_v1','HLT_PFHT650_v2','HLT_PFHT650_v3','HLT_PFHT650_v4','HLT_PFHT650_v5','HLT_PFHT650_v6','HLT_PFHT650_v7','HLT_PFHT650_v8','HLT_PFHT650_v9','HLT_PFHT650_v10',
                               'HLT_PFHT700_v1','HLT_PFHT700_v2','HLT_PFHT700_v3','HLT_PFHT700_v4','HLT_PFHT700_v5','HLT_PFHT700_v6','HLT_PFHT700_v7','HLT_PFHT700_v8','HLT_PFHT700_v9','HLT_PFHT700_v10',
                               'HLT_PFHT750_v1','HLT_PFHT750_v2','HLT_PFHT750_v3','HLT_PFHT750_v4','HLT_PFHT750_v5','HLT_PFHT750_v6','HLT_PFHT750_v7','HLT_PFHT750_v8','HLT_PFHT750_v9','HLT_PFHT750_v10',
                               'HLT_PFNoPUHT650_v1','HLT_PFNoPUHT650_v2','HLT_PFNoPUHT650_v3','HLT_PFNoPUHT650_v4','HLT_PFNoPUHT650_v5','HLT_PFNoPUHT650_v6','HLT_PFNoPUHT650_v7','HLT_PFNoPUHT650_v8','HLT_PFNoPUHT650_v9','HLT_PFNoPUHT650_v10',
                               'HLT_PFNoPUHT700_v1','HLT_PFNoPUHT700_v2','HLT_PFNoPUHT700_v3','HLT_PFNoPUHT700_v4','HLT_PFNoPUHT700_v5','HLT_PFNoPUHT700_v6','HLT_PFNoPUHT700_v7','HLT_PFNoPUHT700_v8','HLT_PFNoPUHT700_v9','HLT_PFNoPUHT700_v10',
                               'HLT_PFNoPUHT750_v1','HLT_PFNoPUHT750_v2','HLT_PFNoPUHT750_v3','HLT_PFNoPUHT750_v4','HLT_PFNoPUHT750_v5','HLT_PFNoPUHT750_v6','HLT_PFNoPUHT750_v7','HLT_PFNoPUHT750_v8','HLT_PFNoPUHT750_v9','HLT_PFNoPUHT750_v10',
                               'HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v1','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v6','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v7','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v8','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v9','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10')

## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## MessageLogger
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10

## HBHENoiseFilterResultProducer
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *

# An empty postfix means that only PF2PAT is run,
# otherwise both standard PAT and PF2PAT are run. In the latter case PF2PAT
# collections have standard names + postfix (e.g. patElectronPFlow)
postfix = "PFlow"
jetAlgo = "AK5"
usePF2PAT(process, runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=options.runOnMC, postfix=postfix,
          jetCorrections=jetCorrectionsAK5PFchs, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

# to run second PF2PAT+PAT with different postfix uncomment the following lines
# and add the corresponding sequence to path
postfix2 = "PFlow2"
jetAlgo2="AK7"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo2, runOnMC=options.runOnMC, postfix=postfix2,
          jetCorrections=jetCorrectionsAK7PFchs, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

# top projections in PF2PAT:
getattr(process,"pfPileUp"+postfix).checkClosestZVertex = False
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

getattr(process,"pfPileUp"+postfix2).checkClosestZVertex = False
getattr(process,"pfNoPileUp"+postfix2).enable = True
getattr(process,"pfNoMuon"+postfix2).enable = True
getattr(process,"pfNoElectron"+postfix2).enable = True
getattr(process,"pfNoTau"+postfix2).enable = False
getattr(process,"pfNoJet"+postfix2).enable = True

# verbose flags for the PF2PAT modules
getattr(process,"pfNoMuon"+postfix).verbose = False
getattr(process,"pfNoMuon"+postfix2).verbose = False

if not options.runOnMC:
    # removing MC matching for standard PAT sequence
    # for the PF2PAT+PAT sequence, it is done in the usePF2PAT function
    removeMCMatching( process, ['All'] )

## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(process,
    cms.InputTag('ak5CaloJets'),
    doJTA = True,
    doBTagging = True,
    jetCorrLabel = jetCorrectionsAK5Calo,
    doType1MET = False,
    genJetCollection = cms.InputTag('ak5GenJetsNoNu'),
    doJetID = True
)
addJetCollection(
    process,
    cms.InputTag('ak7CaloJets'),
    'AK7','Calo',
    doJTA=True,
    doBTagging=True,
    jetCorrLabel=jetCorrectionsAK7Calo,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=True,
    jetIdLabel='ak7',
    genJetCollection=cms.InputTag('ak7GenJetsNoNu')
)

## Keep only jets in the PAT sequence
removeSpecificPATObjects(process,names=['Photons', 'Electrons', 'Muons', 'Taus', 'METs'],postfix='')
removeSpecificPATObjects(process,names=['Photons', 'Electrons', 'Muons', 'Taus', 'METs'],postfix=postfix)
removeSpecificPATObjects(process,names=['Photons', 'Electrons', 'Muons', 'Taus', 'METs'],postfix=postfix2)

## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
process.GlobalTag.globaltag = globalTag + '::All'    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#                                         ##
process.source.fileNames =  ['/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v1/00000/FE7C71D8-DB25-E211-A93B-0025901D4C74.root']
if not options.runOnMC:
    process.source.fileNames =  ['/store/data/Run2012A/Jet/AOD/22Jan2013-v1/20000/30B21345-4172-E211-9EF3-00304867BEC0.root']
#                                         ##
process.maxEvents.input = options.maxEvents
#                                         ##
#   process.out.outputCommands = [ ... ]  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#                                         ##
process.out.fileName = 'patTuple_PATandPF2PAT.root'
#                                         ##
process.options.wantSummary = options.wantSummary


process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_' + ('mc' if options.runOnMC else 'data') + '.root'))

process.ak5 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.InputTag('selectedPatJetsPFlow'),
    calojets        = cms.InputTag('selectedPatJets'),
    genjets         = cms.untracked.InputTag('ak5GenJetsNoNu'),
    ## MET collections ###########################
    pfmet           = cms.InputTag('pfMETPFlow'),
    calomet         = cms.InputTag('met'),
    ## database entry for the uncertainties ######
    PFPayloadName   = cms.string('AK5PFchs'),
    CaloPayloadName = cms.string('AK5Calo'),
    jecUncSrc       = cms.string('Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt'),
    jecUncSrcNames  = cms.vstring('Absolute','HighPtExtra','SinglePionECAL','SinglePionHCAL','FlavorQCD','Time',
                                  'RelativeJEREC1','RelativeJEREC2','RelativeJERHF','RelativePtBB','RelativePtEC1','RelativePtEC2','RelativePtHF','RelativeFSR',
                                  'RelativeStatEC2','RelativeStatHF','PileUpDataMC','PileUpPtBB','PileUpPtEC','PileUpPtHF','PileUpBias',
                                  'SubTotalPileUp','SubTotalRelative','SubTotalPt','SubTotalMC','Total','TotalNoFlavor',
                                  'FlavorZJet','FlavorPhotonJet','FlavorPureGluon','FlavorPureQuark','FlavorPureCharm','FlavorPureBottom'),
    ## calojet extender for the JTA #######
    calojetExtender = cms.InputTag('ak5JetExtender'),
    ## set the conditions for good Vtx counting ##
    offlineVertices = cms.InputTag('goodOfflinePrimaryVertices'),
    goodVtxNdof     = cms.double(4),
    goodVtxZ        = cms.double(24),
    ## rho #######################################
    srcCaloRho      = cms.InputTag('kt6CaloJets','rho'),
    srcPFRho        = cms.InputTag('kt6PFJets','rho'),
    ## MC Generator flags ######################
    useGenInfo      = cms.untracked.bool(False),
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
    triggerName     = triggerNames,
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## Cross section ############################
    Xsec            = cms.untracked.double(1)
)

process.ak7 = process.ak5.clone(
    pfjets           = 'selectedPatJetsPFlow2',
    calojets         = 'selectedPatJetsAK7Calo',
    PFPayloadName    = 'AK7PFchs',
    CaloPayloadName  = 'AK7Calo',
    jecUncSrc        = 'Summer13_V4_DATA_UncertaintySources_AK7PFchs.txt',
    calojetExtender  = 'ak7JetExtender',
    printTriggerMenu = False
)

############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring(
        'HLT_PFJet40_v*','HLT_PFJet80_v*','HLT_PFJet140_v*','HLT_PFJet200_v*','HLT_PFJet260_v*','HLT_PFJet320_v*','HLT_PFJet400_v*','HLT_HT200_v*','HLT_HT250_v*','HLT_HT300_v*','HLT_HT350_v*','HLT_HT400_v*','HLT_HT450_v*','HLT_HT500_v*','HLT_HT550_v*','HLT_HT650_v*','HLT_HT750_v*','HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v*','HLT_PFHT650_v*','HLT_PFHT700_v*','HLT_PFHT750_v*','HLT_PFNoPUHT650_v*','HLT_PFNoPUHT700_v*','HLT_PFNoPUHT750_v*'
    ),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

#-------------------------------------
## Produce a collection of good primary vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
)

## Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

## Filter for good primary vertex
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)
#-------------------------------------

#---------------------------------------
## Optional MET filters:
## https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
process.load("RecoMET.METFilters.metFilters_cff")
process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
#---------------------------------------

#---------------------------------------
## Filter for HCAL laser events in prompt 2012A+B+C, snippet for "Datasets from the 2013 rereco and Multijet parked":
## https://twiki.cern.ch/twiki/bin/view/CMS/PdmVKnowFeatures#HCAL_laser_events_in_prompt_2012
process.load("EventFilter.HcalRawToDigi.hcallaserFilterFromTriggerResult_cff")
#---------------------------------------

## Define event filter sequence
process.filtSeq = cms.Sequence(
    process.noscraping
    * process.primaryVertexFilter
    * process.goodOfflinePrimaryVertices
    * process.HBHENoiseFilter
    * process.CSCTightHaloFilter
    * process.EcalDeadCellTriggerPrimitiveFilter
    * process.eeBadScFilter
    * process.ecalLaserCorrFilter
    * process.trackingFailureFilter
    * process.trkPOGFilters
)
process.filtSeqData = cms.Sequence(
    process.hltFilter
    * process.hcalfilter
    * process.HBHENoiseFilterResultProducer
)
if not options.runOnMC:
    process.filtSeq *= process.filtSeqData

# Let it run
process.p = cms.Path(
    process.filtSeq
    * ( getattr(process,"patPF2PATSequence"+postfix)
    + getattr(process,"patPF2PATSequence"+postfix2)
    + process.patDefaultSequence
    )
    * ( process.ak5
    + process.ak7
    )
)

# Delete predefined output module (needed for running with CRAB)
del process.out
del process.outpath
