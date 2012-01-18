from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

##--------- global tag -------------------------
process.GlobalTag.globaltag = 'GR_R_42_V23::All'

##--------- remove cleaning --------------------
removeCleaning(process)
##--------- remove MC matching -----------------
#removeMCMatching(process)
##--------- add PFMet --------------------------
#addPfMET(process, 'PF')
##--------- jets -------------------------------
process.patJets.embedPFCandidates = False
process.patJets.embedCaloTowers = False
process.patJets.addTagInfos = True
process.patJets.userData.userFloats.src = ['betaAK5PF','betaAK5PFchs']
##--------- good primary vertices ---------------
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    src          = cms.InputTag('offlinePrimaryVertices'),
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) )
)
##--------- value map for betastar variable -----
process.betaAK5PF = cms.EDProducer("JetBetaProducer",
          jets   = cms.InputTag("ak5PFJets")
)
process.betaAK5PFchs = cms.EDProducer("JetBetaProducer",
          jets   = cms.InputTag("pfJetsCHS")
)
##--------- PF2PAT -----------------------------
postfix = 'CHS'
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=False, postfix=postfix,
          jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']))

removeMCMatchingPF2PAT(process,'')

process.pfPileUpCHS.Enable = True
process.pfPileUpCHS.checkClosestZVertex = cms.bool(False)
process.pfPileUpCHS.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfJetsCHS.doAreaFastjet = True
process.pfJetsCHS.doRhoFastjet = False

process.ak5PFJets.doAreaFastjet = True
process.kt6CaloJets.doRhoFastjet = True

process.kt6PFJetsCHS = process.kt6PFJets.clone(
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True)
    )
process.patJetCorrFactorsCHS.rho = cms.InputTag("kt6PFJetsCHS", "rho")

getattr(process,"patPF2PATSequence"+postfix).replace(
    getattr(process,"pfNoElectron"+postfix),
    getattr(process,"pfNoElectron"+postfix)*process.kt6PFJetsCHS)

getattr(process,"patPF2PATSequence"+postfix).replace(
    getattr(process,"pfJets"+postfix),
    getattr(process,"pfJets"+postfix)*process.betaAK5PFchs)

process.patPF2PATseq = cms.Sequence(
    getattr(process,"patPF2PATSequence"+postfix)
    )

##--------- remove MC matching -----------------
removeMCMatching(process)
addPfMET(process, 'PF')
switchJetCollection(process,cms.InputTag('ak5PFJets'),
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'])),
                 doType1MET   = False,
                 doJetID      = False
                 )

addJetCollection(process,cms.InputTag('ak5CaloJets'),
                 'AK5', 'Calo',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5Calo', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute','L2L3Residual'])),
                 doType1MET   = False,
                 jetIdLabel   = "ak5",
                 doJetID      = True
                 )

process.selectedPatJets.cut        = "pt > 20 && abs(eta) < 3.0"
process.selectedPatJetsAK5Calo.cut = "pt > 20 && abs(eta) < 3.0"
process.selectedPatJetsCHS.cut     = "pt > 20 && abs(eta) < 3.0"
##--------- keep only jet and MET PAT objects ---
removeAllPATObjectsBut(process,["Jets","METs"])
##--------- output commands ---------------------
process.out.fileName = 'patuple_multijets.root'
process.out.outputCommands = [
         'keep *_kt6PFJets_rho_PAT',
         'keep *_kt6PFJetsCHS_rho_PAT',
         'keep *_kt6CaloJets_rho_PAT', 
         'keep *_selectedPatJets*__*',
         'keep *_met_*_*',
         'keep *_pfMet_*_*', 
         'keep recoVertexs_*fflinePrimaryVertices_*_*',
         'keep edmTriggerResults_TriggerResults_*_*',
         'keep *_hltTriggerSummaryAOD_*_*',
         'keep L1GlobalTriggerObjectMapRecord_*_*_*',
         'keep L1GlobalTriggerReadoutRecord_*_*_*',
         'keep *_beta*_*_*'
]

process.multiJetFilter = cms.EDFilter('PatMultijetFilter',
    jets     = cms.InputTag('selectedPatJets'),
    minNjets = cms.int32(4),
    minPt    = cms.double(30)
)

process.maxEvents.input = 5000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source.fileNames = [
'/store/data/Run2011B/HT/AOD/PromptReco-v1/000/180/252/6AEEE163-1A05-E111-9DD8-BCAEC518FF50.root'
]

process.options.wantSummary = False

process.p = cms.Path(
   process.ak5PFJets +
   process.kt6PFJets +
   process.kt6CaloJets +
   process.goodOfflinePrimaryVertices +
   process.betaAK5PF +
   process.patPF2PATseq +
   process.patDefaultSequence +
   process.multiJetFilter
)

