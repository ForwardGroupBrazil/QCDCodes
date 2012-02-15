from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

##--------- global tag -------------------------
process.GlobalTag.globaltag = 'GR_R_42_V23::All'

##--------- remove cleaning --------------------
removeCleaning(process)
##--------- jets -------------------------------
process.patJets.embedPFCandidates = False
process.patJets.embedCaloTowers = False
process.patJets.addTagInfos = True
##--------- good primary vertices ---------------
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    src          = cms.InputTag('offlinePrimaryVertices'),
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) )
)

##--------- PF2PAT -----------------------------
postfix = 'CHS'
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=False, postfix=postfix,
          jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute']))

removeMCMatchingPF2PAT(process,'')

process.pfPileUpCHS.Enable = True
process.pfPileUpCHS.checkClosestZVertex = cms.bool(False)
process.pfPileUpCHS.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfJetsCHS.doAreaFastjet = True
process.pfJetsCHS.doRhoFastjet = False

process.ak5PFJets.doAreaFastjet = True
process.kt6PFJets.doRhoFastjet = True

process.kt6PFJetsCHS = process.kt6PFJets.clone(
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True)
    )

process.kt6PFJetsISO = process.kt6PFJets.clone(
    Rho_EtaMax = cms.double(2.4)
    )

process.patJetCorrFactorsCHS.rho = cms.InputTag("kt6PFJetsCHS", "rho")

getattr(process,"patPF2PATSequence"+postfix).replace(
    getattr(process,"pfNoElectron"+postfix),
    getattr(process,"pfNoElectron"+postfix)*process.kt6PFJetsCHS)

process.patPF2PATseq = cms.Sequence(
    getattr(process,"patPF2PATSequence"+postfix)
    )

##--------- remove MC matching -----------------
removeMCMatching(process)
addPfMET(process, 'PF')
switchJetCollection(process,cms.InputTag('ak5PFJets'),
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = False,
                 doJetID      = False
                 )

process.selectedPatJets.cut        = "pt > 10 && abs(eta) < 4.7"
process.selectedPatJetsCHS.cut     = "pt > 10 && abs(eta) < 4.7"

##--------- keep only jet and MET PAT objects ---
removeAllPATObjectsBut(process,["Jets","METs"])
##--------- output commands ---------------------
process.out.fileName = 'test_vbf_patuple_mc.root'
process.out.outputCommands = [
         'keep *_kt6PFJets_rho_PAT',
         'keep *_kt6PFJetsCHS_rho_PAT',
         'keep *_kt6PFJetsISO_rho_PAT',## needed for the QG likelihood
         #'keep *_selectedPatJets__*',
         #'keep *_selectedPatJetsCHS__*',
         'keep *_jetExtender*_*_*', 
         'keep *_addPileupInfo_*_*',  
         'keep recoGenJets_ak5GenJets_*_*',
         'keep *_genParticles_*_*',
         'keep *_HBHENoiseFilterResultProducer_*_*', 
         'keep *_pfMet_*_*', 
         'keep recoVertexs_goodOfflinePrimaryVertices_*_*',
         'keep edmTriggerResults_TriggerResults_*_HLT',
         'keep *_hltTriggerSummaryAOD_*_*'
         #'keep L1GlobalTriggerObjectMapRecord_*_*_*',
         #'keep L1GlobalTriggerReadoutRecord_*_*_*',
]

process.jetExtender = cms.EDProducer("JetExtendedProducer",
    jets    = cms.InputTag('selectedPatJets'),
    payload = cms.string('AK5PF'),
    result  = cms.string('extendedPatJets') 
)

process.jetExtenderCHS = cms.EDProducer("JetExtendedProducer",
    jets   = cms.InputTag('selectedPatJetsCHS'),
    payload = cms.string('AK5PFchs'),
    result = cms.string('extendedPatJetsCHS') 
)

process.multiJetFilter = cms.EDFilter('PatMultijetFilter',
    jets     = cms.InputTag('selectedPatJets'),
    minNjets = cms.int32(4),
    minPt    = cms.double(20)
)

process.maxEvents.input = 100
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source.fileNames = [
'/store/mc/Summer11/QCD_TuneZ2_HT-1000_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/1AB5A492-C4C5-E011-BCD9-90E6BA19A203.root'
]

process.options.wantSummary = False

process.p = cms.Path(
   process.HBHENoiseFilterResultProducer +
   process.ak5PFJets +
   process.kt6PFJets +
   process.kt6PFJetsISO +
   process.goodOfflinePrimaryVertices +
   process.patPF2PATseq +
   process.patDefaultSequence +
   process.jetExtender +
   process.jetExtenderCHS +
   process.multiJetFilter
)

