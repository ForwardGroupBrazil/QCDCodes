from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.load('CMGTools.External.pujetidsequence_cff')

from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

##--------- global tag -------------------------
process.GlobalTag.globaltag = 'GR_R_52_V7::All'

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
          jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']))

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

#----- recommendation from JES: use the standard rho for CHS ------
process.patJetCorrFactorsCHS.rho = cms.InputTag("kt6PFJets", "rho")

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
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'])),
                 doType1MET   = False,
                 doJetID      = False
                 )

process.selectedPatJets.cut        = "pt > 10 && abs(eta) < 4.7"
process.selectedPatJetsCHS.cut     = "pt > 10 && abs(eta) < 4.7"

##---- modify PU jet id -------------
process.puJetId.jets = cms.InputTag('selectedPatJets')
process.puJetId.vertexes = cms.InputTag('goodOfflinePrimaryVertices')
process.puJetMva.jets = cms.InputTag('selectedPatJets')
process.puJetMva.vertexes = cms.InputTag('goodOfflinePrimaryVertices')
process.puJetIdChs.jets = cms.InputTag('selectedPatJetsCHS')
process.puJetIdChs.vertexes = cms.InputTag('goodOfflinePrimaryVertices')
process.puJetMvaChs.jets = cms.InputTag('selectedPatJetsCHS')
process.puJetMvaChs.vertexes = cms.InputTag('goodOfflinePrimaryVertices')

##--------- keep only jet and MET PAT objects ---
removeAllPATObjectsBut(process,["Jets","METs"])
##--------- output commands ---------------------
process.out.fileName = 'patTuple.root'
process.out.outputCommands = [
         'keep *_kt6PFJets_rho_PAT',
         'keep *_kt6PFJetsCHS_rho_PAT',
         'keep *_kt6PFJetsISO_rho_PAT',## needed for the QG likelihood
         'keep *_selectedPatJets__*',
         'keep *_selectedPatJetsCHS__*',
         'keep *_puJet*_*_*',
         'keep *_jetExtender*_*_*', 
         'keep *_ak5SoftTrackJets__*',
         'keep *_HBHENoiseFilterResultProducer_*_*', 
         'keep *_pfMet_*_*', 
         'keep recoVertexs_goodOfflinePrimaryVertices_*_*',
         'keep edmTriggerResults_TriggerResults_*_HLT',
         'keep *_hltTriggerSummaryAOD_*_*'
         #'keep L1GlobalTriggerObjectMapRecord_*_*_*',
         #'keep L1GlobalTriggerReadoutRecord_*_*_*',
]

process.outTracks = cms.EDProducer('PatTracksOutOfJets',
    jets    = cms.InputTag('selectedPatJets'),
    vtx     = cms.InputTag('goodOfflinePrimaryVertices'),
    tracks  = cms.InputTag('generalTracks'),
    btagger = cms.string('combinedSecondaryVertexBJetTags')
)

process.ak5SoftTrackJets = process.ak5TrackJets.clone(src = 'outTracks',jetPtMin = 1.0)

process.jetExtender = cms.EDProducer("JetExtendedProducer",
    jets    = cms.InputTag('selectedPatJets'),
    result  = cms.string('extendedPatJets'),
    payload = cms.string('AK5PF')
)

process.jetExtenderCHS = cms.EDProducer("JetExtendedProducer",
    jets    = cms.InputTag('selectedPatJetsCHS'),
    result  = cms.string('extendedPatJetsCHS'),
    payload = cms.string('AK5PFchs') 
)

############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_QuadJet75_55_35_20_BTagIP_VBF_v*','HLT_QuadJet75_55_38_20_BTagIP_VBF_v*',
                                     'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v*','HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v*'
                                     ),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.maxEvents.input = 100
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source.fileNames = [
'/store/data/Run2012A/MultiJet/AOD/PromptReco-v1/000/190/467/D8C5020C-F980-E111-84F0-003048F118C2.root',
'/store/data/Run2012A/MultiJet/AOD/PromptReco-v1/000/190/519/E65A5017-6B81-E111-9CF6-001D09F291D2.root',
'/store/data/Run2012A/MultiJet/AOD/PromptReco-v1/000/190/517/0820E829-6681-E111-AFD6-003048F110BE.root',
'/store/data/Run2012A/MultiJet/AOD/PromptReco-v1/000/190/492/EEA1A177-2981-E111-AB04-003048F117B6.root',
'/store/data/Run2012A/MultiJet/AOD/PromptReco-v1/000/190/491/BE07ECB3-1C81-E111-96FA-BCAEC5329713.root',
'/store/data/Run2012A/MultiJet/AOD/PromptReco-v1/000/190/490/E8ACC615-2F81-E111-859C-001D09F23D1D.root'
]

process.options.wantSummary = False

process.p = cms.Path(
   #----- produce the HBHE noise flag --------------------------
   process.HBHENoiseFilterResultProducer +
   #----- re-cluster kt6PFJets after activating rho for ISO ----
   process.kt6PFJetsISO +
   #----- create the collection of good PV ---------------------
   process.goodOfflinePrimaryVertices +
   #----- run the PF2PAT sequence: doing CHS -------------------
   process.patPF2PATseq +
   #----- run the default PAT sequence -------------------------
   process.patDefaultSequence +
   #----- pu jet id --------------------------------------------
   process.puJetIdSqeuence +
   process.puJetIdSqeuenceChs +
   #----- extend the PAT jets with additional variables --------
   process.jetExtender +
   #----- extend the CHS PAT jets with additional variables ----
   process.jetExtenderCHS +
   #----- create a collection of tracks out of jets ------------
   process.outTracks +
   #----- reconstruct track jets from the soft tracks ----------
   process.ak5SoftTrackJets
)

