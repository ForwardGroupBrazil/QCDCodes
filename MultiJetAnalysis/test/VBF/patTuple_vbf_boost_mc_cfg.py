from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')


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
process.out.fileName = '/tmp/hinzmann/patTuple2.root'
process.out.outputCommands = [
         #------ MC block ----------------
         'keep *_addPileupInfo_*_*',  
         'keep recoGenJets_ak5GenJets*_*_*',
         'keep *_genParticles_*_*',
         #--------------------------------
         'keep *_kt6PFJets_rho_PAT',
         'keep *_kt6PFJetsCHS_rho_PAT',
         'keep *_kt6PFJetsISO_rho_PAT',## needed for the QG likelihood
         #'keep *_selectedPatJets__*',
         #'keep *_selectedPatJetsCHS__*',
         'keep *_jetExtender*_*_*', 
         'keep *_ca08PrunedPFJets_*_*',
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

process.multiJetFilter = cms.EDFilter('PatMultijetFilter',
    jets     = cms.InputTag('selectedPatJets'),
    minNjets = cms.int32(3),
    minPt    = cms.double(20)
)

process.maxEvents.input = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source.fileNames = [
#'/store/mc/Summer11/QCD_TuneZ2_HT-1000_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/1AB5A492-C4C5-E011-BCD9-90E6BA19A203.root',
'file:///tmp/hinzmann/009BA3E5-1952-E111-BEBE-0015178C4BB0.root',
]

#create pruned CA08 jets
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.ca08PrunedPFJets = cms.EDProducer(
                              "SubJetProducer",
                              PFJetParameters.clone(
                                  src = process.pfJetsCHS.src,
                                  doAreaFastjet = cms.bool(True),
                                  doRhoFastjet = cms.bool(False),
                                  ),
                              AnomalousCellParameters,
                              SubJetParameters,
                              jetAlgorithm = cms.string("CambridgeAachen"),
                              rParam = cms.double(0.8),
                              jetCollInstanceName=cms.string("subjets"),
                             )
process.ca08PrunedPFJets.nSubjets = cms.int32(2)

#add pruned CA08 subjets
addJetCollection(process, 
                 cms.InputTag('ca08PrunedPFJets:subjets'),
                 'CA8', 'Subjets',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel= ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = False,
                 doJetID      = False,
                )
process.patJetCorrFactorsCA8Subjets.rho = cms.InputTag("kt6PFJets","rho")
process.selectedPatJetsCA8Subjets.cut = "pt > 10 && abs(eta) < 4.7"

# extend ca08 Subjets
process.jetExtenderCA8Subjets = cms.EDProducer("JetExtendedProducer",
    jets    = cms.InputTag('selectedPatJetsCA8Subjets'),
    result  = cms.string('extendedPatJetsCA8Subjets'),
    payload = cms.string('AK5PFchs') 
)

# subjet multiplicity filter
process.multiSubjetFilter = cms.EDFilter('PatMultijetFilter',
    jets     = cms.InputTag('selectedPatJetsCA8Subjets'),
    minNjets = cms.int32(2),
    minPt    = cms.double(20)
)

process.options.wantSummary = False

process.p = cms.Path(
   #----- produce the HBHE noise flag --------------------------
   process.HBHENoiseFilterResultProducer +
   #----- re-cluster ak5PFJets after activating the jet area ---
   process.ak5PFJets +
   #----- re-cluster kt6PFJets after activating rho ------------
   process.kt6PFJets +
   #----- re-cluster kt6PFJets after activating rho for ISO ----
   process.kt6PFJetsISO +
   #----- create the collection of good PV ---------------------
   process.goodOfflinePrimaryVertices +
   #----- run the PF2PAT sequence: doing CHS -------------------
   process.patPF2PATseq +
   #----- run the pruned CA8 sequence -------------------------
   process.ca08PrunedPFJets +
   #----- run the default PAT sequence -------------------------
   process.patDefaultSequence +
   #----- extend the PAT jets with additional variables --------
   process.jetExtender +
   #----- extend the CHS PAT jets with additional variables ----
   process.jetExtenderCHS +
   #----- extend the ca8 subjects with additional variables ----
   process.jetExtenderCA8Subjets +
   #----- create a collection of tracks out of jets ------------
   process.outTracks +
   #----- reconstruct track jets from the soft tracks ----------
   process.ak5SoftTrackJets +
   process.multiJetFilter +
   process.multiSubjetFilter
)

