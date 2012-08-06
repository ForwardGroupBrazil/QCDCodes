import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        'root://eoscms//eos/cms/store/cmst3/user/kkousour/VBF_HToBB_M-120_TuneZ2star_8TeV-pythia6/PAT-V3/1e7170a6a5fe2be3759c57372f107e48/patTuple_8_1_wdT.root'
        )
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('PatVBFTree',
    jets        = cms.InputTag('jetExtender','extendedPatJets'),
    met         = cms.InputTag('pfMet'),
    rho         = cms.InputTag('kt6PFJets','rho'),
    rhoQGL      = cms.InputTag('kt6PFJetsISO','rho'),
    gluonJetMva = cms.InputTag('GluonTag'),
    mbbMin      = cms.double(0.0),
    dEtaMin     = cms.double(0.0),
    pu          = cms.untracked.string('addPileupInfo'),
    genjets = cms.untracked.InputTag('ak5GenJets'),
    btagger     = cms.string('combinedSecondaryVertexBJetTags'),
    qglFile     = cms.string('./QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root')
)

process.GluonTag = cms.EDProducer('GluonTagLikelihood',
    jets = cms.InputTag('jetExtender','extendedPatJets'),
    rho  = cms.InputTag('kt6PFJets','rho')
)

process.p = cms.Path(process.GluonTag * process.Hbb)

