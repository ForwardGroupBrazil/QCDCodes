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
        'file://./patTuple_signal.root'
        )
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('PatVBFTree',
    jets    = cms.InputTag('jetExtender','extendedPatJets'),
    met     = cms.InputTag('pfMet'),
    rho     = cms.InputTag('kt6PFJets','rho'),
    rhoQGL  = cms.InputTag('kt6PFJetsISO','rho'),
    mbbMin  = cms.double(50.0),
    dEtaMin = cms.double(1.0),
    pu      = cms.untracked.string('addPileupInfo'),
    btagger = cms.string('combinedSecondaryVertexBJetTags'),
    qglFile = cms.string('./QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root')
)

process.p = cms.Path(process.Hbb)

