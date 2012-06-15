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
        'file://./patTuple.root'
        )
)

############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    #HLTPaths           = cms.vstring('HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v*'),
    HLTPaths           = cms.vstring('HLT_QuadJet75_55_35_20_BTagIP_VBF_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.GluonTag = cms.EDProducer('GluonTagLikelihood',
    jets = cms.InputTag('jetExtender','extendedPatJets'),
    rho  = cms.InputTag('kt6PFJets','rho')
)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('PatVBFTree',
    jets    = cms.InputTag('jetExtender','extendedPatJets'),
    met     = cms.InputTag('pfMet'),
    rho     = cms.InputTag('kt6PFJets','rho'),
    rhoQGL  = cms.InputTag('kt6PFJetsISO','rho'),
    puJetMvaFull = cms.InputTag('puJetMva','fullDiscriminant'),
    puJetMvaSimple = cms.InputTag('puJetMva','simpleDiscriminant'),
    puJetIdCutBased = cms.InputTag('puJetMva','cutbasedId'),
    gluonJetMva = cms.InputTag('GluonTag'),
    mbbMin  = cms.double(0.0),
    dEtaMin = cms.double(0.0),
    btagger = cms.string('combinedSecondaryVertexBJetTags'),
    qglFile = cms.string('./QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root')
)

process.p = cms.Path(process.GluonTag * process.Hbb)

