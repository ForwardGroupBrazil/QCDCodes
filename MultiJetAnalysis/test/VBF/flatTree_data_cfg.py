import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_52_V9::All'

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
'root://eoscms//eos/cms/store/cmst3/user/kkousour/Jet/Run2012A-PromptV1-PAT-V3/6c959e3053513b6fdbecab8620fff14d/patTuple_183_1_Dbm.root',
'root://eoscms//eos/cms/store/cmst3/user/kkousour/Jet/Run2012A-PromptV1-PAT-V3/6c959e3053513b6fdbecab8620fff14d/patTuple_184_1_uwy.root'
        )
)

############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring(
       'HLT_QuadJet75_55_35_20_BTagIP_VBF_v*', 
       'HLT_QuadJet75_55_38_20_BTagIP_VBF_v*',
       'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v*',
       'HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v*',
       'HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v*',
       'HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v*',
       'HLT_DiPFJetAve40_v*',  
       'HLT_DiPFJetAve80_v*', 
       'HLT_L1DoubleJet36Central_v*',  
       'HLT_PFJet40_v*',
       'HLT_PFJet80_v*'
    ),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.GluonTag = cms.EDProducer('GluonTagLikelihood',
    jets = cms.InputTag('jetExtender','extendedPatJets'),
    rho  = cms.InputTag('kt6PFJets','rho')
)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('PatVBFTree',
    jets             = cms.InputTag('jetExtender','extendedPatJets'),
    met              = cms.InputTag('pfMet'),
    rho              = cms.InputTag('kt6PFJets','rho'),
    rhoQGL           = cms.InputTag('kt6PFJetsISO','rho'),
    mbbMin           = cms.double(0.0),
    dEtaMin          = cms.double(0.0),
    ## trigger ###################################
    triggerAlias     = cms.vstring('PF','Calo','DiPFJetAve40','DiPFJetAve80','L1DoubleJet36Central','PFJet40','PFJet80'),
    triggerSelection = cms.vstring(
      'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v* OR HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v* OR HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v* OR HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v*',
      'HLT_QuadJet75_55_35_20_BTagIP_VBF_v* OR HLT_QuadJet75_55_38_20_BTagIP_VBF_v*',
      'HLT_DiPFJetAve40_v*',  
      'HLT_DiPFJetAve80_v*', 
      'HLT_L1DoubleJet36Central_v*',  
      'HLT_PFJet40_v*',
      'HLT_PFJet80_v*'
    ),
    triggerConfiguration = cms.PSet(
      hltResults            = cms.InputTag('TriggerResults','','HLT'),
      l1tResults            = cms.InputTag(''),
      daqPartitions         = cms.uint32(1),
      l1tIgnoreMask         = cms.bool(False),
      l1techIgnorePrescales = cms.bool(False),
      throw                 = cms.bool(False)
    ),
    btagger          = cms.string('combinedSecondaryVertexBJetTags'),
    qglFile          = cms.string('./QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root')
)

process.p = cms.Path(process.hltFilter * process.Hbb)

