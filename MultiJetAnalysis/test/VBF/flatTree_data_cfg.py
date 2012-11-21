import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_53_V13::All'

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1000)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
'root://eoscms//eos/cms/store/cmst3/user/kkousour/BJetPlusX/Run2012B-13Jul2012-PAT-V3/5ecefa428a13706983a82c8246282960/patTuple_99_2_inp.root'
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

############# QuarkGluon tagger ###########################
process.qglAK5PF = cms.EDProducer('QuarkGluonTagger',
    jets     = cms.InputTag('jetExtender','extendedPatJets'),
    rho      = cms.InputTag('kt6PFJetsISO','rho'),
    jec      = cms.string('ak5PFL1FastL2L3Residual'),
    isPatJet = cms.bool(True)
)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('PatVBFTree',
    jets             = cms.InputTag('jetExtender','extendedPatJets'),
    met              = cms.InputTag('pfMet'),
    rho              = cms.InputTag('kt6PFJets','rho'),
    rhoQGL           = cms.InputTag('kt6PFJetsISO','rho'),
    dEtaMin          = cms.double(2.5),
    ptMin            = cms.vdouble(80,60,40,30),
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
    btagger          = cms.string('combinedSecondaryVertexBJetTags')
)

process.p = cms.Path(process.hltFilter * process.qglAK5PF * process.Hbb)

