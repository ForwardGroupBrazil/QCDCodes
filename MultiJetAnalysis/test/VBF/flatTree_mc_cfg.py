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
        input = cms.untracked.int32(10000)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        "root://eoscms//eos/cms/store/cmst3/user/kkousour/VBF_HToBB_M-125_8TeV-powheg-pythia6/PAT-V3/95e858de58ce9ba2c8dcbcc034ee8f47/patTuple_9_1_gzC.root"
        )
)

############# QuarkGluon tagger ###########################
process.qglAK5PF = cms.EDProducer('QuarkGluonTagger',
    jets     = cms.InputTag('jetExtender','extendedPatJets'),
    rho      = cms.InputTag('kt6PFJetsISO','rho'),
    jec      = cms.string('ak5PFL1FastL2L3'),
    isPatJet = cms.bool(True)
)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('PatVBFTree',
    jets        = cms.InputTag('jetExtender','extendedPatJets'),
    met         = cms.InputTag('pfMet'),
    rho         = cms.InputTag('kt6PFJets','rho'),
    rhoQGL      = cms.InputTag('kt6PFJetsISO','rho'),
    dEtaMin     = cms.double(2.5),
    ptMin       = cms.vdouble(80,60,40,30),
    ## trigger ###################################
    triggerAlias     = cms.vstring('PF','Calo','PFJet80'),
    triggerSelection = cms.vstring(
      'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v* OR HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v* OR HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v* OR HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v*',
      'HLT_QuadJet75_55_35_20_BTagIP_VBF_v* OR HLT_QuadJet75_55_38_20_BTagIP_VBF_v*',
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
    pu          = cms.untracked.string('addPileupInfo'),
    genjets     = cms.untracked.InputTag('ak5GenJets'),
    btagger     = cms.string('combinedSecondaryVertexBJetTags')
)

process.p = cms.Path(process.qglAK5PF * process.Hbb)

