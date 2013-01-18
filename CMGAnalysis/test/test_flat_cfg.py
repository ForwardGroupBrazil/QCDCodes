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
        input = cms.untracked.int32(-1)
        )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('file://///afs/cern.ch/work/k/kkousour/private/CMGTools/CMSSW_5_3_7/src/KKousour/CMGAnalysis/test/cmgTuple.root'
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
    jets     = cms.InputTag('cmgPFJetSel'),
    rho      = cms.InputTag('kt6PFJetsForIso','rho'),
    jec      = cms.string('ak5PFL1FastL2L3'),
    isCmgJet = cms.untracked.bool(True)
)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('VbfHbbFlatTreeProducer',
    jets             = cms.InputTag('cmgPFJetSel'),
    genjets          = cms.untracked.InputTag('genJetSel'), 
    genparticles     = cms.untracked.InputTag('genParticlesPruned'),
    met              = cms.InputTag('nopuMet'),
    rho              = cms.InputTag('kt6PFJets','rho'),
    rhoQGL           = cms.InputTag('kt6PFJetsForIso','rho'),
    shiftJES         = cms.double(0.0),
    dEtaMin          = cms.double(2.),
    ptMin            = cms.vdouble(60,40,30,20),
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
    pu               = cms.untracked.string('addPileupInfo')
)

process.HbbJesUp = process.Hbb.clone(shiftJES = 1.0);
process.HbbJesDo = process.Hbb.clone(shiftJES = -1.0);

process.HbbCHS = process.Hbb.clone(jets = cms.InputTag('cmgPFJetSelCHS'))
process.HbbCHSJesUp = process.HbbCHS.clone(shiftJES = 1.0);
process.HbbCHSJesDo = process.HbbCHS.clone(shiftJES = -1.0);

process.p = cms.Path(process.qglAK5PF * process.Hbb * process.HbbJesUp * process.HbbJesDo * process.HbbCHS * process.HbbCHSJesUp * process.HbbCHSJesDo)
