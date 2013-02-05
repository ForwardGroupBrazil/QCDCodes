import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load('FWCore.MessageService.MessageLogger_cfi')

from CMGTools.Production.datasetToSource import *
from CMGTools.Common.Tools.applyJSON_cff import *

process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

process.source = datasetToSource('cmgtools','/JetMon/Run2012B-13Jul2012-v1/AOD/PAT_CMG_V5_12_0','cmgTuple_.*.root')

applyJSON(process,'Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt')

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
       'HLT_PFJet40_v*',
       'HLT_PFJet80_v*'
    ),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('VbfHbbFlatTreeProducer',
    jets             = cms.InputTag('cmgPFJetSel'),
    met              = cms.InputTag('cmgPFMETRaw'),
    rho              = cms.InputTag('kt6PFJets','rho'),
    shiftJES         = cms.double(0.0),
    dEtaMin          = cms.double(2.),
    ptMin            = cms.vdouble(60,40,30,20),
    ## trigger ###################################
    triggerAlias     = cms.vstring('PF','Calo','PF1','PF2','PF3','PF4','Calo1','Calo2','PFJet40','PFJet80','PFnoBtag','CaloNoBtag'),
    triggerSelection = cms.vstring(
      'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v* OR HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v* OR HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v* OR HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v*',
      'HLT_QuadJet75_55_35_20_BTagIP_VBF_v* OR HLT_QuadJet75_55_38_20_BTagIP_VBF_v*',
      'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v*',
      'HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v*',
      'HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v*',
      'HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v*',
      'HLT_QuadJet75_55_35_20_BTagIP_VBF_v*',
      'HLT_QuadJet75_55_38_20_BTagIP_VBF_v*',  
      'HLT_PFJet40_v*',
      'HLT_PFJet80_v*',
      'HLT_QuadPFJet78_61_44_31_VBF_v*',
      'HLT_Quadjet75_55_35_20_VBF_v*'
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

process.p = cms.Path(process.hltFilter * process.Hbb * process.HbbCHS)
