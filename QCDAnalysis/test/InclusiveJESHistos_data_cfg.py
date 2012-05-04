import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('InclusiveJESHistos_data.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.inclusive    = cms.EDAnalyzer('InclusiveJESHistos',
    filenames = cms.vstring(
    'rfio:/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/Apr24th/Jet_Run2011A_May10_ProcessedTree_data.root',
    'rfio:/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/Apr24th/Jet_Run2011A_PromptV4_ProcessedTree_data.root',
    'rfio:/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/Apr24th/Jet_Run2011A_Aug05_ProcessedTree_data.root',
    'rfio:/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/Apr24th/Jet_Run2011A_PromptV6_ProcessedTree_data.root',
    'rfio:/castor/cern.ch/user/k/kkousour/7TeV/data/Jets/Inclusive/Apr24th/Jet_Run2011B_PromptV1_ProcessedTree_data.root'   
    ),
    treename         = cms.string('ProcessedTree'),
    dirname          = cms.string('ak7'),
    yBnd             = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    ptBnd            = cms.vdouble(40,53,67,81,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1684,1784,1890,2000,2116,2238,2366,2500,2640,2787,2941,3103,3273,3500),
    triggers         = cms.vstring('HLT_Jet60_v','HLT_Jet110_v','HLT_Jet190_v','HLT_Jet240_v','HLT_Jet370_v'),
    minPt            = cms.vdouble(114,196,300,362,507,3500), 
    maxMETovSumET    = cms.double(0.3),
    jecUncSrcNames  = cms.vstring('Absolute','HighPtExtra','SinglePion','Flavor','Time',
                                  'RelativeJEREC1','RelativeJEREC2','RelativeJERHF',
                                  'RelativeStatEC2','RelativeStatHF','RelativeFSR',
                                  'PileUpDataMC','PileUpOOT','PileUpPt','PileUpBias','PileUpJetRate',
                                  'SubTotalPileUp','SubTotalRelative','SubTotalPt','SubTotalDataMC','Total'),
    nEvents          = cms.int32(-1)
)

process.p = cms.Path(process.inclusive)

