import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('InclusiveHistos_data_44x.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.inclusive    = cms.EDAnalyzer('InclusiveHistos',
    filenames = cms.vstring(
    '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Mar3rd/Run2011A_ProcessedTree_data.root',
    '/uscms_data/d2/kkousour/7TeV/2011/Jets/data/Mar3rd/Run2011B_ProcessedTree_data.root' 
    ),
    treename         = cms.string('ProcessedTree'),
    dirname          = cms.string('ak7'),
    logname          = cms.string('InclusiveLog_44x.txt'),
    yBnd             = cms.vdouble(0.0,0.5,1.0,1.5,2.0,2.5,3.0),
    ptBnd            = cms.vdouble(40,53,67,81,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1684,1784,1890,2000,2116,2238,2366,2500,2640,2787,2941,3103,3273,3500),
    triggers         = cms.vstring('HLT_Jet30_v','HLT_Jet60_v','HLT_Jet110_v','HLT_Jet190_v','HLT_Jet240_v','HLT_Jet300_v','HLT_Jet370_v'),
    minPt            = cms.vdouble(53,114,196,300,362,395,507), 
    isMC             = cms.bool(False),    
    maxBetaStar      = cms.double(1.0),
    maxMETovSumET    = cms.double(0.3),
    applyHBEHfilter  = cms.bool(False),
    nEvents          = cms.int32(-1)
)

process.p = cms.Path(process.inclusive)

